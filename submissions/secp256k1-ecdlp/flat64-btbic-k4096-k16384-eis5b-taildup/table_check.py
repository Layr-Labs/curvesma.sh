"""Finite table-builder value/count checks, only run in verifier Docker.

The independent formulas use descriptor integers and modular recurrence.
Expanded instruction counts include call, row, tail and setup work. These small
fixtures do not claim to allocate or execute production-size tables.
"""
import json, pathlib, random, sys, time
from collections import defaultdict


def check(root, cert):
    started=time.monotonic()
    root = pathlib.Path(root)
    sys.path.insert(0, str(root/'kernels'))
    import a64
    field = a64.a64_load((root/'kernels/secp_fe.s').read_text())
    route = (root/'kernels/route.s').read_text()
    route = a64.a64_load(route[route.index('.globl tbld\n'):route.index('.globl tinit\n')])
    p, mask = cert.p, (1 << 64)-1
    PB, SP, AB, DST, SENT = 0x600000, 0x500000, 0x300000, 0x700000, 0x51515151
    rng = random.Random(20261123)

    def limbs(v, n=4):
        return [(v >> (64*i)) & mask for i in range(n)]

    def put(mem, addr, v, n=4):
        for i, w in enumerate(limbs(v, n)): mem[addr+8*i] = w

    def get(mem, addr, n=4):
        return sum(mem[addr+8*i] << (64*i) for i in range(n))

    def descriptor(c):
        # Explicit ABI: four field residues, the 192-bit modular inverse,
        # lattice vectors, offsets, bounds, then all candidate translations.
        result = limbs(c['W'][1])+limbs(p-c['U'][1])+limbs(c['kS'])+limbs(c['kT'])
        result += limbs(pow(p, -1, 1 << 192), 3)
        result += limbs(c['U'][0], 3)+limbs(c['W'][0], 3)
        result += limbs(c['U'][1], 2)+limbs(c['W'][1], 2)
        result += limbs(c['o'][0], 3)+limbs(c['o'][1], 2)
        result += limbs(c['C'], 3)+limbs(c['X2'], 2)+[len(c['c'])]
        for a, b in c['c']: result += limbs(a, 3)+limbs(b, 2)
        return result

    def reduction(c, v):
        # Direct integer definition, separate from certificate.eis_red.
        U, W, o = c['U'], c['W'], c['o']
        s = (W[1]*v+c['kS']) % p
        t = (-U[1]*v+c['kT']) % p
        a = o[0]+(s*U[0]+t*W[0])//p
        b = o[1]+(s*U[1]+t*W[1])//p
        for da, db in c['c']:
            if abs(a+da) <= c['C'] and 0 <= b+db < c['X2']:
                return (a+da) % (1 << 192), b+db
        return a % (1 << 192), b

    def run(prog, entry, mem, regs):
        rr = dict(sp=SP, x30=SENT); rr.update(regs)
        count = a64.a64_run(*prog, entry, rr, mem, out_regs=rr, limit=200000)
        assert rr['sp'] == SP and rr['x30'] == SENT
        return count

    code, labels = field
    kstart, kend = labels['tbldk'][0], labels['tbldf'][0]
    rowpc = next(i for i in range(kstart, kend-1)
                 if code[i] == ('ldr', ['x0','[sp,#2032]'])
                 and code[i+1] == ('add', ['x1','x29','#128']))
    # The static row consists of its linear body and five additional passes
    # through the six-way update loop. Calls are counted separately below.
    rowtail = next(i for i in range(rowpc, kend)
                   if code[i] == ('add', ['sp','sp','#2016']))
    inner = next(i for i in range(rowpc, rowtail)
                 if code[i] == ('ldp', ['x10','x11','[x1,#0]']))
    inner_end = next(i for i in range(inner, rowtail)
                     if code[i] == ('b.ne', ['2b']))+1
    static_row = rowtail-rowpc+5*(inner_end-inner)
    assert inner_end-inner == 74 and static_row == 487
    result = {'scope':'Actual expanded assembly on finite full-field rows; no production table allocation.',
              'static_tbldk_row_without_eisr3':static_row, 'tbldk':[], 'tbldf':[], 'tbld':[]}
    controls = {}

    class Trace(defaultdict):
        def __init__(self):
            super().__init__(int); self.calls=0; self.in_row=0; self.krows=0
        def instruction(self, pc, op, args, regs, ev):
            self.calls += (pc == labels['eisr3'][0])
            self.krows += (pc == rowpc)
            self.in_row += (rowpc <= pc < rowtail)

    starts_all = [[0,1,p-1,p-2,p//2,1 << 255], [rng.randrange(p) for _ in range(6)]]
    steps_all = [[p-1,p-1,1,p-2,p//2,1 << 255], [rng.randrange(p) for _ in range(6)]]
    for cl, c in cert.EIS.items():
        series=[]
        eis = cert.EISR3[cl]
        assert eis == 43
        for fixture, (starts, steps) in enumerate(zip(starts_all, steps_all)):
            counts=[]
            for rows in (1,2,3):
                mem=Trace()
                for i, w in enumerate(descriptor(c)): mem[PB+8*i]=w
                # ksel: the 24 words below the descriptor (band thresholds, column thresholds, candidate table)
                for i, w in enumerate(cert.ksel_ext(cl)): mem[PB-192+8*i]=w
                for j in range(6):
                    put(mem, AB+64*j, steps[j]); put(mem, AB+64*j+32, starts[j])
                count=run(field, 'tbldk', mem, dict(x0=DST,x1=AB,x2=rows,x3=PB))
                for row in range(rows):
                    addr=DST+256*row
                    assert mem[addr+120] == addr and mem[addr+248] == addr+128
                    for j, off in enumerate((0,40,80,128,168,208)):
                        want=reduction(c,(starts[j]+row*steps[j]) % p)
                        assert (get(mem,addr+off,3),get(mem,addr+off+24,2)) == want
                assert mem.krows == rows and mem.calls == 6*rows
                assert mem.in_row == 487*rows
                assert count == 4268+rows*(487+6*eis), (cl,rows,count)
                legacy=4268+rows*(486+6*eis)
                assert count-legacy == rows
                assert cert.tbldk_count(rows, cl) >= count
                counts.append(count)
            slope=counts[1]-counts[0]
            assert counts[2]-counts[1] == slope and counts[0]-slope == 4268
            series.append(counts)
        result['tbldk'].append(dict(class_name=cl,counts_by_fixture=series,
            intercept=4268,slope=slope,eisr3=eis,legacy_undercharge_per_row=1))
        # EISR3 at both corners of every ksel cell: the first valid candidate, in 43 instructions.
        corners=0
        for a, b in cert.KSEL[cl]['pts']:
            mem=defaultdict(int)
            for i, w in enumerate(descriptor(c)): mem[PB+8*i]=w
            for i, w in enumerate(cert.ksel_ext(cl)): mem[PB-192+8*i]=w
            put(mem, AB, a % (1 << 192), 3); put(mem, AB+24, b, 2)
            count=run(field, 'eisr3', mem, dict(x0=DST, x1=AB, x2=PB))
            w=cert.eis_first(c, a, b)
            assert count == 43 and (get(mem,DST,3), get(mem,DST+24,2)) == (w[0] % (1 << 192), w[1]), (cl, a, b)
            corners+=1
        result['tbldk'][-1]['eisr3_cell_corners']=corners
    # The pricing assertion really rejects the old expression on executed rows.
    observed=result['tbldk'][0]['counts_by_fixture'][0][0]
    oldprice=4268+486+6*result['tbldk'][0]['eisr3']
    rejected=False
    try: assert oldprice >= observed
    except AssertionError: rejected=True
    assert rejected
    controls['old_486_price_rejected']=rejected

    for variants, stride in ((1, 64), (2, 64), (4, 128)):
        for rows in (1, 2, 3):
            mem=defaultdict(int); starts=starts_all[0][:variants]; steps=steps_all[0][:variants]
            for j in range(variants):
                put(mem,AB+64*j,steps[j]); put(mem,AB+64*j+32,starts[j])
            count=run(field,'tbldf',mem,dict(x0=DST,x1=AB,x2=rows,x3=variants,x26=stride))
            for row in range(rows):
                for j in range(variants):
                    v=(starts[j]+row*steps[j]) % p
                    want=v+((v*cert.C0+(1 << 255)) >> 256)
                    assert get(mem,DST+stride*row+32*j) == want
            for j in range(variants): assert get(mem,AB+64*j+32) == (starts[j]+rows*steps[j]) % p
            assert count == 4+rows*(6+46*variants) == cert.tbldf_count(rows,variants)
            result['tbldf'].append([variants,rows,count])
    # BT2 window 4 is stored reversed at +32 with stride -64; same tbldf count.
    for rows in (1,2,3):
        mem=defaultdict(int); starts=starts_all[0][:1]; steps=steps_all[0][:1]; top=DST+64*3+32
        put(mem,AB,steps[0]); put(mem,AB+32,starts[0])
        count=run(field,'tbldf',mem,dict(x0=top,x1=AB,x2=rows,x3=1,x26=(1<<64)-64))
        for row in range(rows):
            v=(starts[0]+row*steps[0]) % p
            want=v+((v*cert.C0+(1 << 255)) >> 256)
            assert get(mem,top-64*row) == want
        assert all(mem[top-64*rows+8*w] == 0 for w in range(4))
        assert count == 4+rows*(6+46) == cert.tbldf_count(rows,1)
        result['tbldf'].append(['reversed',rows,count])

    # Normal/alternate per-variant row selection is tested in both mask bits,
    # including unequal lengths and sentinels beyond each logical table.
    for normal, alternate, bitmap in [(n,n,0) for n in (1,2,3)]+[(1,3,m) for m in range(4)]:
        mem=defaultdict(int); variants=2
        starts=[p-1,1];steps=[p-2,p-1];sentinel=0xDEADBEEFCAFEF00D
        for j in range(variants):
            put(mem,AB+64*j,steps[j]);put(mem,AB+64*j+32,starts[j])
            for row in range(4):
                for w in range(4):mem[DST+64*row+32*j+8*w]=sentinel
        count=run(route,'tbld',mem,dict(x0=DST,x1=AB,x2=normal,x3=alternate,x4=variants,x5=64,x6=bitmap))
        lengths=[alternate if bitmap>>j&1 else normal for j in range(variants)]
        for j,n in enumerate(lengths):
            for row in range(n): assert get(mem,DST+64*row+32*j) == (starts[j]+row*steps[j]) % p
            for row in range(n,4):
                assert all(mem[DST+64*row+32*j+8*w] == sentinel for w in range(4))
        assert count == 4+11*variants+20*sum(lengths)
        result['tbld'].append([normal,alternate,bitmap,count])
    # P1 window 4 is stored reversed: the builder starts at the top row with stride -64; same count.
    for rows in (1,2,3):
        mem=defaultdict(int); variants=2; starts=[p-1,1]; steps=[p-2,p-1]; top=DST+64*3
        for j in range(variants):
            put(mem,AB+64*j,steps[j]); put(mem,AB+64*j+32,starts[j])
        count=run(route,'tbld',mem,dict(x0=top,x1=AB,x2=rows,x3=rows,x4=variants,x5=(1<<64)-64,x6=0))
        for j in range(variants):
            for row in range(rows): assert get(mem,top-64*row+32*j) == (starts[j]+row*steps[j]) % p
        assert all(mem[top-64*rows+8*w] == 0 for w in range(8))
        assert count == 4+11*variants+20*variants*rows
        result['tbld'].append(['reversed',rows,count])
    result['negative_controls']=controls
    result['ok']=True
    result['seconds']=round(time.monotonic()-started, 3)
    return result


if __name__ == '__main__':
    root=next(p for p in (pathlib.Path(__file__).resolve().parent,pathlib.Path('/snapshot/files'),pathlib.Path.cwd()) if (p/'kernels/a64.py').is_file())
    sys.path.insert(0,str(root))
    import certificate
    print(json.dumps(check(root,certificate),sort_keys=True))
