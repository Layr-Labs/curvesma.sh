"""Local research checks, not a production claim or a full pipeline test."""
import contextlib,io,json,random,runpy,sys
from pathlib import Path

base=Path(__file__).resolve().parent
base=base/'files' if (base/'files'/'certificate.py').is_file() else base
sys.path.insert(0,str(base/'kernels'))
import a64
program,labels=a64.a64_load((base/'kernels'/'scan_ranges.s').read_text())
inputs=json.load(sys.stdin)
rng=random.Random(72193+int(inputs.get("seed",1)))
cases=0
max_ratio=0
for buckets in (0,1,2,4,8,16):
    for trial in range(40):
        groups=[[rng.randrange(4) for _ in range(rng.randrange(18))] for _ in range(buckets)]
        for cap in (0,1,2,8,1000):
            pa,ends,out,records=0x10000,0x20000,0x30000,0x40000
            memory={};n=0;expected=[]
            for bucket,keys in enumerate(groups):
                addresses=[]
                for key in keys:
                    address=records+32*n;memory[pa+8*n]=address
                    memory[address]=key
                    memory[address+8]=(rng.randrange(1<<16)<<48)|bucket
                    addresses.append(address);n+=1
                expected.extend((addresses[i-1],addresses[i]) for i in range(1,len(keys)) if keys[i-1]==keys[i])
                memory[ends+8*bucket]=pa+8*n
            regs={}
            count=a64.a64_run(program,labels,'scan_ranges',dict(x0=buckets,x1=pa,x2=ends,x3=out,x4=out+16*cap,x30=0x5E5E5E5E5E5E5E5E),memory,limit=100000,out_regs=regs)
            got=[(memory[p],memory[p+8]) for p in range(out,regs['x3'],16)]
            assert got==expected[:cap],(buckets,trial,cap,got,expected[:cap])
            bound=(15*n+3)//4+40*buckets+4*(cap+1)+2
            assert count<=bound,(count,bound,n,buckets,cap)
            cases+=1;max_ratio=max(max_ratio,count/bound)
buffer=io.StringIO()
with contextlib.redirect_stdout(buffer):
    certificate=runpy.run_path(str(base/'certificate.py'),run_name='research_certificate')
print(json.dumps({'research_only':True,'scanner_cases':cases,'scanner_correct':True,'instruction_bound_verified_on_cases':True,'max_observed_fraction_of_bound':max_ratio,'certificate_failed_checks':certificate['FAILED'],'q1':certificate['Q1'],'work_log2':certificate['lg'](certificate['W']),'computed_claim':str(certificate['CLAIM']),'per_leaf_instructions':str(certificate['PER_LEAF']),'success_lower':float(certificate['T']['lo']),'certificate_log_tail':buffer.getvalue()[-1500:]}))
