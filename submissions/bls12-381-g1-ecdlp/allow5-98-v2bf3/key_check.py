"""Independent limb semantics and worst-path audit of the changed key block."""
import random
from pathlib import Path

P = int('1a0111ea397fe69a4b1ba7b6434bacd7'
        '64774b84f38512bf6730d2a0f6b0f624'
        '1eabfffeb153ffffb9feffffffffaaab', 16)
MASK = (1 << 64) - 1

def limbs(x):
    return [(x >> (64*i)) & MASK for i in range(6)]

def branch_key(y):
    a, p = limbs(y), limbs(P)
    cost = 1 # TBZ, both outcomes charged
    if a[0] & 1:
        borrow = 0
        for i in range(4):
            v = p[i] - a[i] - borrow
            borrow = int(v < 0)
            if i:
                a[i] = v & MASK
            cost += 1
    cost += 2 # STP and STR
    return a[1] | (a[2] << 64) | (a[3] << 128), cost

def check():
    rng = random.Random(20261031)
    cases = {0,1,2,P-1,P-2,P//2,P//2+1}
    for i in range(1,6):
        for base in [1 << (64*i), P % (1 << (64*i)), P-(1 << (64*i))]:
            cases.update(x for x in range(base-4,base+5) if 0 <= x < P)
    cases.update(rng.randrange(P) for _ in range(12000))
    path_counts = {3:0,7:0}
    for y in cases:
        key,cost=branch_key(y)
        folded=y if y%2==0 else P-y
        assert key == ((folded >> 64) & ((1 << 192)-1))
        assert key == branch_key((P-y)%P)[0]
        assert cost <= 7
        path_counts[cost]+=1
    expected=['tbz x11,#0,.Lparity_key_store','subs x19,x5,x11',
              'sbcs x12,x6,x12','sbcs x13,x7,x13','sbc x14,x8,x14',
              'stp x12,x13,[x2,#96]','str x14,[x2,#112]']
    root=Path(__file__).resolve().parent
    for file in ['measure/side_lanes.S','measure/side_lanes.S.in']:
        s=(root/file).read_text()
        s=s[s.index('\ttbz\tx11,#0,.Lparity_key_store'):]
        s=s[:s.index('\tcmp\tx18,x30')]
        instructions=[' '.join(x.split()) for x in s.splitlines()
                      if x.strip() and not x.strip().endswith(':')]
        assert instructions==expected,(file,instructions)
    return {'cases':len(cases),'path_counts':path_counts,'worst_instructions':7,
            'baseline_instructions':10,'all_equal':True,
            'scope':'limb semantics and instruction paths, not hardware timing'}
