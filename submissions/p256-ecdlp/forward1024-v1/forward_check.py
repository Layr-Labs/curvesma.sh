import json,re
from pathlib import Path
root=Path('/snapshot/files')
s=(root/'measurement/floop.h').read_text()
body=next(x for x in s.splitlines() if x.startswith('#define FLOOP_S'))
loop=next(x for x in s.splitlines() if x.startswith('#define FLOOP '))
copies=[int(x) for x in re.findall(r'FC\((\d+)\)',loop)]
assert copies==list(range(1024))
assert 'subs x26,x26,#1024\\nb.ne 1b' in loop
assert 'b Lf2_%=' in loop
# Each escaped newline accounts for at most one instruction: labels only
# overcount. No assembler alignment, variable directives or literal pools
# are introduced inside this macro.
bound=4*(body.count('\\n')*1024+8)
assert bound < (1<<20)-4
assert not any(x in body for x in ['.align','.p2align','.space','.inst','.word'])
assert not re.search(r'\bx26\b',body)
checks=[]
for B in [1024,2048,4096,8192,1<<16,1<<20,1<<28]:
    old=512-2+512*(B//512-1)
    new=1024-2+1024*(B//1024-1)
    assert old==new==B-2
    saved=2*(B//512-B//1024)
    assert saved==B//512
    if B<=8192:
        def walk(U):
            n=0; out=[]
            for k in range(B//U):
                for j in range(2 if k==0 else 0,U):
                    out.append(n+1); n+=1
            return out
        assert walk(512)==walk(1024)==list(range(1,B-1))
    checks.append({'batch':B,'lanes':new,'saved_instructions':saved})
print(json.dumps({'all_checks':True,'cases':checks,
                  'forward_block_byte_upper_bound':bound,
                  'scope':'Exact schedule and conservative source-level branch-range audit; not a compiled measurement'}))
