"""Deterministic address-specialization witness; reference batch is attributed prior code."""
import re
from pathlib import Path

def check(root):
    text=(root/'p256_fields.s').read_text(); kernels={}
    for name in ['k_mul','k_sqr','k_add','k_sub']:
        body=text.split(name+':\n')[1].split(' ret')[0].strip().splitlines()
        assert all(':' not in line for line in body)
        assert all(not re.search(r'\bx(?:19|2[0-9])\b',line) for line in body)
        kernels[name]=body
    out=[];maps={};saved={'pair_fwd':0,'pair_bwd':0};section=None;calls=0
    for line in (root/'pair_batch_reference.s').read_text().splitlines():
        if line.strip() in ('pair_fwd:','pair_bwd:'):section=line.strip()[:-1]
        m=re.fullmatch(r' (add|sub) (x[012]),(x\d+),#(\d+)',line)
        if m:
            op,reg,ptr,n=m.groups();maps[reg]=(ptr,int(n)*(1 if op=='add' else -1));continue
        m=re.fullmatch(r' bl (k_mul|k_sqr|k_add|k_sub)',line)
        if m:
            name=m[1];needed={'x0','x1'} if name=='k_sqr' else {'x0','x1','x2'}
            assert set(maps)==needed
            for inst in kernels[name]:
                def replace(m):
                    reg,old=m.groups();ptr,n=maps[reg];n+=int(old or 0)
                    assert -512<=n<=504 and n%8==0
                    return '['+ptr+(', #'+str(n) if n else '')+']'
                out.append(re.sub(r'\[(x[012])(?:,\s*#(-?\d+))?\]',replace,inst))
            saved[section]+=len(needed)+2;calls+=1;maps={}
        else:
            assert not maps
            out.append(line)
    def normal(lines):return [re.sub(r'\s+','',x.split('//')[0]) for x in lines if x.split('//')[0].strip()]
    assert normal(out)==normal((root/'pair_batch.s').read_text().splitlines())
    assert saved=={'pair_fwd':10,'pair_bwd':53} and calls==13
    return {'address_specialization_matches':True,'field_calls_inlined':calls,'saved_instructions_per_forward_item':10,'saved_instructions_per_backward_item':53}
