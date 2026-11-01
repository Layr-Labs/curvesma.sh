"""Exact instruction-renaming and dead-scratch-store elimination witness.
Based on our published paired-inline-v1, not a production cost proof by itself.
"""
import re
from pathlib import Path
OUT={'k_mul':[9,10,11,12],'k_sqr':[7,8,9,10],'k_add':[11,12,13,14],'k_sub':[11,12,13,14]}

def generate(root):
    src=(root/'p256_fields.s').read_text();ks={}
    for k in OUT:
        ks[k]=src.split(k+':\n')[1].split(' ret')[0].strip().splitlines()
    assert all(not re.search(r'\b[wx]18\b',line) for body in ks.values() for line in body)
    assert all(not re.match(r'\s*(?!ldp\b)\w+\s+x[78]\s*,',line) for line in ks['k_mul'])
    fwd=(root/'pair_inline_reference.s').read_text().split('.globl pair_bwd')[0]
    old=' ldp x9,x10,[x21]\n ldp x11,x12,[x21,#16]\n orr x9,x9,x10\n orr x11,x11,x12\n orr x9,x9,x11'
    assert fwd.count(old)==1
    fwd=fwd.replace(old,' orr x9,x11,x12\n orr x10,x13,x14\n orr x9,x9,x10')
    # Preserve the four canonical difference limbs in x11..x14.
    # Rename multiplication temporaries away from those live registers.
    start=fwd.index(' // inline k_mul;')
    end=fwd.index(' add x23,x21,#32',start)
    block=fwd[start:end].splitlines(); rewritten=[]; hi=False; removed=0
    for line in block:
        if line.strip()=='ldp x7, x8, [x21]':
            removed+=1;continue
        if line.strip()=='ldp x7, x8, [x21, #16]':
            removed+=1;hi=True;continue
        mapping={11:24,12:25,13:26,14:27,7:(13 if hi else 11),8:(14 if hi else 12)}
        line=re.sub(r'\b([xw])(\d+)\b',lambda m:m[1]+str(mapping.get(int(m[2]),int(m[2]))),line)
        rewritten.append(line)
    assert removed==2
    fwd=fwd[:start]+'\n'.join(rewritten)+'\n'+fwd[end:]
    # Post-index first load advances Y once; second load uses adjusted offset.
    fwd=fwd.replace(' ldp x3, x4, [x20]', ' ldp x3, x4, [x20], #64').replace(' ldp x5, x6, [x20, #16]', ' ldp x5, x6, [x20, #-48]').replace(' add x20,x20,#64\n','')
    specs=[
      ('k_mul',('x23',32),('x23',0),('x21',-32),False,False),
      ('k_mul',('x23',0),('x23',0),('x21',0),True,False),
      ('k_add',('x23',192),('x20',0),('x22',0),True,False),
      ('k_sub',('x23',96),('x20',32),('x22',32),False,False),
      ('k_add',('x23',64),('x20',32),('x22',32),False,'save_plus'),
      ('k_mul',('x23',128),('x23',96),('x23',32),False,True),
      ('k_sqr',('x23',128),('x23',128),None,False,True),
      ('k_sub',('x21',0),('x23',128),('x23',192),True,True),
      ('k_mul',('x23',160),('x23',64),('x23',32),False,'use_plus'),
      ('k_sqr',('x23',160),('x23',160),None,False,True),
      ('k_sub',('x21',32),('x23',160),('x23',192),True,True)]
    lines=['.globl pair_bwd','pair_bwd:',' str x30,[sp,#-16]!',' mov x19,x0',' mov x20,x1',' mov x21,x2',' mov x22,x3',' mov x23,x4',' cbz x19,pb_done','pb_loop:']
    previous=None;previous_location=None;eliminated=0;steps=0;cached_plus=None
    for index,(k,dst,a,b,store,chain) in enumerate(specs):
        mapping={i:i for i in range(19)}
        saved_previous=(previous,previous_location)
        if chain=='use_plus':
            previous=cached_plus;previous_location=('x23',64)
        if chain=='save_plus':
            assert previous==[11,12,13,14] and previous_location==('x23',96)
            mapping.update({11:28,12:29,13:30,14:18,15:0,0:15,18:14})
            assert all(not re.match(r'\s*\w+ x(?:3|4|5|6|7|8|9|10),',x) for x in ks['k_sub'][4:] if not x.strip().startswith('stp'))
        elif chain:
            assert previous is not None and previous_location==a
            mapping=dict(zip([3,4,5,6],previous));unused=set(range(19))-set(previous)
            for i in range(19):
                if i in mapping:continue
                if i in unused:mapping[i]=i;unused.remove(i)
            for i in range(19):
                if i not in mapping:mapping[i]=min(unused);unused.remove(mapping[i])
        if index==0:
            mapping.update({9:24,10:25,11:26,12:27})
        assert len(set(mapping.values()))==19
        loc={'x0':dst,'x1':a,'x2':b}
        right_high=False
        for raw in ks[k]:
            if index in (5,8) and re.search(r'^\s*ldp .*\[x2(?:,|\])',raw):
                right_high='#16' in raw
                eliminated+=1
                continue
            if chain=='save_plus' and re.search(r'^\s*ldp .*\[x[12](?:,|\])',raw):eliminated+=1;continue
            if chain and chain!='save_plus' and re.search(r'^\s*ldp .*\[x1(?:,|\])',raw):eliminated+=1;continue
            if not store and re.search(r'^\s*stp .*\[x0(?:,|\])',raw):eliminated+=1;continue
            # Replace pointers with placeholders before bijective GPR renaming.
            def address(m):
                reg,off=m.groups();p,n=loc[reg];n+=int(off or 0)
                return '[PTR'+p[1:]+(', #'+str(n) if n else '')+']'
            raw=re.sub(r'\[(x[012])(?:,\s*#(-?\d+))?\]',address,raw)
            current=dict(mapping)
            if index in (5,8):current.update({7:(26 if right_high else 24),8:(27 if right_high else 25)})
            raw=re.sub(r'\b([xw])(\d+)\b',lambda m:m[1]+str(current[int(m[2])]) if int(m[2])<19 else m[0],raw)
            raw=re.sub(r'PTR(\d+)',r'x\1',raw)
            lines.append(' '+raw.strip());steps+=1
        previous=[mapping[i] for i in OUT[k]] if not store else None
        previous_location=dst
        if chain=='save_plus':
            cached_plus=previous
            previous,previous_location=saved_previous
    lines+=[' sub x20,x20,#64',' sub x21,x21,#64',' subs x19,x19,#1',' b.ne pb_loop','pb_done:',' ldr x30,[sp],#16',' ret']
    assert eliminated==34
    return fwd+'\n'.join(lines)+'\n',{'removed_backward_memory_instructions':eliminated,'removed_forward_loads':4,'backward_body_instructions':steps+4,'forward_per_item':139,'backward_per_item':710,'calls_reordered':len(specs)}

def check(root):
    text,meta=generate(root)
    assert text==(root/'pair_batch.s').read_text()
    return meta
