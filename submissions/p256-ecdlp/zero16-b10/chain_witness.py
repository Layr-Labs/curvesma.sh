"""Exact instruction-renaming and dead-scratch-store elimination witness.
Based on our published paired-inline-v1, not a production cost proof by itself.
"""
import re
from pathlib import Path
BASE_OUT={'k_mul':[9,10,11,12],'k_sqr':[7,8,9,10],'k_add':[11,12,13,14],'k_sub':[11,12,13,14]}
def generate(root):
    OUT=dict(BASE_OUT)
    src=(root/'p256_fields.s').read_text();ks={}
    for k in OUT:
        ks[k]=src.split(k+':\n')[1].split(' ret')[0].strip().splitlines()
    ks['k_abs']=ks['k_sub'][:8]+[' b.hs abs_store',' subs x11,xzr,x11',' sbcs x12,xzr,x12',' sbcs x13,xzr,x13',' sbc x14,xzr,x14','abs_store:',' stp x11,x12,[x0]',' stp x13,x14,[x0,#16]']
    OUT['k_abs']=[11,12,13,14]
    fused=(root/'fused.s').read_text().split('k_fused:\n')[1].split(' ret')[0]
    fused=fused.split(' lsr x15,x10,#32')[0]
    ks['k_fused']=fused.strip().splitlines()
    OUT['k_fused']=[7,8,9,10]
    removed_constants={}
    for name,body in list(ks.items()):
        tail=False;out=[];removed=0
        for raw in body:
            if re.fullmatch(r'\s*mov w3,\s*#-1',raw):tail=True;removed+=1;continue
            if tail and re.fullmatch(r'\s*mov w4,\s*#-2',raw):removed+=1;continue
            if tail:
                raw=re.sub(r'\bx3\b','CONST_A',raw)
                raw=re.sub(r'\bx4\b','CONST_B',raw)
            out.append(raw)
        ks[name]=out;removed_constants[name]=removed
    for name,body in list(ks.items()):
        out=[]
        for raw in body:
            if 'CONST_B' in raw:
                op=raw.strip().split()[0]
                assert op in ['adc','adcs','sbc','sbcs'],raw
                raw=re.sub(r'\b'+op+r'\b',{'adc':'sbc','adcs':'sbcs','sbc':'adc','sbcs':'adcs'}[op],raw,count=1)
            out.append(raw)
        ks[name]=out
    mc='mov x17, #0xFFFFFFFF00000001'
    assert sum(x.strip()==mc for x in ks['k_mul'])==1
    ks['k_mul']=[re.sub(r'\bx17\b','CONST_B',x) for x in ks['k_mul'] if x.strip()!=mc]
    assert removed_constants=={'k_mul':2,'k_sqr':2,'k_add':2,'k_sub':0,'k_abs':0,'k_fused':4}
    assert all(not re.search(r'\b[wx]18\b',line) for body in ks.values() for line in body)
    assert all(not re.match(r'\s*(?!ldp\b)\w+\s+x(?:7|8|13|14)\s*,',line) for line in ks['k_mul'])
    out=[];locs={};serial=0
    ref=(root/'pair_batch_reference.s').read_text().split('.globl pair_bwd')[0]
    for line in ref.splitlines():
        m=re.fullmatch(r' (add|sub) (x[012]),(x\d+),#(\d+)',line)
        if m:
            op,reg,ptr,n=m.groups();locs[reg]=(ptr,int(n)*(1 if op=='add' else -1));continue
        m=re.fullmatch(r' bl (k_mul|k_sub)',line)
        if m:
            name=m[1];key='k_abs' if name=='k_sub' else name
            out.append(' // inline '+name+'; redundant/absolute research')
            for raw in ks[key]:
                def addr(m):
                    reg,off=m.groups();ptr,n=locs[reg];n+=int(off or 0)
                    return '['+ptr+(', #'+str(n) if n else '')+']'
                raw=re.sub(r'\[(x[012])(?:,\s*#(-?\d+))?\]',addr,raw)
                for label in ['k_mul_store','abs_store']:raw=raw.replace(label,label+'_fw'+str(serial))
                out.append(raw.strip() if raw.strip().endswith(':') else ' '+raw.strip())
            serial+=1;locs={}
        else:out.append(line)
    fwd='\n'.join(out)+'\n'
    old=' ldp x9,x10,[x21]\n ldp x11,x12,[x21,#16]\n orr x9,x9,x10\n orr x11,x11,x12\n orr x9,x9,x11\n cbz x9,pf_fail\n'
    assert fwd.count(old)==1
    fwd=fwd.replace(old,'')
    old_fail='pf_fail:\n mov x0,#0\n'
    assert fwd.count(old_fail)==1
    fwd=fwd.replace(old_fail,'')
    start=fwd.index(' // inline k_mul;')
    end=fwd.index(' add x23,x21,#32',start)
    block=fwd[start:end].splitlines(); rewritten=[]; hi=False; removed=0
    for line in block:
        if line.strip()=='ldp x7, x8, [x21]':
            removed+=1;continue
        if line.strip()=='ldp x13, x14, [x21, #16]':
            removed+=1;hi=True;continue
        mapping={11:24,12:25,13:13,14:14,15:26,16:27,17:1,7:11,8:12,9:15,10:16}
        line=re.sub(r'\b([xw])(\d+)\b',lambda m:m[1]+str(mapping.get(int(m[2]),int(m[2]))),line)
        rewritten.append(line.strip() if line.strip().endswith(':') else line)
    assert removed==2
    fwd=fwd[:start]+'\n'.join(rewritten)+'\n'+fwd[end:]
    fwd=fwd.replace(' ldp x3, x4, [x20]', ' ldp x3, x4, [x20], #64').replace(' ldp x5, x6, [x20, #16]', ' ldp x5, x6, [x20, #-48]').replace(' add x20,x20,#64\n','')
    # b10: hoist the loop-invariant X loads to the prologue (k_mul no longer
    # writes x9/x10, so x7-x10 survive the whole loop); read A from the
    # previous batch entry relative to x21; delay the dead-until-bwd d stores
    # past k_mul and post-increment x21 on the last one. Drops two loads and
    # two pointer adds per item; net prologue delta is zero (see CONST-init).
    old=' ldp x7, x8, [x22]\n ldp x9, x10, [x22, #16]\n'
    assert fwd.count(old)==1
    fwd=fwd.replace(old,'')
    old=' ldp x3, x4, [x23]\n'
    assert fwd.count(old)==1
    fwd=fwd.replace(old,' ldp x3, x4, [x21, #-32]\n')
    old=' ldp x5, x6, [x23, #16]\n'
    assert fwd.count(old)==1
    fwd=fwd.replace(old,' ldp x5, x6, [x21, #-16]\n')
    old=' stp x11,x12,[x21]\n stp x13,x14,[x21, #16]\n // inline k_mul;'
    assert fwd.count(old)==1
    fwd=fwd.replace(old,' // inline k_mul;')
    old=' stp x15,x16,[x21, #32]\n stp x24,x25,[x21, #48]\n add x23,x21,#32\n add x21,x21,#64\n'
    assert fwd.count(old)==1
    fwd=fwd.replace(old,' stp x15,x16,[x21, #32]\n stp x24,x25,[x21, #48]\n stp x13,x14,[x21, #16]\n stp x11,x12,[x21], #64\n')
    for dead in [' mov x22,x3\n',' sub x23,x2,#32\n']:
        assert fwd.count(dead)==1
        fwd=fwd.replace(dead,'')
    assert 'x22' not in fwd.split('pf_loop:')[1].split('pf_success:')[0]
    assert 'x23' not in fwd
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
        if index in (7,10):continue
        if index==3:k='k_abs'
        if index in (6,9):
            k='k_fused';dst=('x21',0 if index==6 else 32);b=('x23',192);store=True
        mapping={i:i for i in range(19)}
        saved_previous=(previous,previous_location)
        if chain=='use_plus':
            previous=cached_plus;previous_location=('x23',64)
        if chain=='save_plus':
            assert previous==[11,12,13,14] and previous_location==('x23',96)
            mapping.update({11:28,12:29,13:30,14:18,15:1,1:15,18:14})
            assert all(not re.match(r'\s*\w+ x(?:3|4|5|6|7|8|9|10),',x) for x in ks['k_abs'][4:] if not x.strip().startswith('stp'))
        elif chain:
            assert previous is not None and previous_location==a
            mapping={0:0,2:2,18:18};mapping.update(zip([3,4,5,6],previous));unused=set([1]+list(range(3,18)))-set(previous)
            for i in [1]+list(range(3,18)):
                if i in mapping:continue
                if i in unused:mapping[i]=i;unused.remove(i)
            for i in [1]+list(range(3,18)):
                if i not in mapping:mapping[i]=min(unused);unused.remove(mapping[i])
        if index==0:
            mapping.update({9:24,10:25,11:26,12:27})
        assert len(set(mapping[i] for i in [1]+list(range(3,18))))==16
        assert not ({mapping[i] for i in [1]+list(range(3,18))}&{0,2})
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
            def address(m):
                reg,off=m.groups();p,n=loc[reg];n+=int(off or 0)
                return '[PTR'+p[1:]+(', #'+str(n) if n else '')+']'
            raw=re.sub(r'\[(x[012])(?:,\s*#(-?\d+))?\]',address,raw)
            current=dict(mapping)
            if index in (5,8):current.update({7:24,8:25,13:26,14:27})
            raw=re.sub(r'\b([xw])(\d+)\b',lambda m:m[1]+str(current[int(m[2])]) if int(m[2])<19 else m[0],raw)
            raw=re.sub(r'PTR(\d+)',r'x\1',raw)
            for name in ['k_mul_store','k_sqr_store','k_add_store','abs_store']:raw=raw.replace(name,name+'_'+str(index))
            raw=raw.replace('fs_store','fs_store_'+str(index)).replace('fs_negative','fs_negative_'+str(index)).replace('fs_high','fs_high_'+str(index)).replace('fs_positive','fs_positive_'+str(index))
            lines.append(raw.strip() if raw.strip().endswith(':') else ' '+raw.strip());steps+=int(not raw.strip().endswith(':'))
        previous=[mapping[i] for i in OUT[k]] if not store else None
        previous_location=dst
        if chain=='save_plus':
            cached_plus=previous
            previous,previous_location=saved_previous
    lines+=[' sub x20,x20,#64',' sub x21,x21,#64',' subs x19,x19,#1',' b.ne pb_loop','pb_done:',' ldr x30,[sp],#16',' ret']
    fwd=fwd.replace(' cbz x19,pf_success\n',' cbz x19,pf_success\n mov w0,#-1\n mov x2,#0xffffffff00000001\n ldp x7, x8, [x3]\n ldp x9, x10, [x3, #16]\n')
    text=fwd+'\n'.join(lines)+'\n'
    text=text.replace(' cbz x19,pb_done\n',' cbz x19,pb_done\n mov w0,#-1\n mov x2,#0xffffffff00000001\n')
    text=text.replace('CONST_A','x0').replace('CONST_B','x2')
    for start,end in [('pf_loop:','pf_success:'),('pb_loop:','pb_done:')]:
        body=text.split(start)[1].split(end)[0]
        for line in body.splitlines():
            z=line.strip().split(None,1)
            if len(z)<2 or z[0].startswith(('b.','cb','tb')) or z[0] in ['cmp','tst','stp','str']:continue
            dst=z[1].split(',')[0].strip()
            assert dst not in ['x0','w0','x2','w2'],line
            if z[0]=='ldp':assert z[1].split(',')[1].strip() not in ['x0','w0','x2','w2'],line
    return text,{'removed_backward_memory_instructions':eliminated,'removed_forward_loads':4,'backward_body_static_instructions':steps+4,'forward_per_item_upper':127,'backward_per_item_upper':633,'calls_reordered':len(specs)}
def normalize_transform(text):
    f,b=text.split('.globl pair_bwd',1)
    old=' b.hs abs_store_fw0\n subs x11,xzr,x11\n sbcs x12,xzr,x12\n sbcs x13,xzr,x13\n sbc x14,xzr,x14\n'
    assert old in f;f=f.replace(old,'')
    hi=f.replace('pair_fwd','pair_fwd_hi');lo=f.replace('pair_fwd','pair_fwd_lo')
    # Four-limb operand swap, one limb per line: pair_fwd_lo must compute X-Y
    # on all four limbs (x3-x6 hold the Y-stream minuend, x7-x10 the X-point
    # subtrahend, so swapping every limb yields X-Y with a continuous borrow).
    lo=lo.replace('subs x11, x3, x7','subs x11, x7, x3')
    lo=lo.replace('sbcs x12, x4, x8','sbcs x12, x8, x4')
    lo=lo.replace('sbcs x13, x5, x9','sbcs x13, x9, x5')
    lo=lo.replace('sbcs x14, x6, x10','sbcs x14, x10, x6')
    for lab in re.findall(r'^([A-Za-z_]\w*):',f,re.M):
        if lab=='pair_fwd':continue
        lo=re.sub(r'\b'+lab+r'\b',lab+'_lo',lo);hi=re.sub(r'\b'+lab+r'\b',lab+'_hi',hi)
    old=' adc x1, xzr, xzr\n cbz x1,k_add_store_4\n adds x28,x28,#1\n sbcs x29,x29,x0\n sbcs x30,x30,xzr\n sbc x18,x18,x2\n'
    assert old in b;b=b.replace(old,'')
    return lo+hi+'.globl pair_bwd'+b
def check(root):
    text,meta=generate(root)
    assert normalize_transform(text)==(root/'pair_batch.s').read_text()
    meta.update(forward_per_item_upper=122,backward_per_item_upper=627,normalization_required=True,sorted_split_required=True)
    return meta
