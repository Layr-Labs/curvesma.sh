"""Exact first-Montgomery-product reduction removal and register permutation."""
import re

def generate(text):
    mark='// ---- inlined blst __mul_mont_384 (verbatim, ret dropped)'
    a=text.index(mark);b=text.index('// ---- end __mul_mont_384',a)
    body=text[a+len(mark):b];lines=body.splitlines()
    code=[i for i,l in enumerate(lines) if l.split('//')[0].strip()]
    cut=code[-13];tail=[lines[i].split('//')[0].strip() for i in code[-13:]]
    assert len(tail)==13 and tail[0].startswith('subs') and all(x.startswith('csel') for x in tail[-6:])
    assert tail[-1].startswith('csel\tx16,')
    # Swap banks 11-16/19-24 so the unreduced output lands in x11..x16.
    mp={i:i for i in range(31)}
    for i in range(6):mp[11+i]=19+i;mp[19+i]=11+i
    assert len(set(mp.values()))==31
    out=[]
    for line in lines[:cut]:
        ins,sep,comment=line.partition('//')
        ins=re.sub(r'\b([xw])(\d+)\b',lambda m:m[1]+str(mp[int(m[2])]),ins)
        out.append(ins+sep+comment)
    prefix=text[:a]
    # y1+y2 goes directly into the renamed input bank.
    for i,op in enumerate(['adds','adcs','adcs','adcs','adcs','adc']):
        old=f'\t{op}\tx{11+i},x{11+i},x{19+i}'
        new=f'\t{op}\tx{19+i},x{11+i},x{19+i}'
        assert prefix.count(old)==1
        prefix=prefix.replace(old,new)
    result=prefix+'// ---- lazy first Montgomery product, exact register permutation\n'+'\n'.join(out)+'\n\t'+text[b:]
    # Fuse the second raw product with the y subtraction.
    a=result.index(mark);b=result.index('// ---- end __mul_mont_384',a)
    body=result[a+len(mark):b];ls=body.splitlines()
    live=[i for i,l in enumerate(ls) if l.split('//')[0].strip()]
    tail=[ls[i].split('//')[0].strip() for i in live[-13:]]
    assert tail[0].startswith('subs') and all(x.startswith('csel') for x in tail[-6:])
    result=result[:a]+'// ---- fused raw second Montgomery product\n'+'\n'.join(ls[:live[-13]])+'\n\t'+result[b:]
    a=result.index('\tldr\tx0,[x18,#-24]',a)
    b=result.index('\tstp\tx11,x12,[x2,#48]',a)
    old=[l.split('//')[0].strip() for l in result[a:b].splitlines() if l.split('//')[0].strip()]
    assert len(old)==23
    block=['\tldr x0,[x18,#-24]','\tldp x11,x12,[x0,#48]','\tldp x13,x14,[x0,#64]','\tldp x15,x16,[x0,#80]']
    for i in range(6):block.append('\t'+('subs' if i==0 else 'sbcs')+f' x{11+i},x{19+i},x{11+i}')
    block.append('\tb.hs .Llazy_y_nonnegative')
    for i in range(6):block.append('\t'+('adds' if i==0 else 'adc' if i==5 else 'adcs')+f' x{11+i},x{11+i},x{5+i}')
    block+=['\tb .Llazy_y_store','.Llazy_y_nonnegative:']
    for i in range(6):block.append('\t'+('subs' if i==0 else 'sbcs')+f' x{19+i},x{11+i},x{5+i}')
    for i in range(6):block.append(f'\tcsel x{11+i},x{11+i},x{19+i},lo')
    block.append('.Llazy_y_store:')
    result=result[:a]+'\n'.join(block)+'\n'+result[b:]
    return result,{'minimum_saved_per_lane':26,'additional_saving_on_negative_y_difference':5,'register_permutation_injective':True,'second_raw_multiply_body_unchanged':True}

def check(root):
    text,meta=generate((root/'measure/side_lanes.S').read_text())
    assert text==(root/'side_lazy.S').read_text()
    return meta
