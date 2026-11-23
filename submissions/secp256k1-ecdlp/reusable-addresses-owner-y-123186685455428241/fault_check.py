"""Static checks of the current expanded fault sites, vectors and stubs.

Inputs are the actual (instructions, labels) pairs produced by a64_load. This
function does not import or execute a candidate and makes no hardware claim.
"""
def run(leafpair,routepair):
    import re
    A='#0x100000000000'
    checks={};counts={}
    def require(name,condition):
        checks[name]=bool(condition)
        assert condition,('fault wiring',name)
    def unique(labels,name):
        require('unique_'+name,name in labels and len(labels[name])==1)
        return labels[name][0]
    def seq(program,pc,want,name):
        require(name,list(program[pc:pc+len(want)])==want)
    def op(mn,*args):return mn,list(args)
    leaf,ll=leafpair[:2];route,rl=routepair[:2]

    sites=[(name,unique(ll,name)) for name in ll if name.startswith('ks')]
    require('leaf_site_count',len(sites)==198)
    emitted={pc+3 for pc,(mn,args) in enumerate(leaf)
             if mn=='ldadd' and args[0]=='x4' and args[-1]=='[x0]'}
    require('leaf_sites_cover_slot_ldadds',emitted=={pc for name,pc in sites})
    stubs=set()
    for name,i in sites:
        require(name+'_two_stores',leaf[i][0]==leaf[i+1][0]=='stp')
        q,p,addr=leaf[i][1];s=addr[1:-1]
        require(name+'_registers',all(re.fullmatch(r'x(?:[0-9]|[12][0-9]|30)',r) for r in (s,p,q)))
        seq(leaf,i,[op('stp',q,p,'['+s+']'),op('stp','x2','x30','['+s+',#16]')],name+'_record')
        seq(leaf,i-3,[op('ldadd','x4',s,'[x0]'),op('bfxil',p,q,'#20','#21')],name+'_allocate')
        require(name+'_key_pack',leaf[i-1][0]=='bfxil' and leaf[i-1][1][0]==q
                and leaf[i-1][1][2:]==['#23','#41'])
        stub='gs'+s+p+q;j=unique(ll,stub);stubs.add(stub)
        seq(leaf,i+0x5000//4,[op('b',stub)],name+'_shadow')
        want=[op('mov','x4',A),op('stp',q,p,'[x4]'),op('stp','x2','x30','[x4,#16]'),
              op('mov','x4','#32'),op('sub',q,s,'#4064'),op('str',q,'[x0]'),
              op('ldr',p,'['+q+',#-16]'),op('add',p,p,'#1'),op('mov',s,'#4096'),
              op('stp',s,p,'['+q+',#-24]'),op('stlur','x0','['+q+',#-32]'),op('dmb','ish'),
              op('mrs',s,'elr_el1'),op('add',s,s,'#8'),op('msr','elr_el1',s),op('eret')]
        seq(leaf,j,want,stub+'_complete_shape')
    require('leaf_stub_inventory',len(stubs)==18 and stubs=={n for n in ll if re.fullmatch(r'gsx\d+x\d+x\d+',n)})
    counts.update(leaf_sites=len(sites),leaf_stubs=len(stubs),leaf_shadow_bytes=0x5000)

    sites=[(name,unique(rl,name)) for name in rl if name.startswith('kr')]
    require('relay_site_count',len(sites)==255)
    emitted={pc+1 for pc,(mn,args) in enumerate(route) if (mn,args)==op('ldadd','x17','x11','[x10]')}
    require('relay_sites_cover_slot_ldadds',emitted=={pc for name,pc in sites})
    stubs=set()
    for name,i in sites:
        require(name+'_two_stores',route[i][0]==route[i+1][0]=='stp')
        a,b,_=route[i][1];c,d,_=route[i+1][1]
        require(name+'_lane',tuple((a,b,c,d)) in (('x20','x21','x22','x23'),('x24','x25','x26','x27')))
        seq(route,i-2,[op('bfi','x10',a,'#3','#42'),op('ldadd','x17','x11','[x10]'),
            op('stp',a,b,'[x11]'),op('stp',c,d,'[x11,#16]')],name+'_record_and_allocate')
        stub='rs'+a;j=unique(rl,stub);stubs.add(stub)
        seq(route,i+0x20000//4,[op('b',stub)],name+'_shadow')
        want=[op('mov','x9',A),op('stp',a,b,'[x9]'),op('stp',c,d,'[x9,#16]'),
            op('sub','x9','x11','#8192'),op('add','x9','x9','#32'),op('str','x9','[x10]'),
            op('orr',c,'x10','#0x400000000000'),op('ldadd','x18','x11','['+c+']'),
            op('cmp','x11','x12'),op('b.lt','rs8'+a),op('str','x12','['+c+']'),op('b','rs9'+a),
            op('stp','x18','x19','[x9,#-16]'),op('str','x11','[x9,#-24]'),
            op('stlur','x10','[x9,#-32]'),op('dmb','ish'),op('mrs','x9','elr_el1'),
            op('add','x9','x9','#8'),op('msr','elr_el1','x9'),op('eret')]
        seq(route,j,want,stub+'_complete_shape')
        require(stub+'_sent_branch',unique(rl,'rs8'+a)==j+12)
        require(stub+'_return_branch',unique(rl,'rs9'+a)==j+16)
    require('relay_stub_inventory',stubs=={'rsx20','rsx24'})
    counts.update(relay_sites=len(sites),relay_stubs=len(stubs),route_shadow_bytes=0x20000)

    i=unique(rl,'kx');j=unique(rl,'xst')
    seq(route,i-1,[op('ldadd','x4','x5','[x0]'),op('stp','x14','x13','[x5]'),
                  op('stp','x2','x26','[x5,#16]')],'xemit_record_and_allocate')
    seq(route,i+0x20000//4,[op('b','xst')],'xemit_shadow')
    seq(route,j,[op('mov','x4',A),op('stp','x14','x13','[x4]'),op('stp','x2','x26','[x4,#16]'),
        op('mov','x4','#32'),op('sub','x16','x5','#4064'),op('str','x16','[x0]'),
        op('ldr','x17','[x16,#-16]'),op('add','x17','x17','#1'),op('mov','x5','#4096'),
        op('stp','x5','x17','[x16,#-24]'),op('stlur','x0','[x16,#-32]'),op('dmb','ish'),
        op('mrs','x5','elr_el1'),op('add','x5','x5','#8'),op('msr','elr_el1','x5'),op('eret')],
        'xemit_complete_stub')
    counts['xemit_sites']=1
    vectors=[]
    for program,labels,name,reg,offset in ((leaf,ll,'vbarg','x4',0x5000),
            (route,rl,'vbarx','x4',0x20000),(route,rl,'vbarr','x9',0x20000)):
        v=unique(labels,name)
        require(name+'_alignment',v%512==0)
        seq(program,v+128,[op('mrs',reg,'elr_el1'),op('add',reg,reg,hex(offset).join(['#',''])),
                          op('br',reg)],name+'_vector')
        vectors.append(dict(name=name,offset=offset,register=reg))
    # Bwd's SP is the immutable Y view. It has no literal store through SP;
    # runtime read-only guards check the broader dynamic address property.
    b0=unique(ll,'bwd');b1=unique(ll,'keys')
    stores={'str','stp','st64b','stlr','stlur','stadd','staddl','ldadd','swp'}
    bad=[pc for pc in range(b0,b1) if leaf[pc][0] in stores
         and any(re.match(r'\[sp(?:,|\])',x) for x in leaf[pc][1])]
    require('bwd_no_literal_sp_store',not bad)
    return dict(ok=True,**counts,vectors=vectors,checks=len(checks),
                bwd_no_literal_sp_store=True)
