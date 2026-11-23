"""Complete scalar callback surcharges on explicit full-size populations."""
def compute(n,gen=1<<20,resolve_reserve=1<<21,request_reserve=4096):
    c=n.c
    dg=gen-c.GEN_MUL
    dr=dg+request_reserve
    old_verify=256*(c.DBL+c.MADD_)+c.INV
    dv=4*(gen-old_verify)+resolve_reserve
    assert min(dg,dr,dv)>=0
    chain_count=sum(c.cd(q,c.M_G)for q in c.Q1G)
    full_chains=2*c.NODES*chain_count*dg
    mine_chains=2*c.NODES*c.NBATCH_MINE*dg
    cap=c.N_HOME*c.CAND_CAP
    candidate=c.ELIST+n.DEFL_PAIRS
    wlt=c.SG_FULL*(c.RHO*c.SG_FULL*c.D_WLT*(2+c.D_M)*c.KB_M+c.DWS*c.KB(132))
    # Ceil(a)+ceil(b)>=ceil(a+b). These additional ceilings deliberately
    # keep all earlier charged exceptional work and round independently.
    rows={
        'initial_Rg':c.NODES*dg,
        'mining_chains':mine_chains,
        'ordinary_chains_full_horizon':full_chains,
        'candidate_resolve':candidate*dv,
        'candidate_request':candidate*2*dr,
        'defective_label_resolve_request':wlt*(dv+2*dr),
        'mining_resolve_request':c.E_CAND_MINE*(dv+2*dr),
        'mining_selected_requests':c.N_POOLS*8*dr,
        'mining_exception':n.ceilf(c._P_MINE_DEF*(cap*(dv+2*dr)+mine_chains)),
        'post_mining_exception':n.ceilf(c.P_EXC*(cap*(dv+2*dr)+full_chains)),
    }
    assert all(v>=0 for v in rows.values())
    # Independently expand the grouped coefficient identity.
    direct=(c.NODES*dg+2*c.NODES*(chain_count+c.NBATCH_MINE)*dg
        +(candidate+wlt+c.E_CAND_MINE)*dv
        +2*(candidate+wlt+c.E_CAND_MINE+4*c.N_POOLS)*dr
        +rows['mining_exception']+rows['post_mining_exception'])
    assert sum(rows.values())==direct
    return rows,dict(gen_mul=gen,chain_add=c.CHAIN_ADD,
        chain_set=c.CHAIN_SET+2*dg,request_serve=c.REQ_SERVE+dr,
        resolve=c.RESOLVE+128*(n.MSG-c.MSG)+dv,
        scalar_verification=gen,ordinary_chains_per_G=chain_count,
        resolve_reserve=resolve_reserve,request_reserve=request_reserve)
