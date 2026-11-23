"""Generic-message repricing on the complete charged populations.

ReferenceC is unchanged; MSG is an internal scalar-service envelope.
"""
def compute(n,msg):
    c=n.c
    assert c.C==1000 and c.MSG==128 and not c.FAILED
    I=n.WS/6
    WQ=n.weighted_int(c.Q1G)
    cap_coeff=c.N_HOME*c.CAND_CAP
    candidate_coeff=c.ELIST+n.DEFL_PAIRS
    wlt_coeff=c.SG_FULL*(c.RHO*c.SG_FULL*c.D_WLT*(2+c.D_M)*c.KB_M+c.DWS*c.KB(132))
    fixed_msg_coeff=(2+2*c.G_GROUPS)*c.ALLN
    mine_fixed_coeff=c.N_POOLS*(1024*5+1+8*4)+c.ALLN+4*c.NODES*c.MX
    delivery_coeff=c.ALLN*(n.PROGRAM_BYTES//32)*2
    mine_det=c.W_MINE_DET_MAX+(n.UNWEIGHTED_VISIT-4)*c.S_MINE_TOT*c.CAP_H+n.CARD54_MINE
    exc_det=c.W_DET_MAX+n.positive_core+(n.UNWEIGHTED_VISIT-4)*c.SG_FULL*c.CAP_H
    

    d=msg-c.MSG
    assert d>=0
    resolve_new=2*(c.CAND_ITEM+64*d)+1
    assert resolve_new-c.RESOLVE==128*d
    barrier_new=max(c.BARRIER,7*(40*c.ALLN+c.BARRIER_MSGS*msg))
    barrier_pernode=n.cd(barrier_new,c.ALLN)
    generator_wave_delta=max(0,2*7*(40+2*msg)-8192)
    batch_delta=d+generator_wave_delta
    extra8_delta=8*max(0,barrier_pernode-4096)*(c.NODES+c.N_REL)
    batch_core_det=c.NODES*c.Q1GD*batch_delta+extra8_delta*c.Q1GD
    assert 4*c.PH4_FIX+4*msg+256<=c.EXCH_ROW//c.ALLN
    delta_mine_det=136*d*cap_coeff
    delta_exc_det=128*d*cap_coeff+batch_core_det
    rows={
        'program_delivery':delivery_coeff*d,
        'group_and_setup_messages':fixed_msg_coeff*d,
        'batch_coin_and_generator_barriers':c.NODES*WQ*batch_delta,
        'extra_eight_generator_relay_waves':extra8_delta*WQ,
        'candidate_resolve':128*d*candidate_coeff,
        'defective_label_resolve':128*d*wlt_coeff,
        'mine_candidate_resolve_and_eight_messages':136*d*c.E_CAND_MINE,
        'mine_tickets_reports_pool_results_x_distribution':mine_fixed_coeff*d,
        # Four inherited MSGs pay the128B reply. Pay the direct request too,
        # once for each selected X, even on paths that abort before fetching.
        'mine_selected_request_messages':8*c.N_POOLS*msg,
        'mine_exception_ceiling':n.ceilf(c._P_MINE_DEF*(mine_det+delta_mine_det))-n.DEF_MINE,
        'post_mine_exception_ceiling':n.ceilf(c.P_EXC*(exc_det+delta_exc_det))-n.EXC_WORK,
        'global_barriers':(c.NBAR_ALLOW+5)*(barrier_new-c.BARRIER),
    }
    assert all(v>=0 for v in rows.values())
    delta=sum(rows.values())
    direct=(c.ALLN*(n.PROGRAM_BYTES//32)*2*msg-delivery_coeff*c.MSG
        +(2+2*c.G_GROUPS)*c.ALLN*msg-fixed_msg_coeff*c.MSG
        +(c.ELIST+n.DEFL_PAIRS)*resolve_new-n.E_CAND_WORK
        +wlt_coeff*(resolve_new-c.RESOLVE)
        +c.E_CAND_MINE*(resolve_new+2*c.REQ_SERVE+8*msg)
          -c.E_CAND_MINE*(c.RESOLVE+2*c.REQ_SERVE+8*c.MSG)
        +c.N_POOLS*(c.PSEL_POOL+1024*(5*msg+32)+msg+4+8*(c.REQ_SERVE+5*msg))
          -c.N_POOLS*(c.PSEL_POOL+1024*(5*c.MSG+32)+c.MSG+4+8*(c.REQ_SERVE+4*c.MSG))
        +c.ALLN*(msg+c.STOP_CHECK)-c.ALLN*(c.MSG+c.STOP_CHECK)
        +c.NODES*c.MX*4*msg-c.NODES*c.MX*4*c.MSG
        +c.NODES*WQ*batch_delta+extra8_delta*WQ
        +rows['mine_exception_ceiling']+rows['post_mine_exception_ceiling']
        +(c.NBAR_ALLOW+5)*(barrier_new-c.BARRIER))
    assert direct==delta,(msg,str(direct-delta))
    return rows,barrier_new,batch_delta
