"""Exact arithmetic only. Does not independently prove the supplied symbolic inequalities."""
from fractions import Fraction as F
from decimal import Decimal,localcontext
import json
r=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
NODES=1<<70;SIDE=1<<35;Q2=1<<28
P1=F(2*(r-1)+1,r*r);RHO=F(r,r-1)**4
KB=lambda f:F(1<<(256-f),r-1)*(F(513,512) if f in (35,96) else F(2)+F(1,1<<30) if f==70 else F(2)+F(1,1<<14) if f in (144,182) else 2)
CAP_P=25*(1<<18)
def dec(x):return Decimal(x.numerator)/Decimal(x.denominator)
def ln_interval(x, terms=64):
    assert 1 <= x <= 2
    t=(x-1)/(x+1);t2=t*t;power=t;lo=F(0)
    for j in range(terms):
        lo+=2*power/(2*j+1);power*=t2
    tail=2*power/((2*terms+1)*(1-t2))
    return lo,lo+tail
def score_upper(total):
    x=F(total,1000*(1<<126));a,b=ln_interval(x);c,d=ln_interval(F(2))
    return F(126)+b/c
def calc(q):
    n=2*q*Q2;S=NODES*n;gam=F(S*S-NODES*n*n,2);mu=gam*P1
    cs=gam*12*(q+Q2)*S*(P1*P1+F(36,(r-1)**2))
    db=Q2*F(2,r-1)+F(1,r)
    qphase=(q+1)//2
    mp=2*qphase*Q2*KB(35)+1;vp=mp+1+60*(qphase+Q2)*KB(35)
    assert CAP_P-12>mp
    dp=vp/(CAP_P-12-mp)**2
    m3=S*KB(96)+1;cap3=-(-m3.numerator//m3.denominator)+(1<<30)
    assert cap3-12>m3
    d3=(m3+1+60*(q+Q2)*KB(96))/(cap3-12-m3)**2
    intr=RHO*(S+12)*KB(144)
    h=RHO*(S*S*KB(70)*(2*KB(182))+(2*S+20)*(2*KB(182)))
    g=2*RHO*n*S*(2*KB(182))
    v=18*RHO*(F(S,SIDE)*S*KB(35)*(2*KB(182))+(2*S+20)*(2*KB(182)))
    dc=h/((1<<28)-4);dr=g/((1<<29)-8);dv=v/((1<<30)-16)
    per=db+dp+d3+dc+dr+dv
    loss=mu*(1+F(1,r))*(F(1,r)+intr+2*per)
    assert F((1<<256)-r,1<<256)<F(1,1<<31)
    glob=NODES*q*db+F(NODES*(q+Q2+1),1<<248)
    mu_lo=F((mu.numerator<<200)//mu.denominator,1<<200)
    term=F(1);e_hi=term
    for k in range(1,41):
        term*=(-mu_lo)/k;e_hi+=term
    assert 0<mu_lo<2
    success=1-e_hi-cs-loss-glob-gam*F(1<<227,r*r)
    return dict(n=n,S=S,mu=mu,success=success,cap3=cap3,mp=mp,cs=cs,
                bucket=cap3*32,intr=intr,loss=loss,global_loss=glob)
base_calc=calc
def calc(q):
    a=base_calc(q)
    inc=q*Q2
    mean=inc*(F(1,2)+F(1,1<<116))
    cap=-(-mean.numerator//mean.denominator)+(1<<50)
    var=F(inc*(q+Q2-1),4)
    anchored=4*(q+Q2)
    assert cap-mean>anchored
    dm=var/(cap-mean-anchored)**2
    a['success']-=a['mu']*(1+F(1,r))*2*dm
    a.update(carry_cap=cap,carry_loss=dm,carry_mean=mean,carry_variance=var,anchored_incidence_bound=anchored)
    return a
with localcontext() as ctx:
    ctx.prec=80
    lo,hi=1<<28,1<<29
    while hi-lo>1:
        m=(hi+lo)//2
        if calc(m)['success']>=F(390000001,10**9):hi=m
        else:lo=m
    assert hi<(1<<29)
    a=calc(hi)
    hot=F(2994,8)+F(610,128)+F(97,4)+F(595,64)
    assert hot==F(6601,16)
    log=lambda x:dec(F(x)).ln()/Decimal(2).ln()
    rows={
      'paired_arithmetic_and_emission':NODES*hi*(F(2994,4)*Q2+39479+1024+34*87),
      'hop1_block_flush_delivery_inflight':F(610,128)*a['S'],
      'block_window_initialization_finishing_and_control':NODES*4096*SIDE,
      'hop2_stream_delivery':F(97,4)*a['S']+NODES*2*(11*SIDE+2),
      'grouping_record_work':F(595,64)*a['S'],
      'two_cursor_initializations':NODES*(6*SIDE+2),
      'two_finishing_sweeps':NODES*(34*SIDE+6),
      'grouping_slots_tables':NODES*((1<<26)*16+F(7,2)*(1<<48)+4),
      'range_scanner_fixed_and_candidates':NODES*(26*(1<<48)+8*((1<<28)+1)+2),
      'program_node_initialization':NODES*(1<<20),
      'all_generator_points_and_headers':NODES*(hi+Q2+1)*23441408,
      'capped_scalar_sampling':NODES*(hi+Q2+1)*256,
      'resolution':NODES*(1<<28)*42410496,
      'source_request_service':NODES*(1<<29)*(1<<10),
      'four_explicit_phase4_stages':NODES*4*(1<<41),
      'barriers':NODES*64*(1<<10),
      'admitted_xsum_corrections':NODES*4*a['carry_cap'],
      'range_tree_construction_and_driver':NODES*(8192*Q2+(1<<17)*hi),
      'normalization_sort_search':NODES*(1280*hi+3968*Q2+4352),
      'extra_phase_receiver_setup':NODES*(1<<14),
      'target_validation_and_output':(1<<28)+NODES*64,
    }
    total=sum(rows.values())
    assert F(total).denominator==1
    assert score_upper(int(total)) <= F(126218640596637371,10**15)
    assert F(1,(1<<34)-108)<=F(1,(1<<30)-16)
    setup=sum(v for k,v in rows.items() if k not in ['paired_arithmetic_and_emission','hop1_block_flush_delivery_inflight','hop2_stream_delivery','grouping_record_work','resolution','source_request_service'])
    assert setup<=1000*(1<<118)
    assert all(x.denominator==1 for x in map(F,rows.values()))
    slotbytes=32*CAP_P*SIDE
    fixed=(1<<51)+2*(1<<39)+2*(1<<41)+hi*96+Q2*96+(1<<34)+(1<<44)+(1<<30)+(4096*SIDE+32)+8*SIDE+64*SIDE+4096
    assert fixed<(1<<52)
    assert (1<<52)+slotbytes+40*(1<<26)*a["cap3"]<(1<<64)
    assert (1<<52)+slotbytes+32*a['n']<(1<<64)
    assert a['cap3']*32<(1<<38)
    assert fixed+(1<<30)+2048+192*Q2+32*hi+(1<<40)<(1<<52)  # added end table fits reservation
    print(json.dumps({'research_only':True,'q1':hi,'q2':Q2,'carry_cap':a['carry_cap'],'carry_loss':str(dec(a['carry_loss'])),
       'conditional_success_formula':str(dec(a['success'])),
       'hot_instruction_price_per_leaf':str(hot),
       'hot_only_log2_units':str(log(a['S']*hot/1000)),
       'full_draft_instructions':str(total),'full_draft_log2_units':str(log(total/1000)),
       'rows':{k:str(v) for k,v in rows.items()},'layout_checks':True,'memory_end':(1<<52)+slotbytes+40*(1<<26)*a['cap3'],
       'S_log2':str(log(a['S'])),'cap3':a['cap3'],
       'group_bucket_bytes_log2':str(log(a['bucket'])),
       'scope':'Fused-emission ledger: 2994/4/inc (113.25 fwd ping-pong + 635.25 fused bwd ldp), hop1 610/128, fixed 34*87; full revised argument review pending, not a score claim.'}))
