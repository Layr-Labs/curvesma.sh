"""Exact arithmetic only. Does not independently prove the supplied symbolic inequalities."""
from fractions import Fraction as F
from decimal import Decimal,localcontext
import json
r=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
NODES=1<<72;SIDE=1<<36;Q2=1<<27
P1=F(2*(r-1)+1,r*r);RHO=F(r,r-1)**4
KB=lambda f:F(2*(1<<(256-f)),r-1)
CAP_P=(1<<22)-(1<<12)
def dec(x):return Decimal(x.numerator)/Decimal(x.denominator)
def calc(q):
    n=2*q*Q2;S=NODES*n;gam=F(S*S-NODES*n*n,2);mu=gam*P1
    # Retain conservative 12-incidence/36-event constants from the source;
    # for two leaves the exact dependency/event bounds are smaller.
    cs=gam*12*(q+Q2)*S*(P1*P1+F(36,(r-1)**2))
    db=Q2*F(2,r-1)+F(1,r)
    mp=n*KB(36)+1;vp=mp+1+60*(q+Q2)*KB(36)
    dp=vp/(CAP_P-12-mp)**2
    m3=S*KB(96)+1;cap3=-(-m3.numerator//m3.denominator)+(1<<32)
    d3=(m3+1+60*(q+Q2)*KB(96))/(cap3-12-m3)**2
    intr=RHO*(S+12)*KB(144)
    h=RHO*(S*S*KB(72)*KB(184)+(2*S+20)*KB(184))
    g=2*RHO*n*S*KB(184)
    v=18*RHO*(F(S,SIDE)*S*KB(36)*KB(184)+(2*S+20)*KB(184))
    dc=h/((1<<28)-4);dr=g/((1<<29)-8);dv=v/((1<<30)-16)
    per=db+2*dp+d3+dc+dr+dv
    loss=mu*(1+F(1,r))*(F(1,r)+intr+2*per)
    # Eight 256-bit rejection draws per scalar. For this r, rejection <2^-31.
    # Charge a deliberately loose 2^-248 cap probability per scalar.
    assert F((1<<256)-r,1<<256)<F(1,1<<31)
    glob=NODES*q*db+F(NODES*(q+Q2+1),1<<248)
    mu_lo=F((mu.numerator<<200)//mu.denominator,1<<200)
    term=F(1);e_hi=term
    for k in range(1,41):
        term*=(-mu_lo)/k;e_hi+=term
    assert 0<mu_lo<2
    success=1-e_hi-cs-loss-glob
    return dict(n=n,S=S,mu=mu,success=success,cap3=cap3,mp=mp,cs=cs,
                bucket=cap3*16,intr=intr,loss=loss,global_loss=glob)
with localcontext() as ctx:
    ctx.prec=80
    lo,hi=1<<27,1<<28
    while hi-lo>1:
        m=(hi+lo)//2
        if calc(m)['success']>=F(3901,10000):hi=m
        else:lo=m
    a=calc(hi)
    hot=F(917,2)+8+F(73,4)+F(145,8)
    assert hot==F(4023,8)
    log=lambda x:dec(F(x)).ln()/Decimal(2).ln()
    rows={
      'paired_arithmetic_and_emission':NODES*hi*(917*Q2+39479+1024),
      'hop1_delivery_and_inflight':8*a['S'],
      'hop2_stream_delivery':F(73,4)*a['S']+NODES*(11*SIDE+3),
      'grouping_record_work':F(145,8)*a['S'],
      'two_cursor_initializations':NODES*(6*SIDE+2),
      'two_finishing_sweeps':NODES*(34*SIDE+6),
      'grouping_slots_tables':NODES*(12*SIDE+4+(1<<24)*15+F(7,2)*(1<<48)+2),
      'range_scanner_fixed_and_candidates':NODES*(24*(1<<48)+4*((1<<28)+1)+2),
      'program_node_initialization':NODES*(1<<20),
      'additional_bookkeeping_reserve':8*a['S'],
      'all_generator_points_and_headers':NODES*(hi+Q2+1)*(1<<27),
      'capped_scalar_sampling':NODES*(hi+Q2+1)*256,
      'resolution':NODES*(1<<28)*(1<<27),
      'source_request_service':NODES*(1<<29)*(1<<10),
      'four_phase4_count_exchanges':NODES*4*(1<<41),
      'barriers':NODES*32*(1<<10),
      'target_validation_and_output':(1<<28)+NODES*64,
    }
    total=sum(rows.values())
    assert F(total).denominator==1
    assert int(total)**10000 <= 1000**10000*(1<<1265279)
    setup=sum(v for k,v in rows.items() if k not in ['paired_arithmetic_and_emission','hop1_delivery_and_inflight','hop2_stream_delivery','grouping_record_work','additional_bookkeeping_reserve','resolution','source_request_service'])
    assert setup<=1000*(1<<118)
    assert all(x.denominator==1 for x in map(F,rows.values()))
    slotbytes=32*CAP_P*SIDE
    fixed=(1<<51)+2*(1<<39)+2*(1<<41)+hi*96+Q2*96+(1<<34)+(1<<44)+(1<<30)
    assert fixed<(1<<52)
    assert (1<<52)+2*slotbytes<(1<<64)
    assert (1<<52)+slotbytes+32*a['n']<(1<<64)
    assert a['cap3']*16<(1<<38)
    assert (1<<24)*(1<<38)+8*(1<<24)*a['cap3']<slotbytes
    print(json.dumps({'research_only':False,'q1':hi,'q2':Q2,
       'conditional_success_formula':str(dec(a['success'])),
       'hot_instruction_price_per_leaf':str(hot),
       'hot_only_log2_units':str(log(a['S']*hot/1000)),
       'full_draft_instructions':str(total),'full_draft_log2_units':str(log(total/1000)),
       'rows':{k:str(v) for k,v in rows.items()},'layout_checks':True,
       'S_log2':str(log(a['S'])),'cap3':a['cap3'],
       'group_bucket_bytes_log2':str(log(a['bucket'])),
       'scope':'Exact consequences of the supplied proof formulas; not independent mathematical certification.'}))
