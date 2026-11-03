from fractions import Fraction as F
from decimal import Decimal,localcontext
import json
r=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
NODES=1<<70;SIDE=1<<35;Q2=1<<28
P1=F(2*(r-1)+1,r*r);RHO=F(r,r-1)**4
KB=lambda f:F(1<<(256-f),r-1)*(F(513,512) if f in (35,96) else F(2)+F(1,1<<30) if f==70 else F(2)+F(1,1<<14) if f in (144,182) else 2)
CAP_P=25*(1<<18)
GF=(21,36,50,64)
CLAIM=F(126081134739098625,10**15)
WORST=F(126238773,10**6)
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
def calc(q,qg):
 n=2*q*Q2;S=NODES*n;gam=F(S*S-NODES*n*n,2);mu=gam*P1
 cs=gam*12*(q+Q2)*S*(P1*P1+F(36,(r-1)**2))
 db=Q2*F(2,r-1)+F(1,r)
 qphase=qg
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
 dc=h/((1<<24)-4);dr=g/((1<<29)-8);dv=v/((1<<30)-16)
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
def calc(q,qg):
 a=base_calc(q,qg)
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
def ends(q):return [q*f//64 for f in GF]
def gmax(q):
 e=[0]+ends(q);return max(e[i+1]-e[i] for i in range(len(GF)))
with localcontext() as ctx:
 ctx.prec=80
 lo,hi=1<<28,1<<29
 while hi-lo>1:
  m=(hi+lo)//2
  if calc(m,gmax(m))['success']>=F(390000001,10**9):hi=m
  else:lo=m
 assert hi<(1<<29)
 E=ends(hi);qg=gmax(hi);a=calc(hi,qg);G=len(GF)
 PHI=[F(1)]+[min(F(1),1-calc(b,qg)['success']) for b in E[:-1]]
 assert all(0<=PHI[g+1]<=PHI[g]<=1 for g in range(G-1))
 hot=F(2942,8)+F(581,128)+F(85,4)+F(595,64)
 log=lambda x:dec(F(x)).ln()/Decimal(2).ln()
 X=F(2942,4)*Q2+39479+1024+1623023+256+(1<<17)+1280+2*Q2*(F(581,128)+F(85,4))
 GRP=2048*SIDE+(11*SIDE+2)+(14*SIDE+10)+(3*SIDE+1)+(17*SIDE+3)+64*1024+(1<<14)
 SCAN_FIX=11*(1<<26)+1+F(5,2)*(1<<48)+1+25*(1<<48)+8*((1<<24)+1)+2+4*(1<<41)
 SCAN_RES=(1<<24)*20592111+(1<<29)*(1<<10)
 REC=F(595,64)*2*Q2
 ONCE=(1<<48)+1+5*(1<<26)+1+(1<<20)+Q2*1577967+1623023+(Q2+1)*256+8192*Q2+3968*Q2+4352+1+4*a['carry_cap']+64
 e0=[0]+E
 W=[(e0[g+1]-e0[g])*X+GRP+e0[g+1]*REC+SCAN_FIX+SCAN_RES+(1024 if g<G-1 else 0) for g in range(G)]
 assert all(F(w).denominator==1 for w in W)
 Ework=NODES*(ONCE+sum(PHI[g]*W[g] for g in range(G)))+(1<<28)
 worst=NODES*(ONCE+sum(W))+(1<<28)
 total=-(-Ework.numerator//Ework.denominator)
 assert score_upper(total)<=CLAIM<score_upper(total)+F(1,10**15)
 assert score_upper(worst)<=WORST
 setup=NODES*(ONCE+sum((e0[g+1]-e0[g])*(1623023+256+(1<<17)+1280)+GRP+SCAN_FIX for g in range(G)))+(1<<28)
 assert setup<=1000*(1<<118)
 assert F(1,(1<<34)-108)<=F(1,(1<<30)-16)
 slotbytes=32*CAP_P*SIDE
 fixed=(1<<51)+2*(1<<39)+2*(1<<41)+hi*96+Q2*96+(1<<34)+(1<<44)+(1<<30)+(1<<30)+(8192*SIDE+32)+(1<<38)+8*SIDE+64*SIDE+4096+4096
 assert fixed<(1<<52)
 mem_end=(1<<52)+slotbytes+40*(1<<26)*a['cap3']+(1<<51)
 assert mem_end<(1<<64)
 assert (1<<52)+slotbytes+32*a['n']<(1<<64)
 assert a['cap3']*32<(1<<38)
 assert fixed+(1<<30)+2048+192*Q2+32*hi+(1<<40)+20480*G<(1<<52)
 print(json.dumps({'research_only':True,'q1':hi,'q2':Q2,'groups':G,'group_ends':E,'largest_group':qg,
       'PHI':[str(dec(p)) for p in PHI],'carry_cap':a['carry_cap'],'carry_loss':str(dec(a['carry_loss'])),
       'success_formula':str(dec(a['success'])),'hot_per_leaf':str(hot),
       'expected_instructions_ceil':str(total),'expected_log2_units':str(log(F(total,1000))),
       'worst_instructions':str(worst),'worst_log2_units':str(log(F(worst,1000))),
       'group_work':[str(w) for w in W],'once_per_node':str(ONCE),'layout_checks':True,'memory_end':mem_end,
       'S_log2':str(log(a['S'])),'cap3':a['cap3'],
       'scope':'Exact grouped expected-work ledger; see proof.md section 15.'}))
