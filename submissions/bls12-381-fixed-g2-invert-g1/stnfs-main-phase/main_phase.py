"""Conditional main-phase arithmetic; recorded toy campaigns are inputs."""
import math,time,json,base64,lzma,hashlib
from fractions import Fraction as F
R=52435875175126190479447740508185965837690552500527637822603658699938581184513

def stats(v):
 v=list(map(lambda x:F(str(x)),v));n=len(v);m=sum(v)/n
 return m,sum((x-m)**2 for x in v)/(n-1),n

def check(M,Ucap,Ulo,ROOT,run_toy=True,run_fields=True):
 t=time.perf_counter();Ucap,Ulo=int(Ucap),int(Ulo);ck={};out={}
 ck['count_arithmetic']=0<Ulo<=Ucap and Ulo==math.floor(F((1<<69)-3)/(69*F(7,10))-15*10**15-72*(1<<35)-(1<<20)-84)
 ck['model_caps']=F(str(M.NU300))==F(333,10000)and F(str(M.NU300_H))==F('0.03329')and M.SD_C==1 and M.P_CAP==1e-9 and M.FILT_EXCESS==1.08 and M.UNIT==20000
 def records(T):
  return {k:stats([F(v[1],1<<k)for v in rows])for k,rows in T.items()}
 fs,fx=records(M.FILT_RUNS),records(M.FILT_X)
 ceiling=F('0.03329')
 ck['filter_moment_support']=all(ceiling>=fx[k][0]and(ceiling-fx[k][0])**2>=4*fx[k][1]/fx[k][2]for k in(19,20))
 ck['filter_sd_support']=all(s*(1<<k)<=1 for T in(fs,fx)for k,(m,s,n)in T.items()if n>=6)
 ck['filter_coefficients']=all(v[2]<=255 for T in(M.FILT_RUNS,M.FILT_X)for rows in T.values()for v in rows)
 ck['large_scale_trend']=max(fs[k][0]for k in(19,20))<=fs[15][0]and max(fx[k][0]for k in(19,20))<=fs[15][0]
 paired={k:stats([F(v[1],v[2])for v in rows])for k,rows in M.FILT_Y.items()}
 weights={k:F(n,1)/s for k,(m,s,n)in paired.items()};w=sum(weights.values());pm=sum(weights[k]*paired[k][0]for k in paired)/w
 ck['changed_norm_support']=pm<=1 and(1-pm)**2>=4/w and all(v[3]<=255 for rows in M.FILT_Y.values()for v in rows)
 ks={k:stats(v)for k,v in M.K5_RUNS.items()};kc=F(str(M.C_ENTRY_K));km=F('9.966');k5_tax=F(16*10**6,Ulo)
 ck['k5_production_mean_support']=all(km>=ks[k][0]and(km-ks[k][0])**2>=4*ks[k][1]/ks[k][2]for k in(16,17))
 ck['k5_scaled_variance_support']=all(v*(1<<k)<256 for k,(m,v,n)in ks.items())
 ck['k5_tax_formula']=F(256,Ulo)/(F('9.97')-km)**2==k5_tax
 ck['k5_support']=all(kc>=ks[k][0]and(kc-ks[k][0])**2>=4*ks[k][1]/ks[k][2]for k in(16,17))
 ck['sample_sizes']=fx[19][2]>=16 and fx[20][2]>=8 and all(ks[k][2]>=8 for k in(16,17))
 mu=F(1087,1000)*Ucap;threshold=math.ceil(F(679,625)*Ucap)
 ck['descent_relation_gate']=threshold**40>=(1<<2543)and threshold>=math.ceil(F(27,25)*Ucap)
 short=mu/(mu-threshold)**2;filter_tax=F(10**10,Ulo);work_tax=F(1,1<<34);overflow=F(1,1<<64)
 ck['relation_mean_interval']=mu<(1<<64)
 ck['relation_overflow_tax']=F(1<<64,((1<<65)-(1<<64))**2)==overflow
 import lattice_count_bounds as LCB
 from pathlib import Path
 root=Path(ROOT);cfg=json.loads((root/'config.json').read_text());pd=json.loads((root/'pool_d.json').read_text())
 th,tp=F(cfg['production_theta_region']),F(cfg['production_theta']);co=LCB.correction(cfg,math.isqrt(int(th*cfg['production_A_squared'])))
 geo=M.Geometry(cfg,th,co['safe_yield_factor'],co['safe_work_population_factor']);nd=sum(pd['draws'][k]for k in pd['judged'])
 rt=M.pool_totals(M.add_sums([pd['sums'][k]for k in pd['judged']]),nd*float((th/tp)**3),nd,geo)
 lower=rt['log2_rels']-M.REL_BAND_LOSS+math.log2(1-2*rt['yield_se_log2']*math.log(2));upper=rt['log2_rels']+2*rt['yield_se_log2']
 ck['recorded_yield_support']=lower>=63.5761 and math.log2(float(mu))<63.5761 and upper<64
 st=json.loads(lzma.decompress(base64.b85decode((root/'stats16.b85').read_bytes())));campaign=st['main_phase_prereg600'];pre=campaign['prereg'];rr=campaign['rows']
 good=lambda row:row['checks_ok']and all(row['checks'].values())and(row['rank_A'],row['rank_O'],row['rank_R'],row['dimension'],row['m'],row['n'])==(79,80,80,80,72,5)
 ck['prereg600_records']=len(rr)==600 and[row['seed']for row in rr]==list(range(40951573,40952173))and all(good(row)for row in rr)
 ck['prereg600_source']=hashlib.sha256(campaign['source_text'].encode()).hexdigest()==pre['source_sha256']and(pre['first_seed'],pre['last_seed'],pre['seeds'],pre['pilot_seed_excluded'])==(40951573,40952172,600,20261107)
 zero=F(99,100)**600;ck['prereg600_probability']=zero<F(1,20)
 out['regularity_campaign']=dict(samples=len(rr),failures=sum(not good(row)for row in rr),source_sha256=pre['source_sha256'],zero_probability_exact=str(zero),zero_probability=float(zero),fresh_seed=40959999,scope=pre['scope'])
 expectation,cap=F('116.92454'),F('116.925')
 ck['work_abort_tax']=M.UNIT>(1<<14)and expectation+14>130 and(cap-expectation)*F(2,3)>F(1,4096)
 ck['required_ledger_rows']=F(str(M.LOG2_SM))>=103 and F(str(getattr(M,'LOG2_MAIN_ADMIN',0)))>=74 and F(str(getattr(M,'COLLECTION_CAP',0)))>=cap
 ck['filter_and_shortage_tax']=mu>threshold and filter_tax+short+work_tax+overflow+k5_tax<F(1,10**9)
 Ncap=F(333,10000)*Ucap;D=math.ceil(Ncap)+21
 sampler=F(5*D+1332,1<<256);eps=F(D,1<<128)+F(D+21,R)+sampler
 ck['projection_coin_tax']=eps<F(1,1<<66)and 18*72+3*12==1332 and sampler<F(1,1<<192)
 ck['field_and_dimensions']=(1<<254)<R<(1<<255)and D<R and M.N_SM==21 and(M.BW_M,M.BW_N)==(72,5)and D<6*(1<<56)
 admin=(1<<88)+(1<<81)+(1<<76)
 ck['administration_row']=admin<int(M.UNIT)*(1<<74)
 ck['projection_accumulator']=math.ceil(F(D,72))*((1<<128)-1)*(R-1)<(1<<448)
 sm_maps=(1<<66)+1;sm_ir=sm_maps*(1<<51)+(1<<70)
 ck['schirokauer_row']=255*72*73//2<(1<<20)and 2*72*72<(1<<14)and sm_ir<int(M.UNIT)*(1<<103)
 ck['descent_map_count_bound']=5*((1<<17)+2)<(1<<20)and F(1,1<<45)==F(1<<20,1<<65)
 ck['entry_cost_function']=F('9.2')*F('.86')+40*F('.14')==F('13.512')and F(str(M.W_ROW))*kc+1215+21*97==F('6252.97')
 def summaries(T):
  return {str(k):dict(n=n,mean=str(m),sample_variance_exact=str(s),variance_times_scale=str(s*(1<<int(k))),se=math.sqrt(float(s/n)),mean_plus_two_se=float(m)+2*math.sqrt(float(s/n)))for k,(m,s,n)in T.items()}
 out.update(filter_samples=summaries(fx),k5_samples=summaries(ks),paired_mean=str(pm),paired_se=math.sqrt(float(1/w)),filter_tax=str(filter_tax),shortage_tax=str(short),work_abort_tax=str(work_tax),relation_overflow_tax=str(overflow),k5_abort_tax=str(k5_tax),combined_cap_tax=str(filter_tax+short+work_tax+overflow+k5_tax),main_probability_lower=str(1-F(1,100)-F(1,1<<66)-F(1,10**9)),main_failure_upper=str(F(1,100)+F(1,1<<66)+F(1,10**9)),relation_gate=threshold,relation_mean_lower=str(mu),regularity_bad_mass_upper="1/100",recorded_yield_lower_log2=lower,recorded_yield_upper_two_se_log2=upper,relation_mean_upper=1<<64,retained_row_cap=1<<65,descent_map_cap=1<<65,global_map_cap=sm_maps,global_map_instruction_bound=sm_ir,descent_map_abort_tax=str(F(1,1<<45)),verification_admin_instructions=admin,sampler_abort_tax=str(sampler),projection_tax=str(eps),Ncap=str(Ncap),Dcap=D,Ulo=Ulo,Ucap=Ucap,admin_charge_log2=74,schirokauer_charge_log2=103)
 if run_fields:
  out['fields']=check_fields();ck['field_prerequisites']=out['fields']['unramified']and out['fields']['disc_h']==58673 and out['fields']['sig']['odd']==21
  out['protected_anchor']=protected_anchor();ck['protected_anchor']=out['protected_anchor']['checks_ok']
 if run_toy:
  out['toy']=toy(seed=40959999);ck['fixed_support_toy']=out['toy']['checks_ok']
 return dict(checks_ok=all(ck.values()),checks=ck,derived=out,requirements={'collection_mean_log2_upper':'116.92454','collection_cap_log2':'116.925','per_cell_random_work_log2_instructions':72,'retained_relation_cap':1<<65,'cold_schirokauer_cap':1<<65,'descent_schirokauer_cap':1<<65,'anchor_maps':1,'global_schirokauer_cap':(1<<66)+1,'descent_node_mean_upper':1<<17,'descent_psi_upper':'0.8','descent_sm_abort_tax':'1/35184372088832','relation_mean_upper':1<<64,'relation_mean_lower_factor':'1.087','relation_gate_factor':'1.0864','descent_input_relations_log2':'63.575','descent_meter_log2':'114.085','success_combination':'union','regularity_bad_mass_upper':'0.01','selected_square_k5_mean_upper':'9.966','selected_square_k5_variance_times_U_upper':256,'schirokauer_log2':103,'main_administration_log2':74},elapsed_seconds=time.perf_counter()-t,scope='Campaign moments and conditional coin bounds recomputed. Production H1/H2 capped means and independent per-cell work, H3-filter/H3c, fixed class/unit prerequisites, and H4u selected-square coefficient moments are hypotheses. The .01 production regularity/H-X bad-mass transport is an explicit distribution premise, not proved by the synthetic campaign or its binomial calculation. Filter moments apply conditional on the new Rneed gate using ALL retained rows; the recorded1.08-excess surrogate is support only. Overall main/descent failures must be combined by union, not by an independence product.')

def toy(seed=20261107,D=80,m=72,n=5):
 import random
 t=time.perf_counter();g=random.Random(seed)
 def uniform():
  for _ in range(256):
   v=g.getrandbits(255)
   if v<R:return v
  raise ValueError('capped sampler abort')
 A=[]
 for i in range(D):
  row=[0]*D
  for j in g.sample([j for j in range(D)if j!=i],6):row[j]=uniform()
  row[i]=-sum(row)%R;A.append(row)
 def mm(A,B):
  out=[];nb=len(B[0])
  for row in A:
   v=[0]*nb
   for k,x in enumerate(row):
    if x:
     for j,y in enumerate(B[k]):v[j]+=x*y
   out.append([x%R for x in v])
  return out
 def rank(A):
  a=[row[:]for row in A];p=0
  for j in range(len(a[0])):
   i=next((i for i in range(p,len(a))if a[i][j]),None)
   if i is None:continue
   a[p],a[i]=a[i],a[p];v=a[p];iv=pow(v[j],-1,R)
   for k in range(j,len(v)):v[k]=v[k]*iv%R
   for i in range(p+1,len(a)):
    f=a[i][j]
    if f:
     for k in range(j,len(v)):a[i][k]=(a[i][k]-f*v[k])%R
   p+=1
   if p==len(a):break
  return p
 a,b=math.ceil(D/m),math.ceil(D/n)
 X=[[0]*D for _ in range(m)]
 for j in range(D):X[m*j//D][j]=g.getrandbits(128)
 Y=[[uniform()for _ in range(n)]for _ in range(D)];powers=[Y]
 for _ in range(a+b-1):powers.append(mm(A,powers[-1]))
 seq=[mm(X,p)for p in powers]
 O=X+mm(X,A);B=[[x for p in powers[:b]for x in p[i]]for i in range(D)]
 AB=[[x for p in powers[1:b+1]for x in p[i]]for i in range(D)]
 H=[[seq[i+j][k][q]for j in range(b)for q in range(n)]for i in range(a)for k in range(m)]
 Hp=[[seq[i+j+1][k][q]for j in range(b)for q in range(n)]for i in range(a)for k in range(m)]
 ra,ro,rb=rank(A),rank(O),rank(B)
 checks=dict(kernel=all(sum(row)%R==0 for row in A),rank_A=ra==D-1,rank_O=ro==D,rank_R=rb==D,hankel=H==mm(O,B),shifted_hankel=Hp==mm(O,AB))
 return dict(checks_ok=all(checks.values()),checks=checks,rank_A=ra,rank_O=ro,rank_R=rb,dimension=D,m=m,n=n,a=a,b=b,seed=seed,seconds=time.perf_counter()-t,scope='One fresh fixed-contiguous-support synthetic sparse row-sum-zero matrix over the actual scalar field; no arithmetic relation, unit, or production-rank claim.')

def check_fields():
	import retain16 as RT,common12 as C
	Hr=[x%R for x in C.H]
	def res(a,b):
		a=RT._trim([x%R for x in a]);b=RT._trim([x%R for x in b]);s=1
		while b and len(b)>1:
			da,db=len(a)-1,len(b)-1;rem=RT.pmod(a,b,R)
			if not rem:
				return 0
			s=s*(-1)**(da*db)*pow(b[-1],da-(len(rem)-1),R)%R;a,b=b,rem
		return s*pow(b[0],len(a)-1,R)%R if b else 0
	def norm_poly(terms,deg):
		n=6*deg+1;xs=list(range(n))
		ys=[res(Hr,[sum(c[i]*pow(k,d,R)for d,c in terms if i<len(c))%R for i in range(6)])for k in xs]
		for j in range(1,n):
			for i in range(n-1,j-1,-1):
				ys[i]=(ys[i]-ys[i-1])*pow(xs[i]-xs[i-j],-1,R)%R
		poly=[ys[-1]]
		for i in range(n-2,-1,-1):
			poly=[(-xs[i]*poly[0]+ys[i])%R]+[(poly[k-1]-xs[i]*poly[k])%R for k in range(1,len(poly))]+[poly[-1]]
		return poly
	out={}
	for name,terms,deg in(('f',C.TERMS,12),('g',RT.G_TERMS,2)):
		poly=norm_poly(terms,deg)
		out['deg_'+name]=len(poly)-1
		out['sqf_'+name]=poly[-1]==1 and len(poly)==6*deg+1 and len(RT.pgcd(poly,RT.pderiv(poly,R),R))==1
	out['unramified']=out['sqf_f']and out['sqf_g']
	out['disc_h']=C.norm([-1,2,3,0,0,6])
	out['sig']=signatures()
	return out
def signatures():
    from fractions import Fraction

    def trim(p):
        p = list(p)
        while p and not p[-1]:
            p.pop()
        return p

    def rem(a, b):
        a = list(a)
        while a and len(a) >= len(b):
            shift = len(a) - len(b)
            q = a[-1] / b[-1]
            for j, bj in enumerate(b):
                a[shift + j] -= q * bj
            a = trim(a)
        return a

    def sturm(p):
        p = trim([Fraction(c) for c in p])
        seq = [p, trim([j * p[j] for j in range(1, len(p))])]
        while seq[-1]:
            r = rem(seq[-2], seq[-1])
            if not r:
                break
            seq.append([-c for c in r])
        assert len(seq[-1]) == 1  # squarefree; all endpoint conventions exact
        return seq

    def variations(signs):
        signs = [s for s in signs if s]
        return sum(a != b for a, b in zip(signs, signs[1:]))

    def sign(x):
        return (x > 0) - (x < 0)

    def counts(p):
        seq = sturm(p)
        minus = variations([sign(q[-1]) * (-1) ** (len(q) - 1) for q in seq])
        zero = variations([sign(q[0]) for q in seq])
        plus = variations([sign(q[-1]) for q in seq])
        assert p[0] != 0
        return minus - zero, zero - plus

    # Increasing-degree coefficient order.
    h = [-1, -1, 1, 1, 0, 0, 1]
    P = [1, 1, 0, 2, 0, -2, 1]
    assert counts(h) == (1, 1)
    assert counts(P) == (0, 0)

    # Verify f(x,y)=y^6 P(x^2/y)-h(y), coefficientwise in x.
    f_by_x_half_degree = {
        0: [1, 1, -1, -1], 1: [0, 0, 0, 0, 0, 1],
        3: [0, 0, 0, 2], 5: [0, -2], 6: [1],
    }
    for k in range(7):
        expected = [0] * 7
        expected[6 - k] = P[k]
        actual = list(f_by_x_half_degree.get(k, [])) + [0] * 7
        actual = actual[:7]
        if k == 0:
            actual = [a + b for a, b in zip(actual, h)]
        assert actual == expected

    # A real K_f or K_f^+ embedding would give a real zero of P; impossible.
    # g=x^2+|u|y has two real roots exactly over the unique negative real y.
    f, fp, g, hs = (0, 36), (0, 18), (2, 5), (2, 2)
    rank = lambda s: s[0] + s[1] - 1
    odd = rank(f) - rank(fp) + rank(g) - rank(hs)
    assert odd == 21
    return dict(f=f, fp=fp, g=g, h=hs, odd=odd)


def protected_anchor():
 import common12 as C,descent_witness as W
 t0=time.perf_counter();a=[1,0,0,0,0,0];b=a[:]
 nf,ng=C.norms(a,b);h=(W.P**12-1)//R;k=pow(h,-1,R)
 t=W.pow12((a,b),h);omega=W.pow12(t,k);sigma=(omega[0],[(-x)%W.P for x in omega[1]])
 p=1399;y=p-40
 hp=lambda z:(z**6+z**3+z*z-z-1)%p
 fp=lambda x,z:(x**12-2*z*x**10+2*z**3*x**6+z**5*x*x+1+z-z*z-z**3)%p
 ck=dict(prime_norm=nf==1399 and all(p%d for d in range(2,math.isqrt(p)+1)),
  ideal_root=hp(y)==fp(-1,y)==0 and fp(1,y)==0,
  cofactor_unit=h%R!=0 and h*k%R==1,
  idempotent=((h*k)**2-h*k)%(W.P**12-1)==0,
  subfield_killed=h%(W.P**6-1)==0,
  generator=omega!=W.ONE and W.pow12(omega,R)==W.ONE,
  sigma_inverse=W.mul12(omega,sigma)==W.ONE,
  omega_checkpoint=omega[1][0]==1823330539042453828028069116870546032896926197478728071363869126167756345638209516204631600716902807865027915471103)
 return dict(checks_ok=all(ck.values()),checks=ck,norm=nf,prime=p,y_root=y,x_root=p-1,cofactor_mod_r=h%R,cofactor_inverse_mod_r=k,omega_x0=omega[1][0],elapsed_seconds=time.perf_counter()-t0,scope='Exact protected alpha=1+x witness and idempotent GT projection. Class numbers and full unit-map injectivity remain conditional; this witness proves nonzero true core advice when its ideal column and all21 map coordinates are retained.')

if __name__=='__main__':
 import sys
 a=json.load(sys.stdin);seed=int(a.get('seed',40959999));n=int(a.get('samples',1))
 if not(0<=seed<=2**31-5 and 1<=n<=4):raise ValueError('seed0..2^31-5,samples1..4')
 out=[toy(seed=seed+i)for i in range(n)]
 print(json.dumps(dict(checks_ok=all(x['checks_ok']for x in out),results=out,scope='Fresh fixed-partition synthetic witnesses only; production regularity transport and class/unit premises are not proved.')))
