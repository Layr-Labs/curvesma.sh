"""expected number of special-q nodes in the two-log descent (DESCENT.md, proof.md section 7), by the same"""
def expected_nodes(Pm,dres):
	bins=dres['bins'];DET=dres['DET'];C={}
	def cmax(side,n):
		return max([C[(b,side)]for b in bins if b[0]<n-1e-9]+[0.0])
	def miss_c():
		return max([C[k]for k in C if k[0][0]<=Pm.THIN_MAX]+[0.0])
	for b in bins:
		for side in('f','g'):
			d=DET[(b,side)]
			sub=d['kA']*(cmax('f',d['nA'])+cmax('g',d['nA']))+d['kB']*(cmax('f',d['n'])+cmax('g',d['n']))+d['miss_children']*miss_c()
			C[(b,side)]=1.0+d['retry']*sub
	bt=dres['boot']
	sub=bt['kA']*cmax('g',bt['nA'])+cmax('g',float(bt['n1']))+bt['miss']*miss_c()
	return 2*bt['retry']*sub,C

# Grenie-Molteni, arXiv:1312.4463v4 Cor.1.4, conditional on GRH
# for the two fixed Dedekind zeta functions. Natural-log discriminant bounds.
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
GRH_DEGREES=(72,12)
GRH_LOG_DISC_UPPER=10000
GRH_ERROR_BITS=47
COUNT_CUTOFF_BITS=18
COUNT_BOUNDARY=sum(GRH_DEGREES)
@lru_cache(None)
def _count_bound():
 t=100
 lo=2*sum((F(1,(2*k+1)*3**(2*k+1))for k in range(t)),F(0))
 hi=lo+F(9,4*(2*t+1)*3**(2*t+1))
 # pi>3, 47<ln(B)<48, sqrt(B)<2^35; drop negative loglog terms.
 C=(F(1,6)+F(58,10*47))*GRH_LOG_DISC_UPPER+(F(48,24)+3)*max(GRH_DEGREES)+F(3,10)+F(133,10*47)
 q=F(2**69)/(69*lo-F(11,10))+2**COUNT_CUTOFF_BITS+2**GRH_ERROR_BITS+COUNT_BOUNDARY
 u=-(-q.numerator//q.denominator)
 checks=dict(ln2_interval=lo<hi<F(7,10),lnB=47<69*lo<69*hi<48,
 integral_cutoff=18*lo>F(121,10),grh_error=C<4096,
 ln3_upper=sum((F(11,10)**j/factorial(j)for j in range(13)),F(0))>3,
 f_root_bound=3**12>4*3**10+16*3**6+32*3**2+15,
 ln6_upper=sum((F(9,5)**j/factorial(j)for j in range(16)),F(0))>6,
 f_discriminant=72*71*F(9,5)<GRH_LOG_DISC_UPPER,
 g_root_bound=15132376222941642752<2**64,
 g_discriminant=4488*F(7,10)<GRH_LOG_DISC_UPPER,
 ceil_width=u-1<q<=u)
 assert all(checks.values())
 return u,checks

def factor_base_upper(bits=69):
 if bits!=69:raise ValueError('GRH factor-base certificate is for B=2^69')
 return _count_bound()[0]

def factor_base_certificate():
 u,checks=_count_bound()
 return dict(Ucap=u,checks=checks,degrees=GRH_DEGREES,
 log_discriminant_upper=GRH_LOG_DISC_UPPER,each_error_bits=GRH_ERROR_BITS)

def repeat_pool_certificate(bins=None):
 # Good absolute-degree-one ideals are a subset of admissible relative-degree-one
 # ideals. Divide their freely paired sigma orbits by two, independently per side.
 bins=tuple(bins)if bins is not None else((69,75),)+tuple((a,a+5)for a in range(75,225,5))
 t=100;lo=2*sum((F(1,(2*k+1)*3**(2*k+1))for k in range(t)),F(0))
 hi=lo+F(9,4*(2*t+1)*3**(2*t+1))
 C=(F(1,6)+F(58,10*47))*10000+(F(158,24)+3)*72+F(3,10)+F(133,10*47)
 lower=[]
 for a,b in bins:
  assert isinstance(a,int)and isinstance(b,int)and 69<=a<b<=225 and 5<=b-a<=6
  sq=1<<(-(-b//2))
  # Two endpoint errors, degree>=2 ideals, then polynomial-bad primes.
  lower.append((F(2**b-2**a,b)/hi-(8192+72)*sq-2**20)/2)
 checks=dict(grh_error_all_bins=C<4096 and 225*hi<158,
 bad_ideals=72*10000/lo<2**20,all_pools=all(q>2**68 for q in lower))
 assert bins and all(checks.values())
 return dict(checks=checks,bins_checked=len(bins),pool_log2_lower=68,
 min_pool_lower_floor=min(q.numerator//q.denominator for q in lower))

if __name__=='__main__':
 import sys,json,random,math
 from pathlib import Path
 h=Path(__file__).resolve().parent;h=h/'files'if(h/'files').is_dir()else h
 sys.path.insert(0,str(h))
 import descent_witness as W,descent16 as D,common12 as C
 a=json.load(sys.stdin);bits=a.get('bits',150);side=a.get('side','f');points=a.get('points',8);seed=a.get('seed',1)
 assert bits in(75,150,225)and side in('f','g')and 1<=points<=16 and 0<=seed<=1000000
 g=random.Random(seed);q,_,_,_,B,_=W.specialq_lattice(bits,side,g);R=D.PLAN16[(bits-5 if bits>75 else 69,side)][4];forms=W.forms();rows=[];c=(math.log2(q)+math.log2(D.CELLS[R]))/12
 for _ in range(points):
  k=W.int_ball_point(g,R);v=[sum(k[i]*B[i][j]for i in range(12))for j in range(12)];nf,ng=C.norms(v[:6],v[6:])
  rows.append(dict(cf=math.log2(nf)-72*c,cg=math.log2(ng)-12*c-W.LOG2U6,f_gap=W.cell_bound(v,forms[0])-math.log2(nf),g_gap=W.cell_bound(v,forms[1])-math.log2(ng)))
 print(json.dumps(dict(checks_ok=all(min(x['f_gap'],x['g_gap'])>=-1e-9 for x in rows),bits=bits,side=side,R=R,seed=seed,points=rows,scope='Fresh production-sized probable-prime lattice; exact integer-cube ball sampling and exact norms, floating reference bound gaps. Evidence only at these points, not a transport proof.')))
