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

# Direct GRH smoothed-explicit-formula proof is reproduced in proof.md.
# No quantitative external prime-ideal or zero-count theorem is assumed.
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
GRH_DEGREES=(72,12)
GRH_LOG_DISC_UPPER=10000
GRH_ERROR_BITS=55
COUNT_AVERAGE_ERROR=F(17523154306377711616,1175)
COUNT_CUTOFF_BITS=18
COUNT_BOUNDARY=sum(GRH_DEGREES)
def prime_count_error_bits(b):
 if not isinstance(b,int)or not 69<=b<=225:raise ValueError('endpoint bits69..225')
 return (3*b+3)//4+3
@lru_cache(None)
def _count_bound():
 t=100
 lo=2*sum((F(1,(2*k+1)*3**(2*k+1))for k in range(t)),F(0))
 hi=lo+F(9,4*(2*t+1)*3**(2*t+1))
 q=F(2**69)/(69*lo-F(11,10))+2**COUNT_CUTOFF_BITS+COUNT_AVERAGE_ERROR+COUNT_BOUNDARY
 u=-(-q.numerator//q.denominator)
 checks=dict(ln2_interval=lo<hi<F(7,10),lnB=47<69*lo<69*hi<48,
 integral_cutoff=18*lo>F(121,10),zero_sum=10003<101**2,
 smoothing_step=F(202,2048)<F(1,10),
 residue_derivative=20000+6+3*72+72*226*F(7,10)<2**15,
 psi_coefficient=101+F(221,200)*101+F(1,2**18)<213,
 prime_powers=F(72*44,2048)*F(7,10)<2,
 theta_coefficient=213+2<256,log_cutoff=44*lo>30,
 pi_coefficient=F(74,256)+256*(F(1,47)+F(4,2700))<8,
 endpoint_error=prime_count_error_bits(69)==GRH_ERROR_BITS,
 field_zero_sums=72*71*F(9,5)+3<96**2 and 4488*F(7,10)+3<57**2,
 field_psi=F(421,200)*96+F(1,2**18)<203 and F(421,200)*57+F(1,2**18)<121,
 field_prime_powers=F(72*44,2048)*F(7,10)<2 and F(12*44,2048)*F(7,10)<1,
 quarter_root=F(27,16)**4>8,
 mean_field_error=44*2**44+F(327,2)*F(27,16)*2**51*(F(1,47)+F(4,2700))==COUNT_AVERAGE_ERROR,
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
 if bits!=69:raise ValueError('factor-base certificate is for B=2^69')
 return _count_bound()[0]

def factor_base_certificate():
 u,checks=_count_bound()
 return dict(Ucap=u,checks=checks,degrees=GRH_DEGREES,
 log_discriminant_upper=GRH_LOG_DISC_UPPER,each_error_bits=GRH_ERROR_BITS,average_factor_error_upper=str(COUNT_AVERAGE_ERROR),
 method='Direct smoothed Dedekind explicit formula; full derivation in proof.md')

def repeat_pool_certificate(bins=None):
 bins=tuple(bins)if bins is not None else((69,75),)+tuple((a,a+5)for a in range(75,225,5))
 t=100;lo=2*sum((F(1,(2*k+1)*3**(2*k+1))for k in range(t)),F(0))
 hi=lo+F(9,4*(2*t+1)*3**(2*t+1))
 lower=[]
 for a,b in bins:
  assert isinstance(a,int)and isinstance(b,int)and 69<=a<b<=225 and 5<=b-a<=6
  sq=1<<(-(-b//2))
  errors=2**prime_count_error_bits(a)+2**prime_count_error_bits(b)
  lower.append((F(2**b-2**a,b)/hi-errors-72*sq-2**20)/2)
 checks=dict(direct_endpoint_bounds=all(prime_count_error_bits(b)==(3*b+3)//4+3 for ab in bins for b in ab),
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

# Absolute-degree-one, simple free-sigma subset gives a separate actual-U lower.
def factor_base_lower():
 q=F(2**69-3,1)/(69*F(7,10))-15*10**15-72*2**35-2**20-84
 u=q.numerator//q.denominator
 assert u==12206442875139313724 and u>10**19
 return u
