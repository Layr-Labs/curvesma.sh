from fractions import Fraction as F
from math import isqrt
r=52435875175126190479447740508185965837690552500527637822603658699938581184513
lam=228988810152649578064853576960394133503
assert lam*lam+lam+1==r
d=2**20;m=(r-1)//d
n=isqrt(11*m//120)
while 120*n*n<11*m:n+=1
assert 120*(n-1)**2<11*m<=120*n*n
X1=2**82;X2=2**70;Z=X1*X2
NA=196*n*d//X2+2;M=(NA+1)//2
ln2lo=sum((F(2,(2*j+1)*3**(2*j+1)) for j in range(64)),F(0))
Pup=(F(X2,1)/(70*ln2lo-F(11,10)))
P=Pup.numerator//Pup.denominator
assert P==24896332478097058409
assert 3*X1*X1<r and 2*X1*X1*X2<r
base=F(36*n**4,m)
D=lambda Y:2*Y*Y//r+1
Amax=(F(228,625)*X1*X1).__floor__()
ZB=(F(57,100)*Z).__floor__()
ZG=(F(16,25)*Z).__floor__()
inner_full=48*P*P*D(X1*X1)
inner_tight=48*P*P*D(Amax)
twoL=96*M*M*D(X2*X2)
oneLoneS=192*M*P*D(Z)
oneLoneS_tight=96*M*P*(D(ZB)+D(ZG))
twoS=48*P*(P+1)*D(X1*X1)
twoS_tight=48*P*(P+1)*D(Amax)
out={'status':'PASS_EXACT_ARITHMETIC_ANALYTIC_PRIMECOUNT_INPUT_DECLARED',
'analytic_input':'Dusart2010 Thm6.9 eq6.6, pi(x)<=x/(log x-11/10) for x>60184; not re-proved',
'prime_count_input_source':'https://arxiv.org/pdf/1002.0442',
'pi_X2_upper':P,'r':r,'lambda':lam,'n':n,'M':M,'X1':X1,'X2':X2,
'inner_pruned_noL_exactly1S':{'full_box_upper':inner_full,'full_box_ratio_exact':str(F(inner_full,1)/base),'full_box_ratio_decimal':float(F(inner_full,1)/base),
'tight_endpoint_upper':inner_tight,'tight_endpoint_ratio_exact':str(F(inner_tight,1)/base),'tight_endpoint_ratio_decimal':float(F(inner_tight,1)/base)},
'multi_wrapped_raw':{
 'twoL_ratio_exact':str(F(twoL,1)/base),'twoL_ratio_decimal':float(F(twoL,1)/base),
 'oneLoneS_ratio_exact':str(F(oneLoneS,1)/base),'oneLoneS_ratio_decimal':float(F(oneLoneS,1)/base),
 'oneLoneS_tight_ratio_exact':str(F(oneLoneS_tight,1)/base),'oneLoneS_tight_ratio_decimal':float(F(oneLoneS_tight,1)/base),
 'twoS_ratio_exact':str(F(twoS,1)/base),'twoS_ratio_decimal':float(F(twoS,1)/base),
 'twoS_tight_ratio_exact':str(F(twoS_tight,1)/base),'twoS_tight_ratio_decimal':float(F(twoS_tight,1)/base),
 'full_total_ratio_exact':str(F(twoL+oneLoneS+twoS,1)/base),'full_total_ratio_decimal':float(F(twoL+oneLoneS+twoS,1)/base),
 'tight_total_ratio_exact':str(F(twoL+oneLoneS_tight+twoS_tight,1)/base),'tight_total_ratio_decimal':float(F(twoL+oneLoneS_tight+twoS_tight,1)/base)},
'combined_inner_plus_multi_tight_ratio_exact':str(F(inner_tight+twoL+oneLoneS_tight+twoS_tight,1)/base),
'combined_inner_plus_multi_tight_ratio_decimal':float(F(inner_tight+twoL+oneLoneS_tight+twoS_tight,1)/base),
'Amax':Amax,'D_Amax':D(Amax),'D_X1sq':D(X1*X1),'D_Z':D(Z),
'scope':'Inner term uses pruned index set and exactly-one-small/no-large gcd. Multi terms raw hence also pruned; coprime and same-side exactly-one-small are not bounded here.'}
