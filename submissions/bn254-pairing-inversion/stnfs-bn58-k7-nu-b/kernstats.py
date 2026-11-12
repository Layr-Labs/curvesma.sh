"""kernstats.py: pure-Python recomputation of the bn1 kernel prices from the shipped archive text sections (no numpy).
  rec7_k(text)   -> [(stratum lower edge, mean + 2 SE of kappa)], decreasing edges (model16.REC7_K form), from 'K7BN k lo n sum sumsq max' lines (values x 1000, rounded up)
  nb_prices(text)-> (C_NB_EVAL, C_NB_ROW parts) from 'NB2BN kappa|eta n sum sumsq max' lines (values x 10, rounded up)
  norms_max(text)-> largest Ir of the exact-norm kernel and the check totals from 'NORMS pairs ... max M' / 'NCHECK ...' lines
  prep_max(text) -> largest basis-preparation Ir from 'PREP7BN st n sum max' lines"""
import math,re
def _m2se(n,s,q):m=s/n;sd=math.sqrt(max(0.,(q-n*m*m)/(n-1)));return m,sd,m+2*sd/math.sqrt(n)
def rec7_k(text):
    out=[]
    for l in text.splitlines():
        t=l.split()
        if t[:1]==['K7BN']:m,sd,u=_m2se(int(t[3]),int(t[4])/1e3,int(t[5])/1e6);out.append((float(t[2]),math.ceil(u*1e4)/1e4,int(t[6])/1e3))
    return sorted(out,reverse=True)
def nb_prices(text):
    r={}
    for l in text.splitlines():
        t=l.split()
        if t[:1]==['NB2BN']:m,sd,u=_m2se(int(t[2]),int(t[3])/10,int(t[4])/100);r[t[1]]=dict(n=int(t[2]),mean=m,sd=sd,up=u,max=int(t[5])/10)
    return r
def norms_max(text):
    mx=max(map(int,re.findall(r"max (\d+)",text)));pairs=sum(map(int,re.findall(r"NCHECK pairs (\d+) mismatches 0 L-violations 0",text)));bad=len(re.findall(r"NCHECK",text))-len(re.findall(r"NCHECK pairs \d+ mismatches 0 L-violations 0 .*recheck-flag errors 0",text))
    return mx,pairs,bad
def prep_max(text):return max(int(l.split()[4])for l in text.splitlines()if l.startswith('PREP7BN'))
