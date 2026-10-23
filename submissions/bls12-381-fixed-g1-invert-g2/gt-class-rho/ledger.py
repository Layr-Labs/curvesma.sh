"""Recomputes every number of proof.md: the step count, the cycle bounds, the itemised charges, the reserve, the claimed score
and the memory bound. Standard library only. Run: python3 ledger.py"""
from math import log2, log, lgamma, ceil, sqrt
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
N = (r-1)//12 + 1
R = 2**16; th = 2.0**-60; M = 2.0**20; c = 1.4; P = 24
MU = 1.20            # explicit extra margin on steps
PEN = sqrt(1+1/R)    # BLS11 additive-walk penalty
RES = 1.07           # explicit unallocated reserve on the total
cap = 20/th
L2 = log2
k0 = sqrt(2*log(1/0.61)*N)
K = MU*PEN*k0
print("log2 r %.5f log2 N' %.5f log2 k0 %.5f log2 K %.5f"%(L2(r),L2(N),L2(k0),L2(K)))
print("model success at margin: 1-0.61^(MU^2) = %.4f ; tolerates N_real/N_eff up to %.3f"%(1-0.61**(MU**2), MU**2))
def lC(n,b): return (lgamma(n+1)-lgamma(b+1)-lgamma(n-b+1))/log(2)
def lrho(Lc):
    best=[lC(R,b) + Lc*L2(b/R) - b*L2(12) for b in range(1, min(Lc//2, R)+1)]
    m=max(best); return m + L2(sum(2**(x-m) for x in best))
rho = {2: 1/(12*R)}
for Lc in range(3,P+1): rho[Lc] = 2**lrho(Lc)
omega_short = sum(rho[Lc]*(Lc+3) for Lc in range(2,P+1))
print("rho2 2^%.3f rho3 2^%.2f rho4 2^%.2f rho24 2^%.2f"%(L2(rho[2]),L2(rho[3]),L2(rho[4]),L2(rho[24])))
print("omega_short 2^%.3f  (L=2 part 2^%.3f)"%(L2(omega_short), L2(5*rho[2])))
worst=-1e9; arg=None
for Lc in list(range(25,400))+[int(2**(x/4)) for x in range(36,200)]:
    v=lrho(Lc)
    if v>worst: worst,arg=v,Lc
print("max over L>=25 (L up to 2^50) of log2 rho_L bound: %.2f at L=%d"%(worst,arg))
print("rho at L=2R: 2^%.0f  L=2^40: 2^%.0f"%(lrho(2*R), lrho(2**40)))
omega_long = cap*(2**worst)*cap
print("omega_long (cap x sum_L rho_L, crude) 2^%.1f"%L2(omega_long))
trails = K*th + M
closure_events = trails*(2/th**2)/(2*N)
print("closure: expected trapped trails 2^%.2f, x cap = 2^%.2f steps"%(L2(closure_events), L2(closure_events*cap)))
omega = omega_short + omega_long
Kexec = K*(1+omega)
steps = Kexec*c
print("Kexec 2^%.5f steps*c 2^%.5f"%(L2(Kexec),L2(steps)))
others = {"starts": (Kexec*th+M)*512, "dp": Kexec*th*0.5, "delay": M*cap*c, "recovery": 2*cap*c*4,
 "closure_waste": closure_events*cap*c, "table": 2.0**45.1, "misc": 2.0**15.1}
for n_,v in others.items(): print(" ", n_, "2^%.2f"%L2(v))
W0 = steps + sum(others.values())
W = W0*RES
print("W before reserve 2^%.5f ; with reserve %.2f 2^%.5f"%(L2(W0),RES,L2(W)))
claim = ceil(L2(W)*100-1e-9)/100
print("claim:", claim, "slack over unreserved ledger: %.4f bits ; over reserved %.4f"%(claim-L2(W0), claim-L2(W)))
print("delay rel 2^%.1f"%L2(M*cap/K))
rec = Kexec*th + M
tab = 12*R*576; win=P*704
pp = tab+R*32+win+4096
print("records 2^%.3f ; 160 B each: 2^%.3f ; per-proc 2^%.3f; x M 2^%.3f ; total 2^%.3f"%(L2(rec), L2(rec*160), L2(pp), L2(M*pp), L2(rec*160+M*pp)))
print("lemma: 12*2^76/(p-1) = 2^%.2f ; union over Kexec: 2^%.2f"%(L2(12*2**76/(p-1)), L2(12*2**76/(p-1))+L2(Kexec)))
print("table: 2^%.3f + 2^%.3f = 2^%.3f"%(L2(2.0**36*512), L2(11*2.0**36), L2(2.0**36*512+11*2.0**36)))
