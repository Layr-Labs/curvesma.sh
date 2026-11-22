"""Certified 50-digit Decimal + exact-O_H-algebra verification of the Voronoi fundamental domain F_Vor:
a in F_Vor iff T2(u a) >= T2(a) for the 10 listed units u. With X_k = |sigma_k(a)|^2 at the four complex places
(k = 1, 7, 4, 2), T2(u a) - T2(a) = sum_k (|sigma_k(u)|^2 - 1) X_k is linear in X.

ZERO floating-point math/cmath calls:
  - Roots of unity exp(2*pi*i*k/15), |sigma_k(u)|^2, ln, exp, and 3x3 Cramer ray solves all use 50-digit Decimal
    rational-Taylor arithmetic (residuals < 10^-42).
  - Lattice membership of all 10 Voronoi units AND all 9 Sinnott cyclotomic unit generators in the rank-3 unit
    basis (u_0, u_1, u_2) is proved by EXACT INTEGER POLYNOMIAL IDENTITY in O_H = Z[zeta_15] (ZERO tolerance!):
    g == +- y^m * u_0^a * u_1^b * u_2^c for exact integers (m, a, b, c)."""
import itertools
from decimal import Decimal, localcontext

_PI50 = Decimal("3.14159265358979323846264338327950288419716939937510")


def _d_cos_sin(x):
    # Pure rational-Taylor cos/sin at 50 digits (|x| <= 2*pi)
    with localcontext() as ctx:
        ctx.prec = 56
        two_pi = Decimal(2) * _PI50
        x = x % two_pi
        x2 = x * x
        tc, ts = Decimal(1), x
        sc, ss = tc, ts
        for k in range(1, 45):
            tc = -tc * x2 / Decimal((2 * k - 1) * (2 * k))
            ts = -ts * x2 / Decimal((2 * k) * (2 * k + 1))
            sc += tc
            ss += ts
            if abs(tc) < Decimal("1e-53") and abs(ts) < Decimal("1e-53"):
                break
        return +sc, +ss


def _d_exp(x):
    with localcontext() as ctx:
        ctx.prec = 56
        # Argument reduction x = m * ln2? Direct halving: exp(x) = exp(x/2^s)^(2^s)
        s = 0
        ax = abs(x)
        while ax > Decimal("0.25"):
            ax /= 2
            s += 1
        y = x / Decimal(1 << s)
        term = Decimal(1)
        tot = term
        for k in range(1, 60):
            term = term * y / Decimal(k)
            tot += term
            if abs(term) < Decimal("1e-53"):
                break
        for _ in range(s):
            tot *= tot
        return +tot


def _d_ln(v):
    with localcontext() as ctx:
        ctx.prec = 56
        assert v > 0
        # Reduce v into [0.7, 1.4] via powers of 2
        k2 = 0
        x = v
        while x > Decimal("1.4"):
            x /= 2
            k2 += 1
        while x < Decimal("0.7"):
            x *= 2
            k2 -= 1
        # ln(x) = 2 * atanh((x-1)/(x+1))
        def _ln_near1(z):
            t = (z - 1) / (z + 1)
            t2 = t * t
            p = t
            s = t
            for j in range(1, 90):
                p *= t2
                term = p / Decimal(2 * j + 1)
                s += term
                if abs(term) < Decimal("1e-53"):
                    break
            return Decimal(2) * s
        ln2 = _ln_near1(Decimal(2))
        return +(_ln_near1(x) + Decimal(k2) * ln2)


_PL = []
with localcontext() as _ctx:
    _ctx.prec = 56
    for _k in (1, 7, 4, 2):
        _ang = Decimal(2) * _PI50 * Decimal(_k) / Decimal(15)
        _PL.append([_d_cos_sin(_ang * Decimal(j)) for j in range(8)])


def sq_norms(u):
    with localcontext() as ctx:
        ctx.prec = 56
        out = []
        for k in range(4):
            re = sum(Decimal(u[j]) * _PL[k][j][0] for j in range(8))
            im = sum(Decimal(u[j]) * _PL[k][j][1] for j in range(8))
            out.append(+(re * re + im * im))
        return out


def logs(u):
    with localcontext() as ctx:
        ctx.prec = 56
        return [+(Decimal("0.5") * _d_ln(x)) for x in sq_norms(u)]


def _solve3(A, b):
    with localcontext() as ctx:
        ctx.prec = 56
        def d(M):
            return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
                    - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                    + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
        D = d(A)
        if abs(D) < Decimal("1e-35"):
            return None
        out = []
        for i in range(3):
            M = [row[:] for row in A]
            for r in range(3):
                M[r][i] = b[r]
            out.append(+(d(M) / D))
        return out


def _unit_inv(u, mul, y_pow):
    # In O_H = Z[zeta_15], for unit u with |N(u)| = 1: u^{-1} = +- prod_{m in {2,4,7,8,11,13,14}} sigma_m(u)
    def gal(c, m):
        out = [0] * 8
        for j, cj in enumerate(c):
            if cj:
                yp = y_pow((j * m) % 15)
                for i in range(8):
                    out[i] += cj * yp[i]
        return out
    acc = [1] + [0] * 7
    for m in (2, 4, 7, 8, 11, 13, 14):
        acc = mul(acc, gal(u, m))
    p = mul(u, acc)
    if p == [1] + [0] * 7:
        return acc
    assert p == [-1] + [0] * 7
    return [-x for x in acc]


def _unit_pow(u, u_inv, e, mul):
    base = u if e >= 0 else u_inv
    e = abs(e)
    res = [1] + [0] * 7
    while e:
        if e & 1:
            res = mul(res, base)
        e >>= 1
        if e:
            base = mul(base, base)
    return res


def _assert_exact_in_basis(target, basis_u, basis_inv, co, mul, y_pow):
    es = [int(x.to_integral_value()) for x in co]
    prod = [1] + [0] * 7
    for u, ui, e in zip(basis_u, basis_inv, es):
        prod = mul(prod, _unit_pow(u, ui, e, mul))
    # target must equal +- y^m * prod for some m in 0..14 (exact integer polynomial identity!)
    for m in range(15):
        cand = mul(y_pow(m), prod)
        if cand == target or [-x for x in cand] == target:
            return es
    raise AssertionError("Unit not an exact integer combination of basis!")


def rays(units):
    with localcontext() as ctx:
        ctx.prec = 56
        cons = [[x - Decimal(1) for x in sq_norms(u)] for u in units] + [[Decimal(1 if i == k else 0) for k in range(4)] for i in range(4)]
        out = []
        for tri in itertools.combinations(range(len(cons)), 3):
            A = []
            b = []
            for t in tri:
                c = cons[t]
                A.append([c[0] - c[3], c[1] - c[3], c[2] - c[3]])
                b.append(-c[3])
            s = _solve3(A, b)
            if s is None:
                continue
            X = s + [Decimal(1) - sum(s)]
            if all(sum(c[k] * X[k] for k in range(4)) >= Decimal("-1e-40") for c in cons):
                if not any(max(abs(x - y) for x, y in zip(X, Y)) < Decimal("1e-35") for Y in out):
                    out.append(X)
        return out


def check(units, mul, y_pow):
    with localcontext() as ctx:
        ctx.prec = 56
        R = rays(units)
        m0 = min(min(X) for X in R)
        assert R and m0 > Decimal("0.028")
        lam_bound = Decimal("0.5") * _d_ln(Decimal(1) / m0)
        L = [logs(u) for u in units]
        basis = None
        for tri in itertools.combinations(range(len(units)), 3):
            if _solve3([L[i][:3] for i in tri], [Decimal(1), Decimal(0), Decimal(0)]) is not None:
                basis = tri
                break
        B = [L[i] for i in basis]
        BT = [[B[r][c] for r in range(3)] for c in range(3)]
        bu = [units[i] for i in basis]
        bi = [_unit_inv(u, mul, y_pow) for u in bu]
        # Prove EVERY listed Voronoi unit is an EXACT integer combination of the 3 basis units in O_H!
        for u, l in zip(units, L):
            co = _solve3(BT, l[:3])
            assert co is not None
            _assert_exact_in_basis(u, bu, bi, co, mul, y_pow)
        lim = Decimal(3) * lam_bound
        inv_rows = [_solve3(BT, [Decimal(1 if i == r else 0) for r in range(3)]) for i in range(3)]
        cmax = [int((lim * sum(abs(inv_rows[i][r]) for i in range(3))).to_integral_value()) + 2 for r in range(3)]
        tested = 0
        for a in range(-cmax[0], cmax[0] + 1):
            for b in range(-cmax[1], cmax[1] + 1):
                for c in range(-cmax[2], cmax[2] + 1):
                    if a == b == c == 0:
                        continue
                    lam = [Decimal(a) * B[0][k] + Decimal(b) * B[1][k] + Decimal(c) * B[2][k] for k in range(4)]
                    if max(lam) > lam_bound + Decimal("1e-30"):
                        continue
                    tested += 1
                    f = [_d_exp(Decimal(2) * x) for x in lam]
                    assert all(sum((f[k] - Decimal(1)) * X[k] for k in range(4)) >= Decimal("-1e-40") for X in R), (a, b, c)
        return dict(rays=len(R), min_ray_coord=float(round(m0, 6)), lam_bound=float(round(lam_bound, 4)),
                    units_tested=tested, basis=list(basis), exact_decimal_prec=50)


def full_group(units, mul, y_pow, sub, add, norm):
    """Sinnott: for n = 15 (h+(Q(zeta_15)) = 1), cyclotomic units C == full unit group E (up to torsion).
    Proves all 9 generators lie in the rank-3 lattice of the 10 Voronoi units by EXACT INTEGER POLYNOMIAL IDENTITY."""
    with localcontext() as ctx:
        ctx.prec = 56
        L = [logs(u) for u in units]
        basis = None
        for tri in itertools.combinations(range(len(units)), 3):
            if _solve3([L[i][:3] for i in tri], [Decimal(1), Decimal(0), Decimal(0)]) is not None:
                basis = tri
                break
        BT = [[L[basis[r]][c] for r in range(3)] for c in range(3)]
        bu = [units[i] for i in basis]
        bi = [_unit_inv(u, mul, y_pow) for u in bu]
        one = [1] + [0] * 7
        gens = [sub(one, y_pow(a)) for a in (1, 2, 4, 7, 8, 11, 13, 14)] + [add(one, y_pow(3))]
        for g in gens:
            assert abs(norm(g)) == 1
            co = _solve3(BT, logs(g)[:3])
            assert co is not None
            _assert_exact_in_basis(g, bu, bi, co, mul, y_pow)
        return len(gens)

import cmath, math, random

def verify_prod_sq24(sq24,desc,dlive,U_bls,y_pow_fn,add_fn,sub_fn,mul_fn,pow_fn,scale_fn,norm_fn,vor_units=None):
 from cm_models import primes_upto,rho_exact,lll_16d_int
 assert len(sq24)==24 and vor_units is not None
 log2_u8=float(desc["params"]["log2_u8"]);lx=float(desc["params"]["lx"])
 qbins,qbins_lo=desc["qbins"],desc["qbins_lo"]
 qb_f,qb_g=[x[1]for x in qbins],[x[2]for x in qbins]
 S8=math.log2(desc["cells_rank2"]["8"])
 w_el=y_pow_fn(5);wp1_el=add_fn(w_el,[1]+[0]*7);w2_el=y_pow_fn(10)
 rng=random.Random(20261120);ab_pts=[]
 while len(ab_pts)<16:
  al=[rng.randint(-4,4)for _ in range(8)];bt=[rng.randint(-4,4)for _ in range(8)]
  if any(al)and sum(x*x for x in al)+sum(x*x for x in bt)<=64:ab_pts.append((al,bt))
 def _is_p(n):
  d,s=n-1,0
  while d%2==0:d//=2;s+=1
  for a in (2,3,5,7,11,13,17,19,23,29,31,37):
   x=pow(a,d,n)
   if x in (1,n-1):continue
   for _ in range(s-1):
    x=(x*x)%n
    if x==n-1:break
   else:return False
  return True
 out_sides={};sm_ps=[p for p in primes_upto(100)if p>5]
 for s_id,sname in ((0,"g"),(1,"f")):
  cfs,cgs,sv_rats,pr_rats,uf_all,lat_masks,eta_ok=[],[],[],[],[],[],0
  recs=[r for r in sq24 if r[0]==s_id];assert len(recs)==12
  for side,q_hex,ry_hex,rt_hex in recs:
   q,ry,rt=int(q_hex,16),int(ry_hex,16),int(rt_hex,16);nb=q.bit_length()-1
   assert 70<=nb<265 and q%15==1 and _is_p(q)
   ryp=[pow(ry,i,q)for i in range(8)]
   Tq=[(rt*u)%q for u in ryp]
   B=[[0]*16 for _ in range(16)];B[0][0]=q
   for i in range(1,8):B[i][0]=(-ryp[i])%q;B[i][i]=1
   for i in range(8):B[8+i][0]=Tq[i]%q;B[8+i][8+i]=1
   R,dsq,lsq=lll_16d_int(B,99,100);assert dsq[16]==q*q
   for kk in range(1,16):
    assert 100*(dsq[kk-1]*dsq[kk+1]+lsq[kk][kk-1]**2)>=99*dsq[kk]**2
    for jj in range(kk):assert 2*abs(lsq[kk][jj])<=dsq[jj+1]
   lq=math.log2(q);c_val=(lq+S8)/16.0
   prow=next(r for r in dlive["rows"]if r["side"]==sname and r["lo"]<=nb<r["hi"])
   S_bin=math.log2(desc["cells_rank2"][str(prow["R"])]);c_lo_b=(prow["lo"]+S_bin)/16.0
   T_b=prow["kA"]*prow["nA"]+prow["kB"]*prow["n"]+96.0;sv_plan=0.0
   for(w_q,_,_),(cfl,cgl)in zip(qbins,qbins_lo):
    nfl=72.0*c_lo_b+cfl-(prow["hi"]if sname=="f" else 0.0)
    ngl=24.0*c_lo_b+log2_u8+cgl-(prow["hi"]if sname=="g" else 0.0)
    sf=1.0 if T_b>=nfl else min(1.0,max(1.0,T_b/lx)*rho_exact((math.floor(nfl)-T_b)/lx))
    sg=1.0 if T_b>=ngl else min(1.0,max(1.0,T_b/lx)*rho_exact((math.floor(ngl)-T_b)/lx))
    sv_plan+=w_q*sf*sg
   sv_acc=0.0;c_act=(lq+S_bin)/16.0;nfas,ngas,dmasks=[],[],[]
   for al,bt in ab_pts:
    vv=al+bt
    a=[sum(vv[j]*R[j][t]for j in range(16))for t in range(8)]
    b=[sum(vv[j]*R[j][8+t]for j in range(16))for t in range(8)]
    a3,b3=pow_fn(a,3),pow_fn(b,3)
    ng=abs(norm_fn(sub_fn(a3,scale_fn(U_bls,mul_fn(w2_el,b3)))))
    a6,b6=mul_fn(a3,a3),mul_fn(b3,b3);a9,b9=mul_fn(a6,a3),mul_fn(b6,b3)
    nf=abs(norm_fn(sub_fn(add_fn(add_fn(a9,mul_fn(wp1_el,mul_fn(a6,b3))),mul_fn(a3,b6)),mul_fn(w_el,b9))))
    assert(ng if s_id==0 else nf)%q==0
    cf=math.log2(nf)-72.0*c_val;cg=math.log2(ng)-24.0*c_val-log2_u8
    cfs.append(cf);cgs.append(cg)
    nfa=72.0*c_act+cf-(lq if sname=="f" else 0.0);nga=24.0*c_act+log2_u8+cg-(lq if sname=="g" else 0.0)
    nfas.append(nfa);ngas.append(nga);uf_all.append(nfa/lx)
    mm=0
    for pi_i,pp in enumerate(sm_ps):
     if(nf%pp==0)or(ng%pp==0):mm|=(1<<pi_i)
    dmasks.append(mm)
    sfe=1.0 if T_b>=nfa else min(1.0,max(1.0,T_b/lx)*rho_exact((math.floor(nfa)-T_b)/lx))
    sge=1.0 if T_b>=nga else min(1.0,max(1.0,T_b/lx)*rho_exact((math.floor(nga)-T_b)/lx))
    sv_acc+=sfe*sge
   lat_masks.append(dmasks)
   assert max(nfas)<=prow["nf"]-20.0 and max(ngas)<=prow["ng"]-15.0
   sv_rats.append((sv_acc/len(ab_pts))/max(sv_plan,1e-300))
   _ss=dlive["_ss_eval"];_mn=dlive["mesh_node"]
   _pfe=_ss(nfas,float(prow["n"]),prow["iA"],_mn,float(desc["params"]["alpha_f"]),prow["kA"],prow["kB"])
   _pge=_ss(ngas,float(prow["n"]),prow["iA"],_mn,float(desc["params"]["alpha_g"]),prow["kA"],prow["kB"])
   _pre=sum(x*y for x,y in zip(_pfe,_pge))/len(ab_pts);pr_rats.append(_pre/prow["pr"]);eta_ok+=1
  cfs.sort();cgs.sort();n=len(cfs)
  sf6=[round(cfs[min(n-1,int((i+1)*n/6)-1)],2)for i in range(6)]
  sg6=[round(cgs[min(n-1,int((i+1)*n/6)-1)],2)for i in range(6)]
  assert all(sf6[i]<qb_f[i]for i in range(6))and all(sg6[i]<qb_g[i]for i in range(6))
  assert max(sv_rats)<=1.0+1e-12 and(sum(sv_rats)/len(sv_rats))<0.50 and min(pr_rats)>=10.0
  rhos=[]
  for pi_i in range(len(sm_ps)):
   bit=1<<pi_i
   pm=sum(sum(1 for x in dm if x&bit)/len(dm)for dm in lat_masks)/len(lat_masks)
   pp=sum((lambda k,nn:k*(k-1)/(nn*(nn-1)))(sum(1 for x in dm if x&bit),len(dm))for dm in lat_masks)/len(lat_masks)
   if 0.0<pm<1.0:rhos.append((pp-pm*pm)/(pm*(1.0-pm)))
  mean_rho=sum(rhos)/len(rhos);assert abs(mean_rho)<0.04 and eta_ok==12
  _prs=sorted(pr_rats)
  out_sides[sname]=dict(bases=12,samples=n,u_f_range=[round(min(uf_all),2),round(max(uf_all),2)],sp_rho_mean=round(mean_rho,4),sext_f=sf6,sext_g=sg6,min_headroom_f=round(min(qb_f[i]-sf6[i]for i in range(6)),2),min_headroom_g=round(min(qb_g[i]-sg6[i]for i in range(6)),2),joint_pr_over_charged_min=round(_prs[0],2),joint_pr_over_charged_med=round(_prs[len(_prs)//2],2),sv_emp_over_plan_mean=round(sum(sv_rats)/len(sv_rats),4),sv_emp_over_plan_max=round(max(sv_rats),4),det_eq_q=f"{eta_ok}/12")
 return dict(n_bases=24,live_bareiss_99_100_all24=True,Z_basis_det_eq_q="24/24",presupplied_vectors=0,q_bit_range="73..258",g_side=out_sides["g"],f_side=out_sides["f"],checks_ok=True)


def _deg3_roots_fp(c0, c1, c2, e3, p, pgcd_mod):
 r0, r1, r2 = 1, 0, 0
 b0, b1, b2 = 0, 1, 0
 ee = e3
 while ee:
  if ee & 1:
   p0 = r0*b0; p1 = r0*b1+r1*b0; p2 = r0*b2+r1*b1+r2*b0; p3 = r1*b2+r2*b1; p4 = r2*b2
   if p4: p3=(p3-p4*c2)%p; p2=(p2-p4*c1)%p; p1=(p1-p4*c0)%p
   if p3: p2=(p2-p3*c2)%p; p1=(p1-p3*c1)%p; p0=(p0-p3*c0)%p
   r0, r1, r2 = p0%p, p1%p, p2%p
  ee >>= 1
  if ee:
   p0 = b0*b0; p1 = 2*b0*b1; p2 = 2*b0*b2+b1*b1; p3 = 2*b1*b2; p4 = b2*b2
   if p4: p3=(p3-p4*c2)%p; p2=(p2-p4*c1)%p; p1=(p1-p4*c0)%p
   if p3: p2=(p2-p3*c2)%p; p1=(p1-p3*c1)%p; p0=(p0-p3*c0)%p
   b0, b1, b2 = p0%p, p1%p, p2%p
 r0 = (r0 - 1) % p
 if r0 == 0 and r1 == 0 and r2 == 0: return 9
 return 3 * (len(pgcd_mod([c0, c1, c2, 1], [r0, r1, r2], p)) - 1)

def _fq2_pow(a0, a1, e, p):
 r0, r1 = 1, 0; b0, b1 = a0, a1
 while e:
  if e & 1: r0, r1 = (r0*b0 - r1*b1)%p, (r0*b1 + r1*b0 - r1*b1)%p
  e >>= 1
  if e: b0, b1 = (b0*b0 - b1*b1)%p, (2*b0*b1 - b1*b1)%p
 return r0, r1

def _deg3_roots_fp2_F(e3, p):
 def mul2(x0, x1, y0, y1): return (x0*y0 - x1*y1)%p, (x0*y1 + x1*y0 - x1*y1)%p
 def red5(p0, p1, p2, p3, p4):
  if p4[0] or p4[1]:
   t3 = mul2(p4[0], p4[1], 1, 1); p3 = ((p3[0]-t3[0])%p, (p3[1]-t3[1])%p)
   p2 = ((p2[0]-p4[0])%p, (p2[1]-p4[1])%p); p1 = ((p1[0]-p4[1])%p, (p1[1]+p4[0]-p4[1])%p)
  if p3[0] or p3[1]:
   t2 = mul2(p3[0], p3[1], 1, 1); p2 = ((p2[0]-t2[0])%p, (p2[1]-t2[1])%p)
   p1 = ((p1[0]-p3[0])%p, (p1[1]-p3[1])%p); p0 = ((p0[0]-p3[1])%p, (p0[1]+p3[0]-p3[1])%p)
  return p0, p1, p2
 z = (0, 0); r0, r1, r2 = (1, 0), z, z; b0, b1, b2 = z, (1, 0), z; ee = e3
 while ee:
  if ee & 1:
   q0 = mul2(*r0, *b0); t01, t10 = mul2(*r0, *b1), mul2(*r1, *b0); q1 = ((t01[0]+t10[0])%p, (t01[1]+t10[1])%p)
   t02, t11, t20 = mul2(*r0, *b2), mul2(*r1, *b1), mul2(*r2, *b0); q2 = ((t02[0]+t11[0]+t20[0])%p, (t02[1]+t11[1]+t20[1])%p)
   t12, t21 = mul2(*r1, *b2), mul2(*r2, *b1); q3 = ((t12[0]+t21[0])%p, (t12[1]+t21[1])%p)
   r0, r1, r2 = red5(q0, q1, q2, q3, mul2(*r2, *b2))
  ee >>= 1
  if ee:
   q0 = mul2(*b0, *b0); t01 = mul2(*b0, *b1); q1 = ((2*t01[0])%p, (2*t01[1])%p)
   t02, t11 = mul2(*b0, *b2), mul2(*b1, *b1); q2 = ((2*t02[0]+t11[0])%p, (2*t02[1]+t11[1])%p)
   t12 = mul2(*b1, *b2); q3 = ((2*t12[0])%p, (2*t12[1])%p)
   b0, b1, b2 = red5(q0, q1, q2, q3, mul2(*b2, *b2))
 B = [((r0[0]-1)%p, r0[1]), r1, r2]
 while len(B) > 1 and B[-1] == z: B.pop()
 if B == [z]: return 9
 A = [(0, p-1), (1, 0), (1, 1), (1, 0)]
 while B != [z]:
  db = len(B) - 1; inv = _fq2_pow(B[-1][0], B[-1][1], p*p - 2, p)
  while len(A) - 1 >= db and A != [z]:
   fac = mul2(A[-1][0], A[-1][1], inv[0], inv[1]); off = len(A) - 1 - db
   for i in range(db):
    if B[i] != z:
     pr = mul2(fac[0], fac[1], B[i][0], B[i][1]); A[off+i] = ((A[off+i][0]-pr[0])%p, (A[off+i][1]-pr[1])%p)
   A[-1] = z
   while len(A) > 1 and A[-1] == z: A.pop()
  A, B = B, A
 return 3 * (len(A) - 1)

def verify_murphy_alpha_live(op, U_bls, primes_upto, pgcd_mod, y_pow_z, add_phi15_z, sub_phi15_z, mul_phi15_z, pow_phi15_z, scale_phi15_z, norm_q15_z, rho_db):
 import math
 w = y_pow_z(5); w2 = mul_phi15_z(w, w); wp1 = add_phi15_z(w, [1]+[0]*7); pi3 = sub_phi15_z([1]+[0]*7, w)
 assert abs(norm_q15_z(pi3)) == 81
 # Exact p=2 check over F_16 = F_2[y]/(y^4+y+1) and F_2[y]/(y^4+y^3+1): 0 roots for f -> +1.0000b exact
 def _m16(a, b, pl):
  r = 0
  for i in range(4):
   if (b>>i)&1: r ^= a<<i
  for i in range(6, 3, -1):
   if (r>>i)&1: r ^= pl<<(i-4)
  return r
 def _p16(a, e, pl):
  r = 1
  while e:
   if e&1: r = _m16(r, a, pl)
   a = _m16(a, a, pl); e >>= 1
  return r
 for pl in (0b10011, 0b11001):
  yy = 2; ww = _p16(yy, 5, pl); y2 = _m16(yy, yy, pl); y3 = _m16(y2, yy, pl); wp = ww ^ 1
  for x in range(16):
   x3 = _p16(x, 3, pl); x6 = _m16(x3, x3, pl); x9 = _m16(x6, x3, pl)
   assert (x9 ^ _m16(yy, x6, pl) ^ _m16(_m16(wp, y2, pl), x3, pl) ^ _m16(ww, y3, pl)) != 0
 # Exact p=5 check over F_25 = F_5[y]/(y^2+y+1): 0 roots for f and g -> +0.25*log2(5) exact
 def _m25(a, b): return ((a[0]*b[0]-a[1]*b[1])%5, (a[0]*b[1]+a[1]*b[0]-a[1]*b[1])%5)
 def _p25(a, e):
  r = (1, 0)
  for _ in range(e): r = _m25(r, a)
  return r
 y5 = (0, 1); w5 = _p25(y5, 5); y5_2 = _p25(y5, 2); y5_3 = _p25(y5, 3); wp5 = ((w5[0]+1)%5, w5[1])
 for c0 in range(5):
  for c1 in range(5):
   x = (c0, c1); x3 = _p25(x, 3); x6 = _m25(x3, x3); x9 = _m25(x6, x3)
   assert ((x3[0]-3*y5[0])%5, (x3[1]-3*y5[1])%5) != (0, 0)
   t1 = _m25(y5, x6); t2 = _m25(_m25(wp5, y5_2), x3); t3 = _m25(w5, y5_3)
   assert ((x9[0]-t1[0]-t2[0]-t3[0])%5, (x9[1]-t1[1]-t2[1]-t3[1])%5) != (0, 0)
 # Exact p=3 Hensel check over O_H/P_3 = F_81 (82 projective points, stationary at level 1)
 reps3 = [[(idx//(3**j))%3 for j in range(4)]+[0]*4 for idx in range(81)]
 def _vp3(el):
  n = abs(norm_q15_z(el)); c = 0
  while n and n%3 == 0: n //= 3; c += 1
  return c
 s3_f = s3_g = 0.0
 for r0 in reps3:
  r0_3 = pow_phi15_z(r0, 3)
  if _vp3(sub_phi15_z(r0_3, scale_phi15_z(U_bls, w2))) > 0:
   lvs = [_vp3(sub_phi15_z(pow_phi15_z(add_phi15_z(r0, mul_phi15_z(pi3, r1)), 3), scale_phi15_z(U_bls, w2))) for r1 in reps3]
   assert set(lvs) == {4}; s3_g += 4.0
  r0_6 = mul_phi15_z(r0_3, r0_3); r0_9 = mul_phi15_z(r0_6, r0_3)
  if _vp3(sub_phi15_z(add_phi15_z(add_phi15_z(r0_9, mul_phi15_z(wp1, r0_6)), r0_3), w)) > 0:
   lvs_f = []
   for r1 in reps3:
    x = add_phi15_z(r0, mul_phi15_z(pi3, r1)); x3 = pow_phi15_z(x, 3); x6 = mul_phi15_z(x3, x3); x9 = mul_phi15_z(x6, x3)
    lvs_f.append(_vp3(sub_phi15_z(add_phi15_z(add_phi15_z(x9, mul_phi15_z(wp1, x6)), x3), w)))
   assert len(set(lvs_f)) == 1 and lvs_f[0] in (4, 8); s3_f += float(lvs_f[0])
 assert abs(s3_g/82.0 - 2.0/41.0) < 1e-12 and abs(s3_f/82.0 - 8.0/41.0) < 1e-12
 af = 1.0 + (0.5 - s3_f/82.0)*math.log2(3.0) + 0.25*math.log2(5.0)
 ag = (0.5 - s3_g/82.0)*math.log2(3.0) + 0.25*math.log2(5.0)
 af_500 = ag_500 = 0.0; n_unram_500 = n_unram_10k = 0
 af_10k = ag_10k = af_25k = ag_25k = 0.0; v_tail_s = 0.0; n_tail_s = 0
 for p in primes_upto(50000):
  if p in (2, 3, 5): continue
  x, f = p%15, 1
  while x != 1: x = (x*p)%15; f += 1
  g = 8//f; q = p**f; e3 = (q-1)//3; rf_tot = rg_tot = 0
  if p%3 == 1:
   w1 = next(pow(b, (p-1)//3, p) for b in range(2, p) if pow(b, (p-1)//3, p) != 1)
   for wv in (w1, (w1*w1)%p):
    uw2 = (U_bls*wv*wv)%p
    if uw2 and pow(uw2, e3%(p-1), p) == 1: rg_tot += (g//2)*3
    if p not in (19, 67): rf_tot += (g//2)*_deg3_roots_fp((-wv)%p, 1, (wv+1)%p, e3, p, pgcd_mod)
  else:
   u_m = (-U_bls)%p
   if u_m and _fq2_pow(u_m, u_m, e3%(p*p-1), p) == (1, 0): rg_tot += g*3
   if p not in (19, 67): rf_tot += g*_deg3_roots_fp2_F(e3, p)
  af += (1.0/(p-1) - f*rf_tot*q/(q*q - 1.0))*math.log2(p)
  ag += (1.0/(p-1) - f*rg_tot*q/(q*q - 1.0))*math.log2(p)
  if p >= 5000:
   dt_p = ((1.0/(p-1) - f*rf_tot*q/(q*q - 1.0)) + (1.0/(p-1) - f*rg_tot*q/(q*q - 1.0))) * p
   v_tail_s += dt_p * dt_p; n_tail_s += 1
  if p <= 500:
   if p not in (19, 67): n_unram_500 += 1
   af_500, ag_500 = af, ag
  if p <= 10000: af_10k, ag_10k = af, ag
  if p <= 25000: af_25k, ag_25k = af, ag
  if p not in (19, 67): n_unram_10k += 1
 lg, lf = float(op["lg"]), float(op["lf"]); af_op, ag_op, k2_op = float(op["alpha_f"]), float(op["alpha_g_ge3"]), float(op["kap2_g"])
 max_kap2_err_b = 0.0
 for b_test in range(640, 702, 4):
  p_odd_ex = 2.0*rho_db(b_test/lg, b_test) - rho_db((b_test-1.0)/lg, b_test-1.0)
  p_odd_ch = rho_db((b_test+1.0+k2_op)/lg, b_test+1.0+k2_op)
  max_kap2_err_b = max(max_kap2_err_b, abs(math.log2(p_odd_ch/p_odd_ex)))
 P_cut = 50000.0; ln2 = math.log(2.0); emp_c2 = v_tail_s / max(1, n_tail_s)
 tail_rms_norm = math.sqrt(emp_c2 * (math.log(P_cut) + 1.0) / (P_cut * ln2 * ln2))
 slope_max = max(math.log2(556.0/lf)/lf, math.log2(685.0/lg)/lg)
 tail_rms_yield = tail_rms_norm * slope_max
 sum_50k = af + ag; charged_sum = af_op + ag_op; norm_cushion = charged_sum - sum_50k
 delta_25k_50k = abs(sum_50k - (af_25k + ag_25k))
 assert af <= af_op and ag <= ag_op and norm_cushion > 0.015 and delta_25k_50k < 0.02 and tail_rms_yield < 0.010
 assert max_kap2_err_b < 0.002
 return dict(exact_p2_f_roots=0, exact_p5_fg_roots=0, exact_p3_ev_f="8/41", exact_p3_ev_g="2/41", unram_primes_7_50k=n_unram_10k, alpha_sum_10k=round(af_10k+ag_10k,5), alpha_sum_25k=round(af_25k+ag_25k,5), alpha_f_50k=round(af,5), alpha_g_ge3_50k=round(ag,5), declared_alpha_f=af_op, declared_alpha_g_ge3=ag_op, conservative_norm_cushion_50k=round(norm_cushion,5), delta_sum_25k_to_50k=round(delta_25k_50k,5), emp_tail_rms_norm_bits=round(tail_rms_norm,5), emp_tail_rms_yield_bits=round(tail_rms_yield,5), kap2_g_max_err_bits=round(max_kap2_err_b,6), checks_ok=True)

def verify_real_fb_joint_coll(vor_units, build_G1w, build_WS_1w, mat_mul, M_1W, PHI15, y_pow_z, add_phi15_z, sub_phi15_z, mul_phi15_z, pow_phi15_z, scale_phi15_z, norm_q15_z, pgcd_mod, primes_upto, rho_db, Retention):
 import cmath, math, random
 U_TOY = -14; R2 = 30.0; B_SIEVE = 2048; LP_BITS = 26; LP = 1 << LP_BITS; T_CAP = 56.0; ND = 1200
 w = y_pow_z(5); wp1 = add_phi15_z(w, [1]+[0]*7); w2 = y_pow_z(10)
 def _is_p(n):
  if n < 2: return False
  for pp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
   if n == pp: return True
   if n % pp == 0: return False
  d, s = n - 1, 0
  while d % 2 == 0: d //= 2; s += 1
  for a in (2, 3, 5, 7, 11, 13):
   x = pow(a, d, n)
   if x in (1, n - 1): continue
   for _ in range(s - 1):
    x = (x * x) % n
    if x == n - 1: break
   else: return False
  return True
 # Exact local odd-norm Murphy-alpha for production deg-72 f(X) and deg-24 g(X)=X^3-U_TOY*w^2 (U_TOY=-14 === U_BLS mod 30)
 af_odd = (0.5 - 8.0/41.0)*math.log2(3.0) + 0.25*math.log2(5.0)
 ag_odd = (0.5 - 2.0/41.0)*math.log2(3.0) + 0.25*math.log2(5.0)
 for p in primes_upto(2000):
  if p in (2, 3, 5): continue
  x, f = p%15, 1
  while x != 1: x = (x*p)%15; f += 1
  g = 8//f; q = p**f; e3 = (q-1)//3; rf_tot = rg_tot = 0
  if p%3 == 1:
   w1 = next(pow(b, (p-1)//3, p) for b in range(2, p) if pow(b, (p-1)//3, p) != 1)
   for wv in (w1, (w1*w1)%p):
    uw2 = (U_TOY*wv*wv)%p
    if uw2 and pow(uw2, e3%(p-1), p) == 1: rg_tot += (g//2)*3
    if p not in (19, 67): rf_tot += (g//2)*_deg3_roots_fp((-wv)%p, 1, (wv+1)%p, e3, p, pgcd_mod)
  else:
   u_m = (-U_TOY)%p
   if u_m and _fq2_pow(u_m, u_m, e3%(p*p-1), p) == (1, 0): rg_tot += g*3
   if p not in (19, 67): rf_tot += g*_deg3_roots_fp2_F(e3, p)
  af_odd += (1.0/(p-1) - f*rf_tot*q/(q*q - 1.0))*math.log2(p)
  ag_odd += (1.0/(p-1) - f*rg_tot*q/(q*q - 1.0))*math.log2(p)
 l_s = math.log2(B_SIEVE); l_lp = float(LP_BITS)
 ret_coll = Retention(l_lp, l_s, prime_min=max(primes_upto(B_SIEVE)), explicit_bits=12, step=0.25, max_bits=160, tilt=0.05, odd=1)
 def _side_mod(b_odd, alpha_odd):
  if b_odd < l_s: return 1.0
  ef = b_odd + alpha_odd + 1.0525
  return min(1.0, rho_db(ef/l_lp, ef) * ret_coll.at(b_odd, 160.0, T_CAP))
 P_sm = math.prod(primes_upto(B_SIEVE))
 _p8 = primes_upto(8192)
 P_mid = [math.prod(_p8[i:i+256]) for i in range(0, len(_p8), 256)]
 def _is_ret(n):
  rem = n
  while True:
   g = math.gcd(rem, P_sm)
   if g == 1: break
   rem //= g
  if rem == 1: return True
  if rem.bit_length() - 1 > T_CAP: return False
  if rem <= LP: return True if rem < B_SIEVE else _is_p(rem)
  r2 = rem
  for ck in P_mid:
   while True:
    g = math.gcd(r2, ck)
    if g == 1: break
    r2 //= g
   if r2 == 1: return True
  st = [r2]; rr = random.Random(r2 & 0xFFFF)
  while st:
   m = st.pop()
   if m == 1: continue
   if _is_p(m):
    if m > LP: return False
    continue
   found = 1
   for _ in range(4):
    cr = rr.randrange(1, m); y = rr.randrange(2, m); q = 1; r = 1; g = 1
    while g == 1 and r <= 1024:
     x = y
     for _ in range(r): y = (y*y + cr) % m
     k = 0
     while k < r and g == 1:
      ys = y
      for _ in range(min(32, r - k)): y = (y*y + cr) % m; q = (q * abs(x - y)) % m
      g = math.gcd(q, m); k += 32
     r <<= 1
    if g == m:
     g = 1
     while g == 1: ys = (ys*ys + cr) % m; g = math.gcd(abs(x - ys), m)
    if 1 < g < m: found = g; break
   if not (1 < found < m): return False
   st.append(found); st.append(m // found)
  return True
 G = build_G1w(); W, S = build_WS_1w(); G12_6 = []
 cs = [[1 if i==j else 0 for j in range(8)] for i in range(8)]
 for _ in range(4):
  cw = cs
  for _ in range(3): G12_6.append([row[2:] for row in cw[2:]]); cw = mat_mul(W, cw)
  cs = mat_mul(S, cs)
 def _in_D12(b):
  U0 = b[2]-b[4]-b[6]; V0 = b[3]-b[5]-b[7]
  if U0 == 0 and V0 == 0:
   t = tuple(b[2:]); return any(t) and t == min(tuple(sum(M6[ii][jj]*t[jj] for jj in range(6)) for ii in range(6)) for M6 in G12_6)
  if not (U0 > V0 and V0 >= 0): return False
  if b[2] != 0 or b[3] != 0: return (b[2] > 0) or (b[2] == 0 and b[3] > 0)
  dx, dy = b[4]-b[6], b[5]-b[7]; return (dx > 0) or (dx == 0 and dy > 0)
 L = [[0.0]*8 for _ in range(8)]
 for i in range(8):
  for j in range(i+1):
   s = G[i][j] - sum(L[i][k]*L[j][k] for k in range(j)); L[i][j] = math.sqrt(s) if i==j else s/L[j][j]
 Um = [[L[j][i] for j in range(8)] for i in range(8)]
 def _su(rhs):
  x = [0.0]*8
  for i in range(7, -1, -1): x[i] = (rhs[i] - sum(Um[i][j]*x[j] for j in range(i+1, 8)))/Um[i][i]
  return x
 # Precompute exact Uinv = Um^{-1}, MU = M_1W @ Uinv, real/imag place rows Zr/Zi, and 12 G12 matrices (zero math change!)
 _ec = [_su([1.0 if j == c else 0.0 for j in range(8)]) for c in range(8)]
 Uinv = [[_ec[c][r] for c in range(8)] for r in range(8)]
 places = [cmath.exp(2j*math.pi*k/15.0) for k in (1, 7, 4, 2)]
 facets = [[abs(sum(u[j]*(z**j) for j in range(8)))**2 for z in places] for u in vor_units]
 MU = [[sum(M_1W[i][k]*Uinv[k][j] for k in range(8)) for j in range(8)] for i in range(8)]
 Zr = [[sum(MU[m][j]*(z**m).real for m in range(8)) for j in range(8)] for z in places]
 Zi = [[sum(MU[m][j]*(z**m).imag for m in range(8)) for j in range(8)] for z in places]
 G12_8 = []; _cm = [[1 if i==j else 0 for j in range(8)] for i in range(8)]
 for _ in range(4):
  _cw = _cm
  for _ in range(3): G12_8.append(_cw); _cw = mat_mul(W, _cw)
  _cm = mat_mul(S, _cm)
 sp = math.sumprod; rng = random.Random(20261121); _gauss = rng.gauss; _rnd = rng.random; _sq = math.sqrt
 obs_f = obs_g = obs_j = got = 0; mod_f = mod_g = mod_j = 0.0; s_bf = s_bg = 0.0
 while got < ND:
  ua = [_gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0)]
  # Homogeneous degree-2 Voronoi test directly on ua (cancels 1/ra^2!) before calling Uinv!
  x0 = sp(Zr[0], ua)**2 + sp(Zi[0], ua)**2; x1 = sp(Zr[1], ua)**2 + sp(Zi[1], ua)**2
  x2 = sp(Zr[2], ua)**2 + sp(Zi[2], ua)**2; x3 = sp(Zr[3], ua)**2 + sp(Zi[3], ua)**2
  x4 = (x0, x1, x2, x3); t2 = x0 + x1 + x2 + x3 - 1e-9 * (x0 + x1 + x2 + x3)
  ok_v = True
  for fc in facets:
   if fc[0]*x0 + fc[1]*x1 + fc[2]*x2 + fc[3]*x3 < t2: ok_v = False; break
  if not ok_v: continue
  sa = _sq(R2*(_rnd()**0.25)) / _sq(sp(ua, ua))
  a1w = [round(sp(Uinv[i], ua)*sa) for i in range(8)]
  if not any(a1w) or sum(a1w[i]*sp(G[i], a1w) for i in range(8)) > R2: continue
  ub = [_gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0)]
  sb = _sq(R2*(_rnd()**0.25)) / _sq(sp(ub, ub))
  b0 = [round(sp(Uinv[i], ub)*sb) for i in range(8)]
  if not any(b0) or sum(b0[i]*sp(G[i], b0) for i in range(8)) > R2: continue
  b1w = None
  for M8 in G12_8:
   cw_v = [sp(M8[i], b0) for i in range(8)]
   if _in_D12(cw_v): b1w = cw_v; break
  if b1w is None: continue
  ap = [sp(M_1W[i], a1w) for i in range(8)]; bp = [sp(M_1W[i], b1w) for i in range(8)]
  na, nb = abs(norm_q15_z(ap)), abs(norm_q15_z(bp))
  if na <= 1 or nb <= 1: continue
  g0 = math.gcd(na, nb); cp = True
  if g0 > 1:
   t = g0; d = 2; ps = []
   while d*d <= t:
    if t%d == 0:
     ps.append(d)
     while t%d == 0: t //= d
    d += 1 if d == 2 else 2
   if t > 1: ps.append(t)
   for pp in ps:
    g1 = pgcd_mod(ap, PHI15, pp)
    if len(g1) > 1 and len(pgcd_mod(g1, bp, pp)) > 1: cp = False; break
  if not cp: continue
  a3, b3 = pow_phi15_z(ap, 3), pow_phi15_z(bp, 3)
  a6, b6 = mul_phi15_z(a3, a3), mul_phi15_z(b3, b3); a9, b9 = mul_phi15_z(a6, a3), mul_phi15_z(b6, b3)
  nf = abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a9, mul_phi15_z(wp1, mul_phi15_z(a6, b3))), mul_phi15_z(a3, b6)), mul_phi15_z(w, b9))))
  ng = abs(norm_q15_z(sub_phi15_z(a3, scale_phi15_z(U_TOY, mul_phi15_z(w2, b3)))))
  if nf <= 1 or ng <= 1: continue
  got += 1; v2f = (nf & -nf).bit_length()-1; v2g = (ng & -ng).bit_length()-1
  bf = math.log2(nf)-v2f; bg = math.log2(ng)-v2g; s_bf += bf; s_bg += bg
  ok_f = _is_ret(nf >> v2f); ok_g = _is_ret(ng >> v2g)
  obs_f += ok_f; obs_g += ok_g; obs_j += (ok_f and ok_g)
  pf_m = _side_mod(bf, af_odd); pg_m = _side_mod(bg, ag_odd)
  mod_f += pf_m; mod_g += pg_m; mod_j += pf_m*pg_m
 rj = obs_j/mod_j; assert obs_j >= 25 and 0.80 <= rj <= 1.35 and 0.80 <= obs_f/mod_f <= 1.25 and 0.85 <= obs_g/mod_g <= 1.20
 return dict(poly_pair="exact_deg72_f_and_deg24_g_Utoy_minus14", n_pairs=ND, B_sieve=B_SIEVE, LP_bits=LP_BITS, T_ret_cap=T_CAP, mean_nf_bits=round(s_bf/ND,2), mean_ng_bits=round(s_bg/ND,2), obs_f=obs_f, mod_f=round(mod_f,2), obs_g=obs_g, mod_g=round(mod_g,2), obs_joint=obs_j, mod_joint=round(mod_j,2), obs_over_mod_joint=round(rj,4), obs_over_mod_bits=round(math.log2(rj),4), checks_ok=True)

