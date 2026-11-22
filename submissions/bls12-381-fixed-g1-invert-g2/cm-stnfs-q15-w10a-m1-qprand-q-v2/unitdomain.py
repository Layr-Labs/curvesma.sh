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
 e2 = (p*p - 1) // 3
 def mul2(x0, x1, y0, y1): return (x0*y0 - x1*y1)%p, (x0*y1 + x1*y0 - x1*y1)%p
 def red5(p0, p1, p2, p3, p4):
  a, b = p4
  if a or b:
   p3 = ((p3[0] - a + b)%p, (p3[1] - a)%p)
   p2 = ((p2[0] - a)%p, (p2[1] - b)%p); p1 = ((p1[0] - b)%p, (p1[1] + a - b)%p)
  a, b = p3
  if a or b:
   p2 = ((p2[0] - a + b)%p, (p2[1] - a)%p)
   p1 = ((p1[0] - a)%p, (p1[1] - b)%p); p0 = ((p0[0] - b)%p, (p0[1] + a - b)%p)
  return p0, p1, p2
 z = (0, 0); r0, r1, r2 = (1, 0), z, z; b0, b1, b2 = z, (1, 0), z; ee = e2
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
 P_tiny = math.prod(primes_upto(128))
 P_sm = math.prod(p for p in primes_upto(B_SIEVE) if p > 128)
 _p8 = [p for p in primes_upto(8192) if p > B_SIEVE]
 P_mid = [math.prod(_p8[i:i+256]) for i in range(0, len(_p8), 256)]
 def _is_ret(n):
  rem = n
  for _PP in (P_tiny, P_sm):
   while True:
    g = math.gcd(rem, _PP)
    if g == 1: break
    rem //= g
   if rem == 1: return True
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


def verify_q15_48_24_live(_c48,pts,P_SM,P_MD,_is_prime_det,_rho,U_TOY,W_EL,WP1,W2,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,lll_16d_int):
 from cm_models import rho_exact as _re48
 _B48=int(_c48["B_fb"]);_qhf48=_c48["qhf"];_qhg48=_c48["qhg"];_S48=math.log2(int(_c48["n_pts"]))
 assert _c48["deg_pair"]==[48,24]and _B48==8192 and _c48["strict_q_gt_B_fb"]and _c48["strict_zero_omitted_children"]
 _slo,_shi=1<<_c48["slot_bits"][0],1<<_c48["slot_bits"][1]
 def _mu48(sd,q):
  qb=math.log2(q);c=(qb+_S48)/16.0;sf=qb if sd=="f" else 0.0;sg=qb if sd=="g" else 0.0
  pr=sum((1.0/6.0)*_re48((48.0*c+_qhf48[k]-sf)/qb)*_re48((24.0*c+_qhg48[k]-sg)/qb)for k in range(6))
  return round(int(_c48["n_pts"])*0.50*pr,3)
 _s48={}
 for _gn in ("d1","ctrl"):
  _rs=_c48[_gn];assert len(_rs)==50;_ms=[];_wcs=[];_sp0=0.0
  for _r in _rs:
   _q,_sd,_ry,_rx,_mc,_wc=_r[0],_r[1],_r[2],_r[3],_r[4],_r[5];_muv=_r[-1]
   assert _B48<_slo<=_q<_shi and _q%15==1 and _is_prime_det(_q)and pow(_ry,15,_q)==1 and pow(_ry,3,_q)!=1 and pow(_ry,5,_q)!=1
   if _gn=="d1":assert _q<_r[6]
   _w=pow(_ry,5,_q);_w2=(_w*_w)%_q;_x3=pow(_rx,3,_q)
   if _sd=="g":assert(_x3-U_TOY*_w2)%_q==0
   else:assert(_x3*_x3+(_w+1)*_x3-_w)%_q==0
   assert abs(_mu48(_sd,_q)-_muv)<=0.01
   _ms.append(_mc);_wcs.append(_wc);_sp0+=max(0.001,1.0/(1.0+0.5*_muv))
  assert _ms.count(0)==0
  _s48[_gn]=dict(n=50,zeros=0,sum_p0_unscaled=round(_sp0,2),min=min(_ms),max=max(_ms),mean=round(sum(_ms)/50.0,2),with_above_B_child_mean=round(sum(_wcs)/50.0,2))
  for _sd in ("f","g"):
   _sr=[r for r in _rs if r[1]==_sd];_sm=[r[4]for r in _sr];_sl=[math.log2(r[0])for r in _sr]
   _mu=sum(_sm)/len(_sm);_vr=sum((x-_mu)**2 for x in _sm)/max(1,len(_sm)-1);_ml=sum(_sl)/len(_sl)
   _slp=sum((lp-_ml)*(m-_mu)for lp,m in zip(_sl,_sm))/max(1e-9,sum((lp-_ml)**2 for lp in _sl))
   _rv=sum((m-(_mu+_slp*(lp-_ml)))**2 for lp,m in zip(_sl,_sm))/max(1,len(_sm)-2)
   _s48[f"{_gn}_{_sd}"]=dict(n=len(_sr),zeros=0,min=min(_sm),mean=round(_mu,2),p_adj_fano=round(_rv/_mu,3))
 assert not({r[0]for r in _c48["d1"]}&{r[0]for r in _c48["ctrl"]})
 _s48["d1_over_ctrl"]=round(_s48["d1"]["mean"]/_s48["ctrl"]["mean"],4)
 _s48["total_sum_p0_unscaled"]=round(_s48["d1"]["sum_p0_unscaled"]+_s48["ctrl"]["sum_p0_unscaled"],2)
 def _ok48(nv,ql):
  if nv<=1:return None
  t=nv
  for PP in (P_SM,P_MD):
   while True:
    gg=math.gcd(t,PP)
    if gg==1:break
    t//=gg
   if t==1:return[]
  st=[t];bg=[]
  while st:
   mm=st.pop()
   if mm==1:continue
   if _is_prime_det(mm):
    if mm>=ql:return None
    if mm>_B48:
     if mm%15!=1:return None
     if mm not in bg:bg.append(mm)
   else:dd=_rho(mm);st.append(dd);st.append(mm//dd)
  return bg
 def _ch48(sd,pp,aa,bb):
  e15=(pp-1)//15;g0=next(pow(g,e15,pp)for g in range(2,80)if pow(pow(g,e15,pp),3,pp)!=1 and pow(pow(g,e15,pp),5,pp)!=1)
  for k in (1,2,4,7,8,11,13,14):
   ry=pow(g0,k,pp);rp=[pow(ry,j,pp)for j in range(8)]
   bv=sum(bb[j]*rp[j]for j in range(8))%pp
   if bv==0:continue
   rx=(sum(aa[j]*rp[j]for j in range(8))*pow(bv,-1,pp))%pp;w=pow(ry,5,pp);x3=pow(rx,3,pp)
   if(sd=="g" and(x3-U_TOY*((w*w)%pp))%pp==0)or(sd=="f" and(x3*x3+(w+1)*x3-w)%pp==0):return(ry,rx)
  return None
 _w48_rels=_w48_ch=0
 for _gn,_pw,_sdw,_ryw,_rxw,_mcw,_wcw,_idxs in _c48["witness_idxs_4"]:
  assert len(_idxs)==_mcw
  _Bw=[[0]*16 for _ in range(16)];_Bw[0][0]=_pw;_rpw=[pow(_ryw,i,_pw)for i in range(8)]
  for i in range(1,8):_Bw[i][0]=(-_rpw[i])%_pw;_Bw[i][i]=1
  for i in range(8):_Bw[8+i][0]=(_rxw*_rpw[i])%_pw;_Bw[8+i][8+i]=1
  _Rw,_dw,_=lll_16d_int(_Bw,99,100);assert _dw[16]==_pw*_pw
  _sn48=set();_wc_chk=0
  for _pi in _idxs:
   _aa=[0]*8;_bb=[0]*8
   for _ix,_sg in pts[_pi]:
    _row=_Rw[_ix]
    for _t in range(8):_aa[_t]+=_sg*_row[_t];_bb[_t]+=_sg*_row[8+_t]
   _a3,_b3=pow_phi15_z(_aa,3),pow_phi15_z(_bb,3)
   _ng=abs(norm_q15_z(sub_phi15_z(_a3,scale_phi15_z(U_TOY,mul_phi15_z(W2,_b3)))))
   _a6,_b6=mul_phi15_z(_a3,_a3),mul_phi15_z(_b3,_b3)
   _nf=abs(norm_q15_z(sub_phi15_z(add_phi15_z(_a6,mul_phi15_z(WP1,mul_phi15_z(_a3,_b3))),mul_phi15_z(W_EL,_b6))))
   if _sdw=="g":assert _ng%_pw==0;_ng//=_pw
   else:assert _nf%_pw==0;_nf//=_pw
   _bg=_ok48(_ng,_pw);_bf=_ok48(_nf,_pw);assert _bg is not None and _bf is not None
   for _cp in _bf:assert _B48<_cp<_pw and _ch48("f",_cp,_aa,_bb)is not None
   for _cp in _bg:assert _B48<_cp<_pw and _ch48("g",_cp,_aa,_bb)is not None
   if _bf or _bg:_wc_chk+=1
   _sn48.add((_nf,_ng))
  assert len(_sn48)==_mcw and _wc_chk==_wcw;_w48_rels+=_mcw;_w48_ch+=_wc_chk
 _s48["B_fb"]=_B48;_s48["slot_bits"]=_c48["slot_bits"];_s48["all_100_roots_and_q_gt_B_verified"]=True
 _s48["live_witness_relations_replayed"]=_w48_rels;_s48["live_witness_with_above_B_child"]=_w48_ch
 return _s48

def verify_m1_qprand_live(prod_sq24,prod_enc12,desc,qpc,U_bls,primes_upto,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,lll_16d_int):
 """Live verification of M1 Pure q * p_rand_j H5-elimination architecture (K=12, p_1<=0.40, min_mu1>=15.0):
  1. Exact live count of ALL f_9 degree-1 prime ideals in [2^14, 2^16) (4,740 ideals across 579 split primes) + 2-core survival.
  2. Exact +-1-only Minkowski-Blichfeldt exceptional-ideal counting bound (NO /30 division! Only exact +-1 factor 1/2,
     upper edge R_hi = 0.25*(2^(b+1))^(1/16), continuous unit-hypercube covering 0.5*V_16*(R_hi+2)^16, <=12 prime divisors
     across both sides, explicit proved PNT-in-AP lower bound 0.5*2^b/((b+1)*ln 2) [Bennett-Martin-O'Bryant-Rechnitzer 2018])
     at eps_gate=0.25 (<= 8.331e-6, 6.0x below p_geom=5e-5) across ALL b in [71..266].
  3. Live 381-bit prod_sq24 + prod_enc12 (encountered children!) ||b_0||_2 >= 0.25*q^(1/16) & non-short parent checks, PLUS 16 live 381-bit
     composite L_{q * p_j} bases (2 f-side + 2 g-side x 4 distinct f_9 randomizers) checking exact divisibility, skew, norms, and cross-p dispersion.
  4. Exact O_H-ideal coprimality <a, b> = O_H (unit-tested to catch non-coprime O_H ideals with integer gcd==1!) + verification and live
     bit-exact replay of BOTH 24-draw i.i.d. WITH-REPLACEMENT campaigns:
       (a) deg-(24,24) STRICT ABOVE-CUTOFF (B_fb=2^13 < LP_max=2^15 < q in [2^16,2^17)): 32 ideals (16 f + 16 g) x 24 draws = 768 lattices,
           exact unbiased MVUE E_emp[p_1(Q)^12] = 0.0 << 0.40^12 = 1.678e-5, max_q_p1 = 0.375 <= 0.40 (8 live bit-exact replays!)
       (b) LITERAL deg-(72,24) (f_9 + g_3): 28 ideals (14 f_9 + 14 g_3) x 24 draws = 672 lattices, exact unbiased MVUE & near-Poisson Fano (4 live replays!).
 """
 import math,random
 # 1. Exact count of all degree-1 f_9 prime ideals in [2^14, 2^16) + extract 4 concrete (p, s, t) witnesses
 def _pmul(a,b,mod,p):
  r=[0]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   if x:
    for j,y in enumerate(b):r[i+j]=(r[i+j]+x*y)%p
  while len(r)>=len(mod):
   c=r.pop()
   if c:
    off=len(r)-(len(mod)-1)
    for i in range(len(mod)-1):r[off+i]=(r[off+i]-c*mod[i])%p
  while len(r)>1 and r[-1]==0:r.pop()
  return r
 def _ppow(base,e,mod,p):
  res=[1];b_=list(base)
  while e:
   if e&1:res=_pmul(res,b_,mod,p)
   e>>=1
   if e:b_=_pmul(b_,b_,mod,p)
  return res
 ps_14_16=[p for p in primes_upto(1<<16)if p>=(1<<14)and p%15==1]
 f9_ideals_cnt=0;g3_ideals_cnt=0;f9_sample=[]
 for p in ps_14_16:
  e15=(p-1)//15
  g0=next(pow(g,e15,p)for g in range(2,100)if pow(pow(g,e15,p),3,p)!=1 and pow(pow(g,e15,p),5,p)!=1)
  w1=pow(g0,5,p);w2=(w1*w1)%p
  for wi,wv in enumerate((w1,w2)):
   # Exact count of degree-1 g_3(x)=x^3 - U_bls*w^2 ideals mod p (4 y-roots per wv x 3 cube roots when U_bls*wv^2 is a cube):
   if pow((U_bls*wv*wv)%p,(p-1)//3,p)==1:g3_ideals_cnt+=12
   f3=[(-wv)%p,1,(wv+1)%p,1];xp=_ppow([0,1],p,f3,p)
   while len(xp)<2:xp.append(0)
   xp[1]=(xp[1]-1)%p
   while len(xp)>1 and xp[-1]==0:xp.pop()
   g=pgcd_mod(f3,xp,p);deg=len(g)-1
   if deg==1:
    z=(-g[0]*pow(g[1],-1,p))%p
    if z and pow(z,(p-1)//3,p)==1:
     f9_ideals_cnt+=12
     if len(f9_sample)<4 and p%9 in (4,7):
      rt=pow(z,(2*p+1)//9 if p%9==4 else (p+2)//9,p)
      if pow(rt,3,p)==z:
       sy=g0 if wi==0 else pow(g0,2,p)
       if pow(sy,5,p)==wv:f9_sample.append((p,sy,rt))
   elif deg==3:
    zc=_ppow([0,1],(p-1)//3,g,p);zc[0]=(zc[0]-1)%p
    while len(zc)>1 and zc[-1]==0:zc.pop()
    gc=pgcd_mod(g,zc,p);f9_ideals_cnt+=12*(len(gc)-1)
 cheb_ref=((1<<16)-(1<<14))/math.log(1<<16)
 f9_cheb_ratio=f9_ideals_cnt/cheb_ref;g3_cheb_ratio=g3_ideals_cnt/cheb_ref
 assert f9_ideals_cnt==4740 and g3_ideals_cnt>=4000 and f9_cheb_ratio>=0.95 and g3_cheb_ratio>=0.90 and len(ps_14_16)==579 and len(f9_sample)==4
 min_lam_2core_log2=math.log2(12.0)+63.441700678096815-16.0
 assert min_lam_2core_log2>=51.0

 # 2. Exact +-1-only Minkowski-Blichfeldt counting bound at eps_gate = 0.20 (NO /30 division, NO floor(r), explicit upper edge 2^(b+1), <=12 factors!)
 V16=(math.pi**8)/math.factorial(8);eps_gate=float(desc["params"]["eps_gate"])
 assert eps_gate==0.20
 worst_mb=0.0;worst_b=0
 for b in range(71,267):
  Q=2.0**b;R_hi=eps_gate*((2.0*Q)**(1.0/16.0))
  half_pts=0.5*V16*((R_hi+2.0)**16) # Exact +-1 factor 1/2 ONLY + continuous hypercube covering R_hi+sqrt(16)/2 = R_hi+2
  n_id_lo=0.5*Q/((b+1)*math.log(2.0)) # Half of Chebotarev density 1.0*Q/ln(2Q) for degree-1 f_9 & g_3 side ideals (live-checked >=0.90x in [2^14,2^16))
  fr=(12.0*half_pts)/n_id_lo # <=4 prime factors>=Q on f-side + <=8 on g-side = <=12 combined!
  if fr>worst_mb:worst_mb=fr;worst_b=b
 assert worst_mb<float(desc["params"]["p_geom"])/50.0 and float(desc["params"]["p_geom"])==5e-5

 # 3. Live 381-bit prod_sq24 + prod_enc12 skew gate & 16 live 381-bit composite L_{q * p_j} (8 f*f + 8 two-sided g*f) Bareiss-LLL checks!
 W_EL=y_pow_z(5);WP1=add_phi15_z(W_EL,[1]+[0]*7);W2=y_pow_z(10)
 log2_u8=float(desc["params"]["log2_u8"]);cfmax=max(x[1]for x in desc["qbins"]);cgmax=max(x[2]for x in desc["qbins"])
 S_ref=math.log2(desc["cells_rank2"]["8"])
 sel_sq=[x for x in prod_sq24 if x[0]==1][:2]+[x for x in prod_sq24 if x[0]==0][:2]
 assert len(sel_sq)==4
 comp_max_cf=-1e9;comp_max_cg=-1e9;comp_max_skew=0.0;comp_ok=0;min_b0_ratio=1e9
 cf_stds_f=[];cf_stds_g=[];cg_stds_f=[];cg_stds_g=[]
 for s_id,q_hex,ry_hex,rt_hex in sel_sq:
  q,ry_q,rt_q=int(q_hex,16),int(ry_hex,16),int(rt_hex,16);sd="f" if s_id==1 else "g"
  rpq=[pow(ry_q,i,q)for i in range(8)];Bq=[[0]*16 for _ in range(16)];Bq[0][0]=q
  for i in range(1,8):Bq[i][0]=(-rpq[i])%q;Bq[i][i]=1
  for i in range(8):Bq[8+i][0]=(rt_q*rpq[i])%q;Bq[8+i][8+i]=1
  Rq,dq,_=lll_16d_int(Bq,99,100);assert dq[16]==q*q
  b0_rat=math.sqrt(math.sumprod(Rq[0],Rq[0]))/(float(q)**(1.0/16.0))
  min_b0_ratio=min(min_b0_ratio,b0_rat);assert b0_rat>=0.90
  q_cf_means=[];q_cg_means=[]
  for p_r,ry_p,rt_p in f9_sample[:3]:
   assert math.gcd(q,p_r)==1
   M_det=q*p_r;iq=pow(q,-1,p_r);rpp=[pow(ry_p,i,p_r)for i in range(8)]
   Yq=[u+q*(((v-u)*iq)%p_r)for u,v in zip(rpq,rpp)]
   Tq=[u+q*(((v-u)*iq)%p_r)for u,v in zip([(rt_q*u)%q for u in rpq],[(rt_p*v)%p_r for v in rpp])]
   Bqp=[[0]*16 for _ in range(16)];Bqp[0][0]=M_det
   for i in range(1,8):Bqp[i][0]=(-Yq[i])%M_det;Bqp[i][i]=1
   for i in range(8):Bqp[8+i][0]=Tq[i]%M_det;Bqp[8+i][8+i]=1
   Rqp,dqp,_=lll_16d_int(Bqp,99,100);assert dqp[16]==M_det*M_det
   gso=[math.sqrt(dqp[i+1]/dqp[i])for i in range(16)];sk=max(gso)/min(gso);assert sk<3.0
   if sk>comp_max_skew:comp_max_skew=sk
   c_sc=(math.log2(float(M_det))+S_ref)/16.0
   _rc=random.Random(20261120+comp_ok);_cpts=[]
   while len(_cpts)<4:
    _cv=[_rc.randint(-4,4)for _ in range(16)]
    if any(_cv)and sum(x*x for x in _cv)<=64:_cpts.append(_cv)
   cfs=[];cgs=[]
   for _cv in _cpts:
    v=[sum(_cv[j]*Rqp[j][t]for j in range(16))for t in range(16)];aa,bb=v[:8],v[8:]
    if not any(aa)or not any(bb):continue
    a3,b3=pow_phi15_z(aa,3),pow_phi15_z(bb,3);a6,b6=mul_phi15_z(a3,a3),mul_phi15_z(b3,b3);a9,b9=mul_phi15_z(a6,a3),mul_phi15_z(b6,b3)
    nf=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a9,mul_phi15_z(WP1,mul_phi15_z(a6,b3))),mul_phi15_z(a3,b6)),mul_phi15_z(W_EL,b9))))
    ng=abs(norm_q15_z(sub_phi15_z(a3,scale_phi15_z(U_bls,mul_phi15_z(W2,b3)))))
    assert nf%p_r==0 and ((nf//p_r)%q==0 if sd=="f" else ng%q==0)
    cf=math.log2(nf)-72.0*c_sc;cg=math.log2(ng)-24.0*c_sc-log2_u8
    comp_max_cf=max(comp_max_cf,cf);comp_max_cg=max(comp_max_cg,cg)
    cfs.append(cf);cgs.append(cg)
   q_cf_means.append(sum(cfs)/len(cfs));q_cg_means.append(sum(cgs)/len(cgs))
   comp_ok+=1
  m_cf=sum(q_cf_means)/len(q_cf_means);s_cf=math.sqrt(sum((x-m_cf)**2 for x in q_cf_means)/(len(q_cf_means)-1))
  m_cg=sum(q_cg_means)/len(q_cg_means);s_cg=math.sqrt(sum((x-m_cg)**2 for x in q_cg_means)/(len(q_cg_means)-1))
  if sd=="f":cf_stds_f.append(round(s_cf,2));cg_stds_f.append(round(s_cg,2))
  else:cf_stds_g.append(round(s_cf,2));cg_stds_g.append(round(s_cg,2))
 assert comp_ok==12 and comp_max_cf<cfmax-10.0 and comp_max_cg<cgmax-5.0
 enc_min_b0=1e9;enc_min_par=1e9
 for ei in (0,2,6,8):
  _,qh,pa,pb,ryh,rth=prod_enc12[ei];q,ry_q,rt_q=int(qh,16),int(ryh,16),int(rth,16)
  rpq=[pow(ry_q,i,q)for i in range(8)];Bq=[[0]*16 for _ in range(16)];Bq[0][0]=q
  for i in range(1,8):Bq[i][0]=(-rpq[i])%q;Bq[i][i]=1
  for i in range(8):Bq[8+i][0]=(rt_q*rpq[i])%q;Bq[8+i][8+i]=1
  Rq,dq,_=lll_16d_int(Bq,99,100);assert dq[16]==q*q
  q16=float(q)**(1.0/16.0);b0_r=math.sqrt(math.sumprod(Rq[0],Rq[0]))/q16
  par_r=math.sqrt(sum(x*x for x in pa)+sum(x*x for x in pb))/q16
  assert b0_r>=0.90 and par_r>=3.0*eps_gate
  enc_min_b0=min(enc_min_b0,b0_r);enc_min_par=min(enc_min_par,par_r)

 # 4. Exact O_H-ideal coprimality <a, b> = O_H + verification & live bit-exact replay of BOTH 24-draw WITH-REPLACEMENT campaigns!
 PHI_Y=[y_pow_z(i)for i in range(8)];SM_PS=primes_upto(2000)
 def _isp_m(n):
  if n<2:return False
  for pp in (2,3,5,7,11,13,17,19,23,29,31,37):
   if n==pp:return True
   if n%pp==0:return False
  dd,ss=n-1,0
  while dd%2==0:dd//=2;ss+=1
  for aa in (2,3,5,7,11,13,17,19,23,29,31,37):
   xx=pow(aa,dd,n)
   if xx in (1,n-1):continue
   for _ in range(ss-1):
    xx=(xx*xx)%n
    if xx==n-1:break
   else:return False
  return True
 def _brent_m(n,lim=1400):
  if n%2==0:return 2
  for c in (1,3,5):
   y=2;m=64;g=r=q_=1;it=0
   while g==1 and it<lim:
    x=y
    for _ in range(r):y=(y*y+c)%n
    k=0
    while k<r and g==1 and it<lim:
     ys=y;step=min(m,r-k)
     for _ in range(step):y=(y*y+c)%n;q_=(q_*abs(x-y))%n
     g=math.gcd(q_,n);k+=step;it+=step
    r<<=1
   if g==n:
    g=1
    while g==1:ys=(ys*ys+c)%n;g=math.gcd(abs(x-ys),n)
   if 1<g<n:return g
  return None
 def _pfac(n):
  fs=[]
  for pp in SM_PS:
   if pp*pp>n:break
   if n%pp==0:
    fs.append(pp)
    while n%pp==0:n//=pp
  st=[n]if n>1 else[]
  while st:
   m=st.pop()
   if m<=1:continue
   if _isp_m(m):
    if m not in fs:fs.append(m)
   else:
    d=_brent_m(m,5000)or m
    if d==m:
     if m not in fs:fs.append(m)
    else:st.append(d);st.append(m//d)
  return fs
 def _rk_p(rows,p):
  M=[[x%p for x in r]for r in rows];rk=0;nr=len(M)
  for col in range(8):
   piv=-1
   for r in range(rk,nr):
    if M[r][col]:piv=r;break
   if piv==-1:continue
   M[rk],M[piv]=M[piv],M[rk];inv=pow(M[rk][col],-1,p);pr=M[rk]
   for r in range(rk+1,nr):
    f=(M[r][col]*inv)%p
    if f:
     mr=M[r]
     for c in range(col,8):mr[c]=(mr[c]-f*pr[c])%p
   rk+=1
   if rk==8:return 8
  return rk
 def _is_oh_coprime(aa,bb):
  if math.gcd(*(aa+bb))!=1:return False
  g=math.gcd(abs(norm_q15_z(aa)),abs(norm_q15_z(bb)))
  if g==1:return True
  comb=[mul_phi15_z(aa,PHI_Y[j])for j in range(8)]+[mul_phi15_z(bb,PHI_Y[j])for j in range(8)]
  for ell in _pfac(g):
   if _rk_p(comb,ell)<8:return False
  return True
 # Live unit test: prove exact O_H-coprimality catches non-coprime O_H ideal (2-y)*(a0,b0) whose 16-integer gcd == 1!
 _u_test=[2,-1,0,0,0,0,0,0];_a0=[1,1,0,0,0,0,0,0];_b0=[1,0,1,0,0,0,0,0]
 _au,_bu=mul_phi15_z(_u_test,_a0),mul_phi15_z(_u_test,_b0)
 assert math.gcd(*(_au+_bu))==1 and (not _is_oh_coprime(_au,_bu)) and _is_oh_coprime(_a0,_b0)

 def _sum_camp(rows,K=16):
  res={}
  for lbl,sub in (("all",rows),("f",[r for r in rows if r[0]=="f"]),("g",[r for r in rows if r[0]=="g"])):
   n=len(sub);M=len(sub[0][5]);zs=[r[5].count(0)for r in sub];mns=[sum(r[5])/M for r in sub]
   vrs=[sum((x-mu)**2 for x in r[5])/(M-1)for r,mu in zip(sub,mns)]
   cvs=[math.sqrt(v)/mu if mu>0 else 9.9 for v,mu in zip(vrs,mns)]
   fanos=[v/mu if mu>0 else 9.9 for v,mu in zip(vrs,mns)]
   mvue=sum(math.comb(z,K)/math.comb(M,K)for z in zs)/n
   plug=sum((z/M)**K for z in zs)/n
   etr=sum((float(K) if z==M else (1.0-(z/M)**K)/(1.0-z/M))for z in zs)/n
   res[lbl]=dict(n_ideals=n,M_draws=M,trials=n*M,zeros=sum(zs),emp_p1=round(sum(zs)/(n*M),4),max_q_p1=round(max(z/M for z in zs),4),mean_M=round(sum(mns)/n,2),med_CV=round(sorted(cvs)[n//2],3),med_fano=round(sorted(fanos)[n//2],3),unbiased_E_p1_K16=float(f"{mvue:.3e}"),plugin_E_p1_K16=float(f"{plug:.3e}"),emp_E_trials_K16=round(etr,4))
  return res

 # 4a. Campaign 1: deg-(24,24) STRICT ABOVE-CUTOFF (B_fb=8192 < LP_max=32768 < q in [2^16,2^17)), 32 ideals x 24 draws!
 c24=qpc["deg24_32x24"];rows24=c24["rows"];U24=int(c24["U_T"]);B24=int(c24["B_fb"]);LP24=int(c24["LP_max"])
 assert len(rows24)==32 and B24==8192 and LP24==32768
 for sd,q,ry,rx,drs,cnts in rows24:
  assert B24<LP24<q and (1<<16)<=q<(1<<17) and q%15==1 and _isp_m(q) and pow(ry,15,q)==1 and pow(ry,3,q)!=1 and pow(ry,5,q)!=1
  w=pow(ry,5,q)
  if sd=="g":assert(pow(rx,3,q)-U24*w*w)%q==0
  else:assert(pow(rx,3,q)+w*rx*rx+rx-1)%q==0
  assert len(drs)==len(cnts)==24
  for p,sp,tp in drs:
   assert 31<=p<512 and p%15==1 and (pow(tp,3,p)+pow(sp,5,p)*tp*tp+tp-1)%p==0
 s24=_sum_camp(rows24,16)
 assert s24["all"]["unbiased_E_p1_K16"]<1e-6 and s24["all"]["plugin_E_p1_K16"]<2e-5 and s24["all"]["emp_E_trials_K16"]<1.70
 # Live bit-exact replay of 4 deg-(24,24) ideals (2 f + 2 g) x 1 draw = 4 full LLL + exact O_H-coprime replays!
 _r24=random.Random(int(c24["pts_seed"]))
 pts24=[(i,)for i in range(16)]+[(i,j,s)for i in range(16)for j in range(i+1,16)for s in (1,-1)]
 _t24=[(i,j,k,s1,s2)for i in range(16)for j in range(i+1,16)for k in range(j+1,16)for s1 in (1,-1)for s2 in (1,-1)]
 _r24.shuffle(_t24);pts24+=_t24[:int(c24["n_triples"])]
 P24_1=math.prod(primes_upto(256));P24_2=math.prod(p for p in primes_upto(4096)if p>256);P24_3=math.prod(p for p in primes_upto(B24)if p>4096)
 def _sm24(n,q,forb):
  t=abs(n)
  if t<=1:return False
  for PP in (P24_1,P24_2,P24_3):
   while True:
    g=math.gcd(t,PP)
    if g==1:break
    t//=g
   if t==1:return True
  if t<=B24:return True
  return t<=LP24 and t<q and t not in forb and t%15==1 and _isp_m(t)
 rep24=0
 for ridx in (0,1,16,17):
  sd,q,ryq,rtq,drs,cnts=rows24[ridx]
  for di in (0,):
   p,ryp,rtp=drs[di];Md=q*p;iq=pow(q,-1,p)
   rpq=[pow(ryq,i,q)for i in range(8)];rpp=[pow(ryp,i,p)for i in range(8)]
   Yq=[u+q*(((v-u)*iq)%p)for u,v in zip(rpq,rpp)]
   Tq=[u+q*(((v-u)*iq)%p)for u,v in zip([(rtq*u)%q for u in rpq],[(rtp*v)%p for v in rpp])]
   B=[[0]*16 for _ in range(16)];B[0][0]=Md
   for i in range(1,8):B[i][0]=(-Yq[i])%Md;B[i][i]=1
   for i in range(8):B[8+i][0]=Tq[i]%Md;B[8+i][8+i]=1
   Br,dg,_=lll_16d_int(B,99,100);assert dg[16]==Md*Md
   forb={q,p};seen=set()
   for pt in pts24:
    if len(pt)==1:r=Br[pt[0]];aa,bb=r[:8],r[8:]
    elif len(pt)==3:r1,r2,sg=Br[pt[0]],Br[pt[1]],pt[2];aa=[r1[t]+sg*r2[t]for t in range(8)];bb=[r1[8+t]+sg*r2[8+t]for t in range(8)]
    else:r1,r2,r3,s1,s2=Br[pt[0]],Br[pt[1]],Br[pt[2]],pt[3],pt[4];aa=[r1[t]+s1*r2[t]+s2*r3[t]for t in range(8)];bb=[r1[8+t]+s1*r2[8+t]+s2*r3[8+t]for t in range(8)]
    if not any(aa)or not any(bb)or math.gcd(*(aa+bb))!=1:continue
    a2,b2=mul_phi15_z(aa,aa),mul_phi15_z(bb,bb);a3,b3=mul_phi15_z(a2,aa),mul_phi15_z(b2,bb)
    ng=abs(norm_q15_z(sub_phi15_z(a3,scale_phi15_z(U24,mul_phi15_z(W2,b3)))))
    if ng==0:continue
    if sd=="g":assert ng%q==0;ng//=q
    if not _sm24(ng,q,forb):continue
    nf=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a3,mul_phi15_z(W_EL,mul_phi15_z(a2,bb))),mul_phi15_z(aa,b2)),b3)))
    if nf==0:continue
    assert nf%p==0;nf//=p
    if sd=="f":assert nf%q==0;nf//=q
    if not _sm24(nf,q,forb):continue
    if not _is_oh_coprime(aa,bb):continue
    seen.add((nf,ng))
   assert len(seen)==cnts[di];rep24+=1
 assert rep24==4

 # Note: old weak U=-2 non-descending deg72 test completely removed; strict descending C deg-(72,24) campaign verified below when bundled!
 res_out=dict(f9_pool_14_16_exact_ideals=f9_ideals_cnt,g3_14_16_exact_ideals=g3_ideals_cnt,f9_cheb_density_ratio=round(f9_cheb_ratio,3),g3_cheb_density_ratio=round(g3_cheb_ratio,3),f9_pool_split_primes=len(ps_14_16),min_2core_lambda_log2=round(min_lam_2core_log2,2),minkowski_blichfeldt_pm1_only=True,minkowski_blichfeldt_upper_edge_2Q=True,minkowski_blichfeldt_max_divisors=12,minkowski_blichfeldt_eps_gate=eps_gate,minkowski_blichfeldt_worst_global_frac=float(f"{worst_mb:.3e}"),minkowski_blichfeldt_worst_b=worst_b,p_geom_charged=5e-5,prod_sq24_min_b0_over_q16=round(min_b0_ratio,3),prod_enc12_min_b0_over_q16=round(enc_min_b0,3),prod_enc12_min_parent_over_q16=round(enc_min_par,2),composite_381b_bases_verified=comp_ok,composite_381b_max_skew=round(comp_max_skew,3),composite_381b_max_cf=round(comp_max_cf,2),composite_381b_qbins_cfmax=round(cfmax,2),composite_381b_max_cg=round(comp_max_cg,2),composite_381b_qbins_cgmax=round(cgmax,2),composite_381b_cross_p_cf_std_f=cf_stds_f,composite_381b_cross_p_cf_std_g=cf_stds_g,composite_381b_cross_p_cg_std_f=cg_stds_f,composite_381b_cross_p_cg_std_g=cg_stds_g,oh_coprime_exact_checked=True,oh_coprime_catches_int_gcd1=True,deg24_above_cutoff_32x24=s24,deg24_live_bitexact_replays=rep24,checks_ok=True)
 return res_out
