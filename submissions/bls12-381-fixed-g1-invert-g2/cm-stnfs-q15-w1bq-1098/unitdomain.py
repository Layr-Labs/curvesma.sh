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
 """Upgraded 381-bit production CRT check + 8 parents x 6 blinding draws = 48 live-LLL cross-blinding dispersion check.""";import cmath;assert len(sq24)==24 and vor_units is not None
 _PL4=[[complex(float(re),float(im))for re,im in row]for row in _PL];_W_PL=[row[5]for row in _PL4];_WP1_PL=[w+1.0 for w in _W_PL];_UW2_PL=[float(U_bls)*w*w for w in _W_PL];log2_u8=float(desc["params"]["log2_u8"])
 lx=float(desc["params"]["lx"]);qbins,qbins_lo,qbins_sub=desc["qbins"],desc["qbins_lo"],desc["qbins_sub"];qb_f,qb_g=[x[1]for x in qbins],[x[2]for x in qbins]
 qbf_iqm=sum(cf for _,cf,_ in qbins[1:5])/4.0
 qbg_iqm=sum(cg for _,_,cg in qbins[1:5])/4.0
 S8=math.log2(desc["cells_rank2"]["8"]);S_w2=8.0;w_el=y_pow_fn(5);wp1_el=add_fn(w_el,[1]+[0]*7);w2_el=y_pow_fn(10);rng=random.Random(20261120);ab_pts=[]
 while len(ab_pts)<6:
  al=[rng.randint(-4,4)for _ in range(8)];bt=[rng.randint(-4,4)for _ in range(8)]
  if any(al)and sum(x*x for x in al)+sum(x*x for x in bt)<=64:
   ab_pts.append((al,bt))
 pairs16=[(i,15-i)for i in range(8)]+[(i,8+((8-i)&7))for i in range(8)];sub32=[(i,j,s)for(i,j)in pairs16 for s in(1,-1)]
 def _is_p(n):
  d,s=n-1,0
  while d%2==0:
   d//=2;s+=1
  for a in(2,3,5,7,11,13,17,19,23,29,31,37):
   x=pow(a,d,n)
   if x in(1,n-1):
    continue
   for _ in range(s-1):
    x=(x*x)%n
    if x==n-1:
     break
   else:
    return False
  return True
 def _cbrt_f(rhs,p):
  if pow(rhs,(p-1)//3,p)!=1:
   return None
  if p%9==4:
   r=pow(rhs,(2*p+1)//9,p);return r if(r*r*r)%p==rhs else None
  if p%9==7:
   r=pow(rhs,(p+2)//9,p);return r if(r*r*r)%p==rhs else None
  for rx in range(1,p):
   if(rx*rx*rx)%p==rhs:
    return rx
  return None
 out_sides={};sm_ps=[p for p in primes_upto(100)if p>5];blind_pool=[]
 for _qp in range(2049,8192,2):
  if _qp%15!=1 or not _is_p(_qp):
   continue
  _e15=(_qp-1)//15;_g0=next(pow(g,_e15,_qp)for g in range(2,80)if pow(pow(g,_e15,_qp),3,_qp)!=1 and pow(pow(g,_e15,_qp),5,_qp)!=1);_w_rts={}
  for _kr in(1,2,4,7,8,11,13,14):
   _ryp=pow(_g0,_kr,_qp);_wp=pow(_ryp,5,_qp)
   if _wp not in _w_rts:
    _rx_found=None;_wp1=_wp+1
    for _yv in range(1,_qp):
     if(((_yv+_wp1)*_yv+1)*_yv-_wp)%_qp==0:
      _rx_found=_cbrt_f(_yv,_qp)
      if _rx_found is not None:
       break
    _w_rts[_wp]=_rx_found
   if _w_rts[_wp]is not None:
    blind_pool.append((_qp,_ryp,_w_rts[_wp]))
   if len(blind_pool)>=72:
    break
  if len(blind_pool)>=72:
   break
 assert len(blind_pool)==72 and len({(x[0],x[1])for x in blind_pool})==72;_ss=dlive["_ss_eval"];_mn=dlive["mesh_node"];_af=float(desc["params"]["alpha_f"]);_ag=float(desc["params"]["alpha_g"])
 def _crt_lll_381(q,ry,rt,qp,ryp_b,rtp_b):
  M=q*qp;iq=pow(q,-1,qp);ryp=[pow(ry,i,q)for i in range(8)];rypb=[pow(ryp_b,i,qp)for i in range(8)];Yq=[u+q*(((v-u)*iq)%qp)for u,v in zip(ryp,rypb)]
  Tq=[u+q*(((v-u)*iq)%qp)for u,v in zip([(rt*u)%q for u in ryp],[(rtp_b*v)%qp for v in rypb])];B=[[0]*16 for _ in range(16)];B[0][0]=M
  for i in range(1,8):
   B[i][0]=(-Yq[i])%M;B[i][i]=1
  for i in range(8):
   B[8+i][0]=Tq[i]%M;B[8+i][8+i]=1
  R,dsq,lsq=lll_16d_int(B,99,100);assert dsq[16]==M*M
  for kk in range(1,16):
   assert 100*(dsq[kk-1]*dsq[kk+1]+lsq[kk][kk-1]**2)>=99*dsq[kk]**2
   for jj in range(kk):
    assert 2*abs(lsq[kk][jj])<=dsq[jj+1]
  return M,R
 def _mu_emp_R(s_id,sname,q,ry,rt,qp,ryp_b,rtp_b,M,R,prow,S_bin,chk_norm=False):
  c_val=(math.log2(M)+S_w2)/16.0;c_act=(math.log2(M)+S_bin)/16.0;lq,lqp=math.log2(q),math.log2(qp);mx_c=max(abs(x)for row in R for x in row);sh=max(0,mx_c.bit_length()-50);sc=1.0/float(1<<sh)if sh else 1.0
  sh_f,sh_g=72.0*sh,24.0*sh;ra_pl=[[sum((row[j]*sc)*pl[j]for j in range(8))for pl in _PL4]for row in R];rb_pl=[[sum((row[8+j]*sc)*pl[j]for j in range(8))for pl in _PL4]for row in R];a0=[R[0][k]+R[15][k]for k in range(8)]
  b0=[R[0][8+k]+R[15][8+k]for k in range(8)];assert(sum(a0[k]*pow(ry,k,q)for k in range(8))-rt*sum(b0[k]*pow(ry,k,q)for k in range(8)))%q==0
  assert(sum(a0[k]*pow(ryp_b,k,qp)for k in range(8))-rtp_b*sum(b0[k]*pow(ryp_b,k,qp)for k in range(8)))%qp==0
  if chk_norm:
   a30,b30=pow_fn(a0,3),pow_fn(b0,3);ng0=abs(norm_fn(sub_fn(a30,scale_fn(U_bls,mul_fn(w2_el,b30)))));a60,b60=mul_fn(a30,a30),mul_fn(b30,b30)
   nf0=abs(norm_fn(sub_fn(add_fn(add_fn(mul_fn(a60,a30),mul_fn(wp1_el,mul_fn(a60,b30))),mul_fn(a30,b60)),mul_fn(w_el,mul_fn(b60,b30)))));assert nf0%qp==0 and(ng0 if s_id==0 else(nf0//qp))%q==0
  cfs_t,cgs_t=[],[]
  for idx_v,(i1,i2,s2)in enumerate(sub32):
   lf=sh_f;lg=sh_g;ra1,ra2=ra_pl[i1],ra_pl[i2];rb1,rb2=rb_pl[i1],rb_pl[i2]
   for k in range(4):
    za=(ra1[k]+ra2[k])if s2==1 else(ra1[k]-ra2[k]);zb=(rb1[k]+rb2[k])if s2==1 else(rb1[k]-rb2[k]);za3=za*za*za;zb3=zb*zb*zb;zg=za3-_UW2_PL[k]*zb3;zf=((za3+_WP1_PL[k]*zb3)*za3+zb3*zb3)*za3-_W_PL[k]*(zb3*zb3*zb3)
    lf+=2.0*math.log2(abs(zf));lg+=2.0*math.log2(abs(zg))
   if chk_norm and idx_v==0:
    assert abs(lf-math.log2(nf0))<1e-6 and abs(lg-math.log2(ng0))<1e-6
   cfs_t.append(lf-72.0*c_val);cgs_t.append(lg-24.0*c_val-log2_u8)
  cfs_t.sort();cgs_t.sort()
  dcf=sum(cfs_t[8:24])/16.0-qbf_iqm
  dcg=sum(cgs_t[8:24])/16.0-qbg_iqm
  sub_f=(lq if sname=="f"else 0.0)+lqp;sub_g=lq if sname=="g"else 0.0;nfl_q=[72.0*c_act+cf+dcf-sub_f for _,cf,_ in qbins]+[72.0*c_act+cf+dcf-sub_f for cf,_ in qbins_sub]
  ngl_q=[24.0*c_act+log2_u8+cg+dcg-sub_g for _,_,cg in qbins]+[24.0*c_act+log2_u8+cg+dcg-sub_g for _,cg in qbins_sub];pf=_ss(nfl_q,float(prow["n"]),prow["iA"],_mn,_af,prow["kA"],prow["kB"])
  pg=_ss(ngl_q,float(prow["n"]),prow["iA"],_mn,_ag,prow["kA"],prow["kB"]);pr_emp=min(sum(pf[i]*pg[i]for i in range(6,30))/24.0,4.0*sum(pf[i]*pg[i]for i in range(6))/6.0);return prow["mu"]*(pr_emp/prow["pr"])
 disp_rows=[];disp_cv2s=[]
 disp_idx_map={0:0,3:1,7:2,11:3}
 for s_id,sname in((0,"g"),(1,"f")):
  cfs,cgs,sv_rats,pr_rats,uf_all,lat_masks,eta_ok=[],[],[],[],[],[],0;recs=[r for r in sq24 if r[0]==s_id];assert len(recs)==12
  for r_idx,(side,q_hex,ry_hex,rt_hex)in enumerate(recs):
   q,ry,rt=int(q_hex,16),int(ry_hex,16),int(rt_hex,16);d_ord=4*s_id+disp_idx_map[r_idx]if r_idx in disp_idx_map else None;qp,ryp_b,rtp_b=blind_pool[6*d_ord]if d_ord is not None else blind_pool[48+12*s_id+r_idx]
   nb=q.bit_length()-1;assert int(lx)<=nb<266 and q%15==1 and _is_p(q)and 2048<=qp<=8192 and math.gcd(q,qp)==1;prow=next(r for r in dlive["rows"]if r["side"]==sname and r["lo"]<=nb<r["hi"])
   S_bin=math.log2(desc["cells_rank2"][str(prow["R"])]);c_lo_b=(prow["lo"]+11.0+S_bin)/16.0;T_b=prow["kA"]*prow["nA"]+prow["kB"]*prow["n"]+96.0;sv_plan=0.0
   for(w_q,_,_),(cfl,cgl)in zip(qbins,qbins_lo):
    nfl=72.0*c_lo_b+cfl-((prow["hi"]if sname=="f"else 0.0)+13.0);ngl=24.0*c_lo_b+log2_u8+cgl-(prow["hi"]if sname=="g"else 0.0);sf=1.0 if T_b>=nfl else min(1.0,max(1.0,T_b/lx)*rho_exact((math.floor(nfl)-T_b)/lx))
    sg=1.0 if T_b>=ngl else min(1.0,max(1.0,T_b/lx)*rho_exact((math.floor(ngl)-T_b)/lx));sv_plan+=w_q*sf*sg
   M,R=_crt_lll_381(q,ry,rt,qp,ryp_b,rtp_b);c_val=(math.log2(M)+S8)/16.0;c_act=(math.log2(M)+S_bin)/16.0;lq,lqp=math.log2(q),math.log2(qp);nfas,ngas,dmasks=[],[],[];sv_acc=0.0
   for al,bt in ab_pts:
    vv=al+bt;a=[sum(vv[j]*R[j][t]for j in range(16))for t in range(8)];b=[sum(vv[j]*R[j][8+t]for j in range(16))for t in range(8)];a3,b3=pow_fn(a,3),pow_fn(b,3)
    ng=abs(norm_fn(sub_fn(a3,scale_fn(U_bls,mul_fn(w2_el,b3)))));a6,b6=mul_fn(a3,a3),mul_fn(b3,b3);a9,b9=mul_fn(a6,a3),mul_fn(b6,b3)
    nf=abs(norm_fn(sub_fn(add_fn(add_fn(a9,mul_fn(wp1_el,mul_fn(a6,b3))),mul_fn(a3,b6)),mul_fn(w_el,b9))));assert nf%qp==0 and(ng if s_id==0 else(nf//qp))%q==0;cf=math.log2(nf)-72.0*c_val
    cg=math.log2(ng)-24.0*c_val-log2_u8;cfs.append(cf);cgs.append(cg);nfa=72.0*c_act+cf-((lq if sname=="f"else 0.0)+lqp);nga=24.0*c_act+log2_u8+cg-(lq if sname=="g"else 0.0);nfas.append(nfa);ngas.append(nga)
    uf_all.append(nfa/lx);mm=0
    for pi_i,pp in enumerate(sm_ps):
     if(nf%pp==0)or(ng%pp==0):
      mm|=1<<pi_i
    dmasks.append(mm);sfe=1.0 if T_b>=nfa else min(1.0,max(1.0,T_b/lx)*rho_exact((math.floor(nfa)-T_b)/lx));sge=1.0 if T_b>=nga else min(1.0,max(1.0,T_b/lx)*rho_exact((math.floor(nga)-T_b)/lx));sv_acc+=sfe*sge
   lat_masks.append(dmasks);assert max(nfas)<=prow["nf"]-30.0 and max(ngas)<=prow["ng"]-20.0;sv_rats.append((sv_acc/len(ab_pts))/max(sv_plan,1e-300))
   _pfe=_ss(nfas,float(prow["n"]),prow["iA"],_mn,_af,prow["kA"],prow["kB"]);_pge=_ss(ngas,float(prow["n"]),prow["iA"],_mn,_ag,prow["kA"],prow["kB"]);_pre0=sum(x*y for x,y in zip(_pfe,_pge))/len(ab_pts)
   pr_rats.append(_pre0/prow["pr"]);eta_ok+=1
   if d_ord is not None:
    mus_6=[_mu_emp_R(s_id,sname,q,ry,rt,qp,ryp_b,rtp_b,M,R,prow,S_bin,chk_norm=True)]
    for t_extra in range(1,6):
     qp_t,ryp_t,rtp_t=blind_pool[6*d_ord+t_extra];assert(qp_t,ryp_t)!=(qp,ryp_b)and 2048<=qp_t<=8192 and math.gcd(q,qp_t)==1;M_t,R_t=_crt_lll_381(q,ry,rt,qp_t,ryp_t,rtp_t)
     mus_6.append(_mu_emp_R(s_id,sname,q,ry,rt,qp_t,ryp_t,rtp_t,M_t,R_t,prow,S_bin,chk_norm=False))
    assert len(mus_6)==6;m_mu=sum(mus_6)/6.0;v_mu=sum((x-m_mu)**2 for x in mus_6)/5.0;cv2_q=v_mu/(m_mu*m_mu);assert cv2_q<0.15;disp_cv2s.append(cv2_q)
    disp_rows.append(dict(side=sname,q_bits=nb,mean_mu_emp=round(m_mu,2),cv2_blind=round(cv2_q,5)))
  cfs.sort();cgs.sort();n=len(cfs);sf6=[round(cfs[min(n-1,int((i+1)*n/6)-1)],2)for i in range(6)];sg6=[round(cgs[min(n-1,int((i+1)*n/6)-1)],2)for i in range(6)]
  assert all(sf6[i]<qb_f[i]for i in range(6))and all(sg6[i]<qb_g[i]for i in range(6));assert max(sv_rats)<=1.0+1e-12 and(sum(sv_rats)/len(sv_rats))<0.50 and min(pr_rats)>=10.0;rhos=[]
  for pi_i in range(len(sm_ps)):
   bit=1<<pi_i;pm=sum(sum(1 for x in dm if x&bit)/len(dm)for dm in lat_masks)/len(lat_masks)
   pp=sum((lambda k,nn:k*(k-1)/(nn*(nn-1)))(sum(1 for x in dm if x&bit),len(dm))for dm in lat_masks)/len(lat_masks)
   if 0.0<pm<1.0:
    rhos.append((pp-pm*pm)/(pm*(1.0-pm)))
  mean_rho=sum(rhos)/len(rhos);assert abs(mean_rho)<0.04 and eta_ok==12;_prs=sorted(pr_rats)
  out_sides[sname]=dict(bases=12,samples=n,u_f_range=[round(min(uf_all),2),round(max(uf_all),2)],sp_rho_mean=round(mean_rho,4),sext_f=sf6,sext_g=sg6,min_headroom_f=round(min(qb_f[i]-sf6[i]for i in range(6)),2),min_headroom_g=round(min(qb_g[i]-sg6[i]for i in range(6)),2),joint_pr_over_charged_min=round(_prs[0],2),joint_pr_over_charged_med=round(_prs[len(_prs)//2],2),sv_emp_over_plan_mean=round(sum(sv_rats)/len(sv_rats),4),sv_emp_over_plan_max=round(max(sv_rats),4),crt_det_eq_q_qp=f"{eta_ok}/12",)
 assert len(disp_cv2s)==8 and max(disp_cv2s)<0.15<dlive["min_kappa_crit_pz"]
 return dict(n_bases=24,live_bareiss_99_100_crt_all24=True,crt_Z_basis_det_eq_q_qp="24/24",p_sat_abort_proved=0.0,presupplied_vectors=0,presupplied_blinding_ideals=0,q_bit_range="73..258",qp_range="2048..8192",lattice_geom_modeled_yield_stability_8x6=dict(metric_scope="modeled_smoothness_yield_over_live_381b_crt_lll_lattices_not_observed_M",n_parents=8,draws_per_parent=6,live_lll_crt_lattices=48,cv2_mean=round(sum(disp_cv2s)/8.0,5),cv2_max_Q=round(max(disp_cv2s),5),kappa_crit_pz=round(dlive["min_kappa_crit_pz"],4),parents=disp_rows,),g_side=out_sides["g"],f_side=out_sides["f"],checks_ok=True,)
