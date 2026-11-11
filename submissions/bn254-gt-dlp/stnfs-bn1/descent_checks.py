"""Descent checks for the evidence / audit scripts (pure Python). Data: descent_stats.json next to this file
(keys prime_structure, omitted_law, pools, nodes, numbers) or the same keys inside common12.stats()."""
import json,os,math,random
from math import log2,ceil,floor
import descent16 as DS,retain16 as RT,counts16
_HERE=os.path.dirname(os.path.abspath(__file__));_ST=None
def stats():
    global _ST
    if _ST is None:
        p=os.path.join(_HERE,'descent_stats.json')
        if os.path.exists(p):_ST=json.load(open(p))
        else:import common12;_ST=common12.stats()
        RT.load_table(_ST)
    return _ST
N_MAX_LOG2=20.
def tree():
    """[7] full evaluation of the shipped plan: two-log work, conditions, widths, node count, success terms."""
    st=stats();Pm=DS.P();res=DS.evaluate(Pm,DS.PLAN16,DS.BOOT16);bt=res['boot'];det=res['DET'].values();sv=[d for d in det if d['how']=='sieve']
    EN,_=counts16.expected_nodes(Pm,res);small=log2(2*2**DS.log2li(res['bins'][0][1]))
    out=dict(two_log=res['two_log_units'],one_log_instr_log2=log2(res['one_log_instr']),boot_own_log2=bt['own_log2'],boot_sub_log2=bt['sub_log2'],boot_retry=bt['retry'],boot_psi=bt['psi'],boot_phi=bt['phi'],boot_balls_log2=log2(bt['balls']),
        nodes=len(res['DET']),sieve_nodes=len(sv),mu_min=min(d['mu']for d in det),phi_max=max(d['phi']for d in det),psi_max=max(d['psi']for d in det),ret_min=min(d['ret']for d in sv),loss_max=max(max(d['loss_f'],d['loss_g'])for d in sv),
        R_max=max(d['R']for d in sv),store_bytes_log2=log2(4*max(DS.CELLS[d['R']]for d in sv)),score_max=16*max(max(d['nf'],d['ng'])for d in sv),log2_EN=log2(EN),abort_log2=log2(EN)-N_MAX_LOG2,repeat_log2=2*N_MAX_LOG2-1+6-small,max_node_units_log2=max(log2(w/DS.UNIT)for w in res['W'].values()))
    out['ok']=out['mu_min']>=Pm.MIN_RELS and out['phi_max']<=.05 and max(out['psi_max'],bt['psi'])<=Pm.PSI_MAX and out['score_max']<2**16 and out['R_max']<=128 and all(d['R']>=8 for d in sv)
    return out,res
def recurrence():
    """[7r] fast: the work/failure recurrence recomputed from the stored per-node (own, mu, window) records (no integration)."""
    st=stats();Pm=DS.P();nd={(tuple(r[0]),r[1]):r for r in st['nodes']};bins=DS.bins_of(Pm);W={};PHI={}
    wmax=lambda s,n:max([W[b,s]for b in bins if b[0]<n-1e-9]+[.0]);phimax=lambda s,n:max([PHI[b,s]for b in bins if b[0]<n-1e-9]+[.0])
    for b in bins:
        for side in'fg':
            _,_,own_fixed,surv,mu,miss_ch=nd[b,side];n,iA,kA,kB,R,thin=DS.PLAN16[b[0],side];h=(n-Pm.LX)/Pm.MESH if n>Pm.LX else.0;nA=Pm.LX+iA*h
            low=[k for k in W if k[0][0]<=Pm.THIN_MAX];pm=max([PHI[k]for k in low]+[.0]);wm=max([W[k]for k in low]+[.0])
            psi=kA*(phimax('f',nA)+phimax('g',nA))+kB*(phimax('f',n)+phimax('g',n))+miss_ch*pm;sub=kA*(wmax('f',nA)+wmax('g',nA))+kB*(wmax('f',n)+wmax('g',n))+miss_ch*wm
            window=(kA if nA>Pm.LX+1e-9 else 0)+(kB if n>nA+1e-9 else 0);r=Pm.NB_R1 if window<=1 else Pm.NB_R2
            W[b,side]=own_fixed+surv*DS.lazy_fraction(mu,r,psi)+sub/(1-psi);PHI[b,side]=DS.phi_nb(mu*(1-psi),Pm,window)
    bo=st['numbers']['boot_record'];n1,iA,kA,Rb=DS.BOOT16;nA=Pm.LX+iA*(n1-Pm.LX)/Pm.MESH;low=[k for k in W if k[0][0]<=Pm.THIN_MAX];pm=max(PHI[k]for k in low);wm=max(W[k]for k in low)
    psi=kA*phimax('g',nA)+phimax('g',float(n1))+bo['miss']*pm;balls=(1-psi**Pm.BOOT_RELS)/(1-psi)*(bo['per_rel']/DS.CELLS[Rb]+Pm.BALL_PAIR);sub=kA*wmax('g',nA)+wmax('g',float(n1))+bo['miss']*wm
    one=Pm.LLL_NODE*balls+balls*DS.CELLS[Rb]*bo['trial']+sub/(1-psi);return log2(2*one/DS.UNIT),psi**Pm.BOOT_RELS,max(PHI.values())
def probes(keys=((56,'f'),(56,'g'),(102,'f'),(197,'g'))):
    """[7p] single nodes recomputed by the model (integration, retention) against the stored records."""
    st=stats();Pm=DS.P();nd={(tuple(r[0]),r[1]):r for r in st['nodes']};ok=True;out=[]
    for b in DS.bins_of(Pm):
        for side in'fg':
            if(b[0],side)in keys:
                d=DS.node_cost(Pm,b[0],b[1],side,float(DS.PLAN16[b[0],side][0]),*DS.PLAN16[b[0],side][1:5],bool(DS.PLAN16[b[0],side][5]));r=nd[b,side]
                good=math.isclose(d['mu'],r[4],rel_tol=1e-9)and math.isclose(DS.P.LLL_NODE+2**d['sieve_log2'],r[2],rel_tol=1e-9)and math.isclose(2**d['survcost_log2'],r[3],rel_tol=1e-9);ok&=good;out.append((b[0],side,round(d['mu'],2),round(d['own_log2'],3),round(d['ret'],4),good))
    return ok,out
def laws():
    """[11b] omitted-mass law summaries recomputed from the local structure table; residue-class samples at 2, 3, 5 against the local bounds."""
    st=stats();stored=st['omitted_law'];RT.LAWS=None;RT.omitted_pdf.cache_clear();ok=True;out={}
    for sd in'fg':
        got=RT.law(sd);ok&=got==stored[sd];pdf,_=RT.omitted_pdf(sd,True);out[sd]=(round(sum(x*v for x,v in pdf.items()if x<RT.SENT)/RT.GRID,2),RT.cell_loss(sd,96))
    RT.LAWS=stored;return ok,out
def residues(seed=1,n=300):
    import descent_witness as W;st=stats();tab=st['prime_structure']['table'];ok=True;out={}
    for p in(2,3,5,7,13):
        rs=W.residue_valuations(p,n,seed+p)
        for side in'fg':
            bound=1-RT.local_pdf(p,tab[str(p)],side,False,True).get(0,.0);frac=1-rs[side].get(0,0)/rs['n'];good=frac<=bound+3*math.sqrt(bound*(1-bound)/rs['n'])+1e-12;ok&=good;out[f"{p}{side}"]=(round(frac,4),round(bound,4))
    return ok,out
def structure_spot(ps=(2,3,5,7,13,29,31,257)):
    st=stats();tab=st['prime_structure']['table'];return all(RT.structure(p)==tab[str(p)]for p in ps)
def pools(seed=1):
    """[11] one fresh special-q lattice per side and one fresh certified boot lattice against the model bins."""
    import descent_witness as W;Pm=DS.P;cfmax=max(b[1]for b in Pm.QBINS);cgmax=max(b[2]for b in Pm.QBINS);ok=True;out=[]
    for bits,side in((90,'f'),(150,'g')):
        r=W.specialq_ball_stats(bits,side,seed,points=12,radii=(16,96))
        for s in r['stats']:ok&=s['cf_max']<=cfmax and s['cg_max']<=cgmax;out.append((bits,side,s['R'],round(s['cf_max'],1),round(s['cg_max'],1)))
    rng=random.Random(seed);z=([rng.randrange(W.P)for _ in range(6)],[rng.randrange(W.P)for _ in range(6)]);z=W.pow12(z,(W.P**12-1)//W.R);B,_=W.reduce(W.boot_basis(z),W.QB);W.certify_boot(z,B)
    for s in W.boot_ball_stats(z,B,seed,points=24,radii=(4,8)):ok&=s['cg_max']<=max(x[1]for x in Pm.BOOT_BINS);out.append(('boot',s['R'],round(s['cg_mean'],1),round(s['cg_max'],1)))
    return ok,out
def campaigns():
    """[7c] stored campaign records: certificate ratios (DISC) and unit duplicates (UDUP)."""
    st=stats();c=st['pools'];Pm=DS.P;n=c['cert']['n'];ok=c['cert']['below_one']==0 and 1-0.05**(1/n)<=Pm.DISC and c['udup']['boot_yr_lo']>=Pm.UDUP and c['udup']['node_yr_min']>=Pm.UDUP;return ok,dict(cert_n=n,cert_min_log2=c['cert']['min_log2'],bound95=1-0.05**(1/n),udup=c['udup'])
