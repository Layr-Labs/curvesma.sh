import json,os,math,random
from math import log2,ceil,floor
import descent16 as DS,retain16 as RT
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
    import hybrid_model
    return hybrid_model.evaluate()
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
