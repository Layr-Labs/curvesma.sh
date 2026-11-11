"""Whole-ball finite-pool retries plus independent single-cell boot.
Run only through llm-verifier. Actual finite-family positivity remains explicit H1.
"""
import json,sys,math
from pathlib import Path
H=Path(__file__).resolve().parent
sys.path.insert(0,str(H/'files' if(H/'files').is_dir()else H))
import common12
import descent16 as D,descent_checks as DC

def fixed_modulus_ecm(fbits, nbits, success):
    import ecm16 as E
    z=E.schedule(fbits,nbits,success)
    ops=1.1*(z['stage1_M']+z['stage2_M'])
    mul=z['mul_instr']; gd=z['gcd_instr']
    outputs=E.pi_upper(z['B1'])+z['giant']
    # Original-modulus curves; charge a full replay on EVERY curve, plus
    # output gcd and conservative remainder/partition refinement work.
    return z['curves']*(2*ops*mul+2*gd+5*outputs*gd)
D.ecm_instr=fixed_modulus_ecm

def evaluate(pool=None):
    pools={"f":12,"g":36} if pool is None else {"f":pool,"g":pool}
    DC.stats();P=D.P;bins=D.bins_of(P);W={};details={};node_counts={}
    def mx(s,n):return max([v for(b,t),v in W.items()if t==s and b[0]<n-1e-9]+[0.])
    for b in bins:
        for side in 'fg':
            class Shift(P):pass
            dx,dy=(25,8) if side=='f' else(32,1)
            Shift.QBINS=[(w,cf+dx,cg+dy)for w,cf,cg in P.QBINS]
            Shift.QBINS_LO=[(cf+dx,cg+dy)for cf,cg in P.QBINS_LO]
            Shift.QBINS_SUB=[(cf+dx,cg+dy)for cf,cg in P.QBINS_SUB]
            Shift.MIN_RELS=0
            n,iA,kA,kB,R,thin=D.PLAN16[b[0],side]
            d=D.node_cost(Shift,b[0],b[1],side,float(n),iA,kA,kB,R,bool(thin))
            if d is None:raise ValueError('nonpositive smoothness model')
            nA=d['nA']
            low=[v for (bb,s),v in W.items()if bb[0]<=P.THIN_MAX]
            missing=d['miss_children']*max(low+[0.])
            sub=kA*(mx('f',nA)+mx('g',nA))+kB*(mx('f',n)+mx('g',n))+missing
            # Charge EVERY potential cofactor survivor, no negative-binomial lazy discount.
            # Factor two prices simultaneous factorization forest and all bounded overhead.
            own=2**40+P.LLL_NODE+(2**d['sieve_log2'] if d['how']=='sieve' else 0)+2**d['survcost_log2']
            delta=.5;d['det']=delta;pool_side=pools[side];balls=pool_side/delta
            if delta<.5:raise ValueError('finite schedule detection below half')
            work=balls*own+sub;W[b,side]=work
            details[b,side]=dict(d,own=own,W_log2=math.log2(work),balls=balls,sub_log2=math.log2(sub)if sub else None,phi=0.,psi=0.,retry=1.,pool=pool_side,finite_population_mean=d['mu'])
    n,iA,kA,R=D.BOOT16;nA=P.LX+iA*(n-P.LX)/P.MESH;c=(P.LOG2P6+math.log2(D.CELLS[R]))/12
    m=12*c+P.LOG2U6+max(x[1]for x in P.BOOT_BINS)
    pr=sum(w*D.ss_lower2(12*c+P.LOG2U6+cg,56,n,iA,kA,1,P.MESH,P.ALPHA_G,0)for w,cg in P.BOOT_BINS)
    pri=sum(D.ss_lower2(12*c+P.LOG2U6+cg,56,n,iA,kA,1,P.MESH,P.ALPHA_G,0)for cg in P.BOOT_SUB)/len(P.BOOT_SUB)
    pr=min(pri,P.PRCAP_BOOT*pr)
    det=.5
    import retain16
    prob=pr*det*P.UDUP*retain16.counted_fraction()
    trial=2**40+P.NORMS+D.ecm_instr(nA,m,.999)+P.PRIME+P.IDEAL
    low=max([v for (bb,s),v in W.items()if bb[0]<=P.THIN_MAX]+[0.])
    sub=kA*mx('g',nA)+mx('g',n)+D.missing(P,0,m)*low
    one=trial/prob+sub
    bt=dict(n1=n,nA=nA,kA=kA,R=R,m=m,pr_log2=math.log2(prob),own_log2=math.log2(trial/prob),sub_log2=math.log2(sub),retry=1.,psi=0.,phi=0.,balls=1/prob)
    res=dict(W=W,DET=details,bins=bins,boot=bt,two_log_units=math.log2(2*one/D.UNIT),one_log_instr=one,PHI={k:0. for k in W})
    vals=list(details.values());sv=[d for d in vals if d['how']=='sieve']
    t=dict(two_log=res['two_log_units'],one_log_instr_log2=math.log2(one),boot_own_log2=bt['own_log2'],boot_sub_log2=bt['sub_log2'],boot_retry=1.,boot_psi=0.,boot_phi=0.,boot_balls_log2=math.log2(bt['balls']),nodes=len(vals),sieve_nodes=len(sv),mu_min=min(d['mu']for d in vals),phi_max=0.,psi_max=0.,ret_min=min(d['ret']for d in sv),loss_max=max(max(d['loss_f'],d['loss_g'])for d in sv),R_max=max(d['R']for d in vals),store_bytes_log2=math.log2(4*max(D.CELLS[d['R']]for d in sv)),score_max=16*max(max(d['nf'],d['ng'])for d in sv),log2_EN=0.,abort_log2=-1000.,repeat_log2=-1000.,max_node_units_log2=max(math.log2(v/D.UNIT)for v in W.values()),pool_size=pools,mean_min=min(d['mu']for d in vals))
    t['ok']=t['score_max']<65536 and t['R_max']<=128 and t['mu_min']>0
    return t,res

def main(req):
    t,res=evaluate(int(req['pool']) if 'pool' in req else None)
    t['rows']=[dict(bin=b[0],side=s,own_log2=math.log2(d['own']),work_log2=d['W_log2'],mu=d['mu'],balls=d['balls'],R=d['R'],nf=d['nf'],ng=d['ng'])for(b,s),d in res['DET'].items()]
    return t
if __name__=='__main__':print(json.dumps(main(json.load(sys.stdin))))
