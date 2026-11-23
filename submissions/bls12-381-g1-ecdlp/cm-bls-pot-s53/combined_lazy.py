import contextlib,io,json,math,runpy,sys
sys.set_int_max_str_digits(0)
from pathlib import Path
root=Path(__file__).resolve().parent
if (root/'files').is_dir():root=root/'files'
sys.path.insert(0,str(root))
import kernel_check, nocarry_witness, batch_audit
raw=sys.stdin.read();cfg=json.loads(raw);sys.stdin=io.StringIO(raw)
buf=io.StringIO()
with contextlib.redirect_stdout(buf):runpy.run_path(str(root/'evidence.py'),run_name='__main__')
out=json.loads(buf.getvalue())
import a64,gc;a64._a64_parse.cache_clear();gc.collect()
from evidence import R_PRODUCTION as r, D_PRODUCTION as d
from math import gcd
m=(r-1)//d; omega=pow(7,m,r); k=256
assert gcd(m,d)==4096 and pow(omega,d//2,r)==r-1
assert pow(omega,d,r)==1 and (m*k)%d==0
assert pow(pow(7,k,r),m,r)==pow(omega,k,r)!=1
out['character_base']='passed: log base omega, not the full-field generator'
out['kernel_check']=kernel_check.run({'seed':113,'samples':200})
out['kernel_check_plus']=kernel_check.run({'seed':113,'samples':200,'plus':True})
out['leaf_kernel']=kernel_check.run_leaf({'seed':113,'samples':200})
out['leaf_kernel_plus']=kernel_check.run_leaf({'seed':113,'samples':200,'plus':True})
out['beta_keys']=batch_audit.beta_run({'seed':11})
out['table_level']=batch_audit.run({'table':True,'seed':7})
out['nocarry_witness']=nocarry_witness.run({'seed':5,'samples':2000})
if 'ledger' in out:
 from fractions import Fraction
 import log_bound
 L=out['ledger'];S=L['score_numerator'];E=L['expected_numerator_x2pow64'];R=L['run_instructions']
 assert E==64*2**64+L['gate_K']*L['expected_run_instructions'] and R==L['run_total_instructions']+64
 out['exact_log_bound']=log_bound.certify(E,2000*2**64,Fraction(S,10**12))
 out['opened_run_log_bound']=log_bound.certify(R,2000,Fraction(L['run_score_numerator'],10**12))
 for shift in (-1,1):
  try: log_bound.certify(E,2000*2**64,Fraction(S+shift,10**12))
  except AssertionError: pass
  else: raise AssertionError('wrong log grid accepted')
 out['log_grid_negative_controls']=2
 if 'KT' in cfg or 'ratio' in cfg: out['manifest_check']='skipped: input overrides KT or ratio'
 else:
  # the manifest's declared figures must be the ledger's; every quoted hypothesis bound must hold and stand in the text
  M=json.loads((root/'submission.json').read_text());c=M['claims'][0];st=c['statement']
  g=lambda n:'%d.%012d'%divmod(n,10**12);lo=Fraction(L['run_success_lower'])*10**60;assert lo.denominator==1
  dec={'score':'W(n)<=2^%s G1'%g(S),'every_run':'at most 2^%s units'%g(L['run_score_numerator']),'gate_K':'K=%d is'%L['gate_K'],
   'lo_run':'lo_run=0.%s being'%str(lo.numerator).zfill(60),'K_T':'K_T = %s ('%L['KT'],'E_run':'E_run=%d instructions'%L['expected_run_instructions'],
   'Q':'Q=%s '%L['Q'],'preprocessing':'preprocessing<2^%.4f units'%L['preprocessing_log2_ceil_1e4'],'memory':'memory<2^%.4f bytes'%L['memory_log2_ceil_1e4']}
  assert Fraction(int(c['score']['numerator']),int(c['score']['denominator']))==Fraction(S,10**12) and g(S) in M['title'],'declared score'
  for k_,v in dec.items(): assert st.count(v)==1,k_
  T=(root/'source-proof.md').read_text();TP=(root/'proof.md').read_text()
  up=lambda x,n:'%.*f'%(n,math.ceil(x*10**n)/10**n)
  fm,pb,pg=L['false_match_mean_log2'],*L['pool_mean_log2'];Kf=Fraction(L['K']);nu=L['h5_conditional_mean_nu']
  H=[('H3 pools','means below 2^%s and 2^%s against 2^56 and 2^57'%(up(pb,4),up(pg,4)),float(up(pb,4))<L['pool_threshold_log2'][0] and float(up(pg,4))<L['pool_threshold_log2'][1]),
   ('H3 false matches','3 N_b N_g/2^192 < 2^%s against the 2^45 cap'%up(fm,2),float(up(fm,2))<45),
   ('H3 owner candidates','above 2^16 (mean 2^%s)'%up(fm-62,2),fm-62<-17),
   ('H3 candidate links','(mean below 2^%s)'%up(fm-82,1),fm-82<-37),
   ('H5 total mean','N_b N_g/M < 2^%s + %s (proof.md s43)'%(up(fm,2),up(Kf,2)),Kf>=Fraction(L['KT'])),
   ('H5 conditional mean','<= %s at the claimed coverage'%up(nu,2),nu<8),
   ('H2 bases','Bases 2^%s in all, below 2^63'%up(L['bases_log2'],2),L['bases_log2']<63),
   ('H2 prime margin','primes: margin above %.2f'%(math.floor(L['prime_margin']*100)/100),L['prime_margin']>3),
   ('slot load','3 N_g/2^119 < %s'%up(L['slot_load'],4),L['slot_load']<1)]
  for k_,q,ok in H: assert ok and (TP.count(q)==2 if 'margin' in k_ else q in (TP if k_[:2]=='H2' or k_[0]=='s' else T)),k_
  out['manifest_check']={'declared_equals_ledger':sorted(dec),'score':g(S),'every_run':g(L['run_score_numerator'])}
  out['hypothesis_bounds']=[[k_,q] for k_,q,ok in H]
import table_trade_audit
out['table_trade_audit']=table_trade_audit.run()
import range_audit
out['range_audit']=range_audit.run()
out['scope_new']='s45 interleaved rounds (owner bodies executed by audit_il, per-round comparison budgets, two-sided expected work), s46 row products (batch round and kernels executed with rows), s47 table reduction (leaf and batch inverse executed against exact models), s41 manifest and hypothesis-bound checks, s49 zero-high-word and raw-inverse consumer range audits; finite checks.'
import init_audit
out['owner_initialization_audit']=init_audit.run()
if 'ledger' in out:
 assert out['ledger']['rows_instructions']['owner initialisation']==(1<<62)*init_audit.production_cost(1<<61)
 if 'KT' not in cfg and 'ratio' not in cfg:
  assert 'DCZID_EL0.BS=6,DZP=0' in st and 'C=2000' in st
out['scope_new']+=' s50 declared 256-byte zero-block initialization/count/capability audit; s53 1-temp SQR with live out/buf, branch-canonical full-child and beta kernels, and post-indexed batch-round updates executed word-for-word by a64.py.'
print(json.dumps(out))
