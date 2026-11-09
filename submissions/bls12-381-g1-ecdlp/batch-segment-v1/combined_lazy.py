import contextlib,io,json,runpy,sys
from pathlib import Path
root=Path(__file__).resolve().parent
if (root/'files').is_dir():root=root/'files'
sys.path.insert(0,str(root))
import kernel_check, nocarry_witness, batch_audit
buf=io.StringIO()
with contextlib.redirect_stdout(buf):runpy.run_path(str(root/'evidence.py'),run_name='__main__')
out=json.loads(buf.getvalue())
if out['mode']=='scan':
 import evidence
 out['ratio_scan']=[]
 for q in ('23/20','231/200','29/25','115917/100000','233/200','117/100','47/40'):
  z=evidence.ledger({'ratio':q})
  out['ratio_scan'].append({'q':q,'instructions':z['W_instructions'],'score_numerator':z['score_numerator']})
 assert min(out['ratio_scan'],key=lambda z:z['instructions'])['q']=='115917/100000'
 
 print(json.dumps(out));sys.exit(0)
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
out['h7_lanes']=batch_audit.xge_run({})
out['table_level']=batch_audit.run({'table':True,'seed':7})
out['nocarry_witness']=nocarry_witness.run({'seed':5,'samples':2000})
from fractions import Fraction
import log_bound
out['exact_log_bound']=log_bound.certify(out['ledger']['W_instructions'],2000,Fraction(115886487960384670,10**15))
for shift in (-1,1):
 try: log_bound.certify(out['ledger']['W_instructions'],2000,Fraction(115886487960384670+shift,10**15))
 except AssertionError: pass
 else: raise AssertionError('wrong log grid accepted')
out['log_grid_negative_controls']=2
out['scope_new']='Coordinate-only full kernels; strict log grid; finite checks.'
print(json.dumps(out))
