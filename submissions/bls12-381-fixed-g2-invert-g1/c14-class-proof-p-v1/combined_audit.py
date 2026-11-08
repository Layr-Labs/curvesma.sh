"""Single declared run: inherited attack audit plus fixed-field certificate."""
import contextlib,io,json,runpy,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
if (ROOT/'files').is_dir():ROOT=ROOT/'files'
sys.path.insert(0,str(ROOT))
def run(name,request):
    saved=sys.stdin;out=io.StringIO()
    try:
        sys.stdin=io.StringIO(json.dumps(request))
        with contextlib.redirect_stdout(out):namespace=runpy.run_path(str(ROOT/name),run_name='__main__')
    finally:sys.stdin=saved
    return namespace,json.loads(out.getvalue())
request=json.load(sys.stdin)
_,attack=run('audit_refinement.py',request)
if isinstance(attack["descent"],str):attack["descent"]=attack["descent"].replace('{"descent": true}','{"descent": "tree"}')
certificate_globals,certificate=run('class_checker.py',{})
import common12 as C,descent_witness as W
basis=certificate_globals['CERTIFICATE']['basis']
binding=(C.H==certificate_globals['H'] and C.U==certificate_globals['U']
         and C.WP==[1,1,0,0,0,0] and C.BAS2==basis and W.R==certificate_globals['R'])
F=certificate_globals['F']
expected_f=[(2*i,[F[i]*v for v in C._pw(C.W,6-i)]) for i in range(6,-1,-1) if F[i]]
binding=binding and C.TERMS==expected_f and C.UW==[C.U*v for v in C.W]
assert binding and attack['checks_ok'] and certificate['status']=='ok' and certificate['k_class_number']==1
print(json.dumps({'checks_ok':True,'field_binding':binding,'attack_audit':attack,
                 'class_number_certificate':certificate,
                 'scope':'Unchanged attack bound; only the class-number part of H3 is certified.'},sort_keys=True))
