"""Run all finite witnesses and exact accounting for the final candidate bytes.
An executable certificate of identities/counts; not a proof of H_T or real hardware.
"""
import contextlib,hashlib,io,json,os,runpy,struct,subprocess,sys,tempfile
from pathlib import Path
root=Path(__file__).resolve().parent
if(root/'files').is_dir():root=root/'files'
sys.path.insert(0,str(root));sys.dont_write_bytecode=True
out=io.StringIO();saved=sys.stdin;sys.stdin=io.StringIO('{"seed":20261105}')
with contextlib.redirect_stdout(out):runpy.run_path(str(root/'probe.py'),run_name='__main__')
sys.stdin=saved;d=json.loads(out.getvalue());assert d['all_checks']
L=d['ledger'];assert L['node_bits']==72 and L['side_bits']==36 and L['groups']==1600 and L['L']==69250000
assert L['partial_key_bits']==45 and L['id_recovery_instructions']==1<<22
assert d['source_recovery']['all_checks'] and d['metadata_boundary'] and d['phase0']['all_checks']
assert d['fixed_layout']['nonoverlapping'] and d['fixed_layout']['peak_exp']==136
import a64,flip_witness as fw
src=(root/'p256_fields.s').read_text()
full=src+(root/'p256_inv.s').read_text()+fw.fwd_listing(src)+fw.bwd_listing(src)+fw.ychain_f_listing(src)+fw.ychain_b_listing(src)+fw.yprep_listing(src)+fw.tbld_listing()+fw.group_listing()
inputs={p.name:p.read_text() for p in root.glob('*.s')}
inputs['expanded.s']=full;inputs['expanded-block.s']=fw.block_production(root)
assembly={}
def text_size(raw):
 assert raw[:4]==b'\x7fELF' and raw[4]==2 and raw[5]==1
 shoff=struct.unpack_from('<Q',raw,40)[0];ents,num,idx=struct.unpack_from('<HHH',raw,58)
 heads=[struct.unpack_from('<IIQQQQIIQQ',raw,shoff+ents*k) for k in range(num)]
 h=heads[idx];names=raw[h[4]:h[4]+h[5]]
 return next(h[5] for h in heads if names[h[0]:].split(b'\0',1)[0]==b'.text')
# Assembly is a separate optional native check; hosted Python experiments do not invoke this certificate.
with tempfile.TemporaryDirectory(prefix='p256-grid72-') as td:
 for name,s in inputs.items():
  pp=Path(td)/name;obj=pp.with_suffix('.o');pp.write_text(s)
  r=subprocess.run(['clang','-target','aarch64-linux-gnu','-march=armv8.7-a+ls64+lse','-c',str(pp),'-o',str(obj)],capture_output=True,text=True)
  assert r.returncode==0,(name,r.stderr)
  n=len(a64.a64_load(s)[0]);size=text_size(obj.read_bytes());assert size==4*n,(name,size,n)
  assembly[name]={'instructions':n,'text_bytes':size}
assert assembly['expanded.s']['instructions']==7782
w=(root/'wolfram.txt').read_text();assert 'Exit code: 0' in w and '7197421090497760703, True, True' in w and '475392, True' in w
proof=(root/'proof.md').read_text();assert 'home70 bits' not in proof and '*2^128/(r-1)' not in proof
# Identity inventory for reproducible final bytes, excluding generated caches.
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(root.iterdir()) if p.is_file()}
cn=int(L['claim'].split('/')[0])
print(json.dumps({'all_checks':True,'zero_failed_checks':True,'score':L['claim'],'score_decimal':'%d.%015d'%(cn//10**15,cn%10**15),'instructions':L['E_total_instructions'],'gate_k':L['gate_k'],'run_success':L['success_formula'],'source_recovery':d['source_recovery'],'layout':d['fixed_layout'],'assembly':assembly,'manifest_matches_ledger':d['manifest_matches_ledger'],'hashes':hashes,'scope':'Exact accounting, finite arithmetic/protocol checks and native encodability. H_T(0.35) and ideal-machine premises remain explicit.'}))
