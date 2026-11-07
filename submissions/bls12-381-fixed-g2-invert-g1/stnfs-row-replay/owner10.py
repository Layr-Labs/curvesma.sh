"""Seven active vector owners (N < 7 x 2^56 at merge weight 300); blocks (72, 5)."""
import math
ACTIVE_OWNERS=7
C_ROW=1350.0
LA_BOUND=116.036935623
def row_tail(send=15,combine=62,finish=92,fold9=180,reduce6=183,owners=ACTIVE_OWNERS):
	return(owners-1)*(send+18)+(owners-1)*(combine+2)+finish+fold9+reduce6+2+10*owners+40+28*owners
def install(model):
	assert model.C_ROW==2710.0
	assert(model.BW_M,model.BW_N)==(72,5)
	model.C_ROW=C_ROW
def check(model,result):
	owners=math.ceil(2**result['log2_N']/(1<<56))
	row=row_tail()
	return(
	owners<=ACTIVE_OWNERS
	and 2**result['log2_N']<7*(1<<56)
	and row==1345
	and row<=model.C_ROW==C_ROW
	and(model.BW_M,model.BW_N)==(72,5)
	and result['log2_la']<=LA_BOUND
	)

# Lossless packaging only: verify the full archived source before executing it.
def _load_source(namespace,name,archive='stats16.b85'):
 import base64,hashlib,json,lzma
 from pathlib import Path
 here=Path(namespace['__file__']).resolve().parent
 root=here/'files'if(here/'files').is_dir()else here
 entry=json.loads(lzma.decompress(base64.b85decode((root/archive).read_bytes())))['source_archive'][name]
 source=entry['source'];assert hashlib.sha256(source.encode()).hexdigest()==entry['sha256']
 exec(compile(source,namespace['__file__'],'exec'),namespace)
