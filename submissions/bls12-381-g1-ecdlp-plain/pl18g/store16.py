"";VX=30;V6=16;VV=60
def VY(cell,chunk_bits=VX):assert 0<=cell<1<<VV;return cell>>chunk_bits,cell&(1<<chunk_bits)-1
def route(VU,chunk_bits=VX,radix_bits=V6):
	out=list(VU)
	for shift in range(0,VV-chunk_bits,radix_bits):
		V4={}
		for VW in out:chunk,_=VY(VW[0],chunk_bits);digit=chunk>>shift&(1<<radix_bits)-1;V4.setdefault(digit,[]).append(VW)
		out=[VW for digit in sorted(V4)for VW in V4[digit]]
	return out
def scatter(VU,chunk_bits=VX):
	tails={}
	for VW in VU:tails.setdefault(VY(VW[0],chunk_bits)[0],[]).append(VW)
	return[VW for chunk in sorted(tails)for VW in tails[chunk]]
def V3(VU,V0,cells,chunk_bits):
	V5=[0]*cells
	for cell,ideal,weight in VU:V5[cell]+=weight
	routed=route(VU,chunk_bits,4);chunks={}
	for cell,ideal,weight in routed:ch,off=VY(cell,chunk_bits);chunks.setdefault(ch,[]).append((off,ideal,weight))
	got=[0]*cells;recovered=[];peak=0
	for ch,rows in chunks.items():
		V2=[0]*min(1<<chunk_bits,cells-(ch<<chunk_bits));peak=max(peak,len(V2))
		for off,ideal,weight in rows:V2[off]+=weight
		for off,v in enumerate(V2):got[(ch<<chunk_bits)+off]=v
		recovered.extend(((ch<<chunk_bits)+off,ideal,weight)for off,ideal,weight in rows if V2[off]>=V0)
	assert got==V5,'chunked apply differs from the full-array oracle';assert sorted(recovered)==sorted(r for r in VU if V5[r[0]]>=V0),'recovery differs';assert peak<=1<<chunk_bits;return dict(chunks=len(chunks),peak_cells=peak,recovered=len(recovered))
def VZ():
	import math;from fractions import Fraction;from ideal_canonical import determinant;import common12;cfg=common12.stats()['config'];m=[[Fraction(x)for x in row]for row in cfg['matrix_b']];den=1
	for row in m:
		for x in row:den=den*x.denominator//math.gcd(den,x.denominator)
	M=[[int(x*den)for x in row]for row in m];d=determinant(M);rad=Fraction(cfg['production_region_lam'])*Fraction(103,100)*cfg['production_A_squared']*den;return math.log2(math.prod(2*math.sqrt(float(rad*determinant([[M[a][b]for b in range(6)if b!=i]for a in range(6)if a!=i])/d))+1 for i in range(5)))
def check():
	import random,math;rng=random.Random(761);V1=VZ();assert V1<VV,'a line box exceeds 2^60 cells';VU=[(rng.randrange(1<<VV),rng.randrange(1<<63),rng.randrange(1,256))for _ in range(3000)];VU+=[(c,1,17)for c in(0,(1<<30)-1,1<<30,(1<<60)-1)];V7=route(VU);assert V7==sorted(VU,key=lambda r:r[0]>>VX),'routing is not the stable chunk order';assert scatter(VU)==V7,'scatter differs from the radix route';assert max(VY(r[0])[0]for r in VU)<1<<30;toy=[(rng.randrange(1031),rng.randrange(100),rng.randrange(1,256))for _ in range(4000)];eq=V3(toy,600,1031,5);half=6863095027560050329435265332068;assert half<1<<103 and half-1>>40<1<<63
	for idx in 0,(1<<40)-1,1<<40,(1<<94)-1,1<<94,half-1:own,off=idx>>40,idx&(1<<40)-1;assert own<1<<64 and off<1<<40 and own<<40|off==idx
	assert 16*69<1<<15 and 16*2800<1<<16;n_model=2.**59.7306;blocks=math.ceil(n_model/(1<<56));assert blocks==14 and 20*(1<<61)>=n_model*32
	return dict(line_cells_log2_max=round(V1,3),routing_stable=True,production_boundaries=[dict(cell=str(c),chunk=VY(c)[0],offset=VY(c)[1])for c in(0,(1<<30)-1,1<<30,(1<<60)-1)],toy_chunked_pass=eq,descent_owner_bits=round(math.log2(half-1>>40),2),la_blocks_required=blocks,la_blocks=20)
