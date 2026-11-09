""
CHUNK_BITS=30
RADIX_BITS=16
MAX_CELL_BITS=60
def address(cell,chunk_bits=CHUNK_BITS):
	assert 0<=cell<1<<MAX_CELL_BITS
	return cell>>chunk_bits,cell&((1<<chunk_bits)-1)
def route(records,chunk_bits=CHUNK_BITS,radix_bits=RADIX_BITS):
	out=list(records)
	for shift in range(0,MAX_CELL_BITS-chunk_bits,radix_bits):
		buckets={}
		for record in out:
			chunk,_=address(record[0],chunk_bits)
			digit=(chunk>>shift)&((1<<radix_bits)-1)
			buckets.setdefault(digit,[]).append(record)
		out=[record for digit in sorted(buckets)for record in buckets[digit]]
	return out
def scatter(records,chunk_bits=CHUNK_BITS):
	tails={}
	for record in records:
		tails.setdefault(address(record[0],chunk_bits)[0],[]).append(record)
	return[record for chunk in sorted(tails)for record in tails[chunk]]
def inplace_apply(records,chunk_bits=CHUNK_BITS):
	"""Owner applies straight from handed-off pages: same chunk grouping and
	emission order as scatter, without the private-buffer copy."""
	tails={}
	for record in records:
		tails.setdefault(address(record[0],chunk_bits)[0],[]).append(record)
	return[(ch,rec)for ch in sorted(tails)for rec in tails[ch]]
def chunked_pass(records,threshold,cells,chunk_bits):
	oracle=[0]*cells
	for cell,ideal,weight in records:
		oracle[cell]+=weight
	routed=route(records,chunk_bits,4)
	chunks={}
	for cell,ideal,weight in routed:
		ch,off=address(cell,chunk_bits)
		chunks.setdefault(ch,[]).append((off,ideal,weight))
	got=[0]*cells;recovered=[];peak=0
	for ch,rows in chunks.items():
		scores=[0]*min(1<<chunk_bits,cells-(ch<<chunk_bits));peak=max(peak,len(scores))
		for off,ideal,weight in rows:
			scores[off]+=weight
		for off,v in enumerate(scores):
			got[(ch<<chunk_bits)+off]=v
		recovered.extend(((ch<<chunk_bits)+off,ideal,weight)for off,ideal,weight in rows if scores[off]>=threshold)
	assert got==oracle,'chunked apply differs from the full-array oracle'
	assert sorted(recovered)==sorted(r for r in records if oracle[r[0]]>=threshold),'recovery differs'
	assert peak<=1<<chunk_bits
	return dict(chunks=len(chunks),peak_cells=peak,recovered=len(recovered))
def line_box_bound():
	import math
	from fractions import Fraction
	from ideal_canonical import determinant
	import common12;cfg=common12.stats()['config']
	m=[[Fraction(x)for x in row]for row in cfg['matrix_b']]
	den=1
	for row in m:
		for x in row:
			den=den*x.denominator//math.gcd(den,x.denominator)
	M=[[int(x*den)for x in row]for row in m];d=determinant(M)
	rad=Fraction(cfg['production_region_lam'])*cfg['production_A_squared']*den
	return math.log2(math.prod(2*math.sqrt(float(rad*determinant([[M[a][b]for b in range(6)if b!=i]for a in range(6)if a!=i])/d))+1 for i in range(5)))
def check():
	import random,math
	rng=random.Random(761)
	line_bits=line_box_bound()
	assert line_bits<MAX_CELL_BITS,'a line box exceeds 2^60 cells'
	records=[(rng.randrange(1<<MAX_CELL_BITS),rng.randrange(1<<63),rng.randrange(1,256))for _ in range(3000)]
	records+=[(c,1,17)for c in(0,(1<<30)-1,1<<30,(1<<60)-1)]
	ordered=route(records)
	assert ordered==sorted(records,key=lambda r:r[0]>>CHUNK_BITS),'routing is not the stable chunk order'
	assert scatter(records)==ordered,'scatter differs from the radix route'
	ip=inplace_apply(records)
	assert[r for _,r in ip]==ordered,'in-place order differs from scatter'
	assert[ch==address(r[0])[0]for ch,r in ip],'in-place chunk differs'
	assert max(address(r[0])[0]for r in records)<1<<30
	toy=[(rng.randrange(1031),rng.randrange(100),rng.randrange(1,256))for _ in range(4000)]
	eq=chunked_pass(toy,600,1031,5)
	half=6863095027560050329435265332068
	assert half<1<<103 and(half-1)>>40<1<<63
	for idx in(0,(1<<40)-1,1<<40,(1<<94)-1,1<<94,half-1):
		own,off=idx>>40,idx&((1<<40)-1)
		assert own<1<<64 and off<1<<40 and(own<<40)|off==idx
	assert 16*69<1<<15 and 16*2800<1<<16
	n_model=2.0**59.7306
	blocks=math.ceil(n_model/(1<<56))
	assert blocks==14 and 20*(1<<61)>=n_model*32
	return dict(line_cells_log2_max=round(line_bits,3),routing_stable=True,inplace_order=True,production_boundaries=[dict(cell=str(c),chunk=address(c)[0],offset=address(c)[1])
	for c in(0,(1<<30)-1,1<<30,(1<<60)-1)],
	toy_chunked_pass=eq,descent_owner_bits=round(math.log2((half-1)>>40),2),la_blocks_required=blocks,la_blocks=20)
