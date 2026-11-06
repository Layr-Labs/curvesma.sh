"""Executable reference of the chunked record protocol (STORAGE.md, where the"""
CHUNK_BITS=30
RADIX_BITS=16
MAX_CELL_BITS=59
def address(cell,chunk_bits=CHUNK_BITS):
	assert 0<=cell<1<<MAX_CELL_BITS
	return cell>>chunk_bits,cell& ((1<<chunk_bits) -1)
def route(records,chunk_bits=CHUNK_BITS,radix_bits=RADIX_BITS):
	out=list(records)
	for shift in range(0,MAX_CELL_BITS-chunk_bits,radix_bits):
		buckets= {}
		for record in out:
			chunk,_=address(record[0],chunk_bits)
			digit= (chunk>>shift) & ((1<<radix_bits) -1)
			buckets.setdefault(digit, []).append(record)
		out= [record for digit in sorted(buckets)for record in buckets[digit]]
	return out
def scatter(records,chunk_bits=CHUNK_BITS):
	tails= {}
	for record in records:
		tails.setdefault(address(record[0],chunk_bits)[0], []).append(record)
	return[record for chunk in sorted(tails)for record in tails[chunk]]
def chunked_pass(records,threshold,cells,chunk_bits):
	oracle= [0] *cells
	for cell,ideal,weight in records:
		oracle[cell] +=weight
	routed=route(records,chunk_bits,4)
	chunks= {}
	for cell,ideal,weight in routed:
		ch,off=address(cell,chunk_bits)
		chunks.setdefault(ch, []).append((off,ideal,weight))
	got= [0] *cells;recovered= [];peak=0
	for ch,rows in chunks.items():
		scores= [0] *min(1<<chunk_bits,cells- (ch<<chunk_bits));peak=max(peak,len(scores))
		for off,ideal,weight in rows:
			scores[off] +=weight
		for off,v in enumerate(scores):
			got[(ch<<chunk_bits) +off] =v
		recovered.extend(((ch<<chunk_bits) +off,ideal,weight)for off,ideal,weight in rows if scores[off] >=threshold)
	assert got==oracle,'chunked apply differs from the full-array oracle'
	assert sorted(recovered) ==sorted(r for r in records if oracle[r[0]] >=threshold),'recovery differs'
	assert peak<=1<<chunk_bits
	return dict(chunks=len(chunks),peak_cells=peak,recovered=len(recovered))
def line_cells_bound():
	import json,math,os
	from fractions import Fraction
	from ideal_canonical import determinant
	cfg=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'config.json')))
	m= [[Fraction(x)for x in row]for row in cfg['matrix_b']]
	den=1
	for row in m:
		for x in row:
			den=den*x.denominator//math.gcd(den,x.denominator)
	det5=Fraction(determinant([[int(m[i][j] *den)for j in range(5)]for i in range(5)]),den**5)
	a2=cfg['production_A_squared']
	vol=8*math.pi**2/15*a2**2.5/math.sqrt(float(det5))
	a_min=math.sqrt(a2/float(sum(m[i][i]for i in range(5))))
	return math.log2(vol* (1+math.sqrt(5) / (2*a_min)) **5)
def check():
	import random,math
	rng=random.Random(761)
	line_bits=line_cells_bound()
	assert line_bits<MAX_CELL_BITS,'a line exceeds 2^59 cells'
	records= [(rng.randrange(1<<MAX_CELL_BITS),rng.randrange(1<<63),rng.randrange(1,256))for _ in range(3000)]
	records+= [(c,1,17)for c in(0, (1<<30) -1,1<<30, (1<<59) -1)]
	ordered=route(records)
	assert ordered==sorted(records,key=lambda r:r[0] >>CHUNK_BITS),'routing is not the stable chunk order'
	assert scatter(records) ==ordered,'scatter differs from the radix route'
	assert max(address(r[0])[0]for r in records) <1<<29
	toy= [(rng.randrange(1031),rng.randrange(100),rng.randrange(1,256))for _ in range(4000)]
	eq=chunked_pass(toy,600,1031,5)
	half=6863095027560050329435265332068
	assert half<1<<103 and(half-1) >>40<1<<63
	for idx in(0, (1<<40) -1,1<<40, (1<<94) -1,1<<94,half-1):
		own,off=idx>>40,idx& ((1<<40) -1)
		assert own<1<<64 and off<1<<40 and(own<<40) |off==idx
	assert 16*69<1<<15 and 16*2800<1<<16
	la=la_check()
	return dict(line_cells_log2_max=round(line_bits,3),routing_stable=True,toy_chunked_pass=eq,
	            descent_owner_bits=round(math.log2((half-1)>>40),2),linear_algebra=la)

def la_address(index,next_buffer=False):
	assert 0<=index<13*(1<<55)
	owner,local=index>>55,index&((1<<55)-1)
	base=(1<<20)+((3<<61) if next_buffer else 0)
	return owner,base+192*local

def la_check():
	capacity=13*(1<<55);size=3<<61;base=1<<20
	for i in(0,(1<<55)-1,1<<55,(12<<55)-1,12<<55,capacity-1):
		owner,a=la_address(i);owner2,b=la_address(i,True)
		assert owner==owner2<13 and 0<=a-base<size and b-a==size
		assert b+192<=base+2*size<1<<64
	assert (1<<64)-(base+2*size)>1<<61
	for ordinal in(0,(1<<40)-1,1<<40,(1<<68)-1):
		shard,offset=ordinal>>40,(ordinal&((1<<40)-1))*8
		assert shard<1<<28 and offset+8<=1<<43
	for row in(0,(1<<40)-1,1<<40,(1<<59)-1):
		shard,local=row>>40,row&((1<<40)-1)
		assert shard<1<<19 and local*32+32<=1<<45 and local*1024+1024<=1<<50
	for local in(0,1,(1<<54)-1,1<<54,(1<<55)-1):
		for coefficient in(1,2,127,128,253,255):
			word=(coefficient<<55)|local
			assert word<1<<63 and (word&((1<<55)-1))==local and word>>55==coefficient
	assert (1<<34)*1024 <1<<61
	return dict(owners=13,block_bits=55,stride=192,buffers=2,local_buffers_bytes=str(2*size),
	            local_reserved_lower=str((1<<64)-(base+2*size)),entry_shard_bits=40,entry_record_bytes=8)
