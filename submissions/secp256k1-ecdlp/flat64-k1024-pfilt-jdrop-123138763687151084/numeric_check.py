"""Independent definition checks for Card66; no production-scale execution."""
import json, pathlib, hashlib, sys
HERE = next(p for p in (pathlib.Path(__file__).resolve().parent, pathlib.Path("/snapshot/files"), pathlib.Path.cwd()) if (p/"certificate.py").is_file())
sys.path.insert(0,str(HERE))
from fractions import Fraction as F
import certificate as c
assert not c.FAILED, c.FAILED
# Definition-bound schedule inequalities, independent of float seed/root loops.
D=22000**17; T=(256*106000000000)**32
for i,m in enumerate(c.MQ1):
 target=T*i**17
 assert m**32*D<=target<(m+1)**32*D
qq=[b-a for a,b in zip(c.MQ1,c.MQ1[1:])]
assert qq==c.Q1G and len(qq)==22000 and sum(qq)==27136000000000
assert min(qq)==655279709 and max(qq)==133854748188
assert max(qq)<1024*(1<<27) and min(qq)>1024*(1<<19)
assert c.BANK_SLOTS==1<<27 and c.N_BANK==(1<<48)*1024
assert all((q+(1<<27)-1)//(1<<27)<=1024 for q in qq)
# Unchanged proof anchors and signed cap gaps are independently recomputed.
home_mean=c.S_MINE_TOT*F(6*(1<<174),c.r-1)
home_threshold=c.CAP_H//32736
home_exponent=F(82)+F(3,2)*home_mean-home_threshold
assert home_exponent < -200 and home_exponent==c._MINE_HOME_EXPONENT
assert c.P_SAME_A_ALL<F(1,1<<39) and c.D_XMINE_ALL<F(1,1<<140) and c.P_MINE_ZERO_ALL<F(1,1<<140)
assert c.TFULL["cs"]<F(1,1<<82) and c.D_DROP<F(1,1<<26)
assert c.TH_RELAY>c.M_VG_MAX
assert c.CAP_H-c.FIXED_LEAVES-(c.TFULL["Sp"]*c.cf+1)-c.CAP_H//4>0
# Independent exact interval log: dyadic arithmetic, atanh power series.
S=1<<320

def cdiv(n,d):return (n+d-1)//d

def ln_range(value):
 assert 1<=value<=2
 t=(value-1)/(value+1); lo=t.numerator*S//t.denominator; hi=cdiv(t.numerator*S,t.denominator)
 lo2=lo*lo//S;hi2=cdiv(hi*hi,S)
 p0,p1=lo,hi;z0=z1=0
 for j in range(180):
  z0+=p0//(2*j+1);z1+=cdiv(p1,2*j+1)
  p0=p0*lo2//S;p1=cdiv(p1*hi2,S)
 # Remaining atanh terms after exponent359. t<=1/3 =>1/(1-t²)<=9/8.
 tail=cdiv(9*p1,8*361)
 return F(2*z0,S),F(2*(z1+tail),S)

l2,h2=ln_range(F(2))
x=c.E_TOTAL/1000
k=x.numerator.bit_length()-x.denominator.bit_length()
while x<F(2)**k:k-=1
while x>=F(2)**(k+1):k+=1
lo,hi=ln_range(x/F(2)**k)
loglo=F(k)+lo/h2;loghi=F(k)+hi/l2
assert loghi<=c.CLAIM and loglo>c.CLAIM-F(1,10**15)
assert c.CLAIM<F(123148512898914260,10**15)
# Gate is minimal for the stated success lower bound, no probabilistic inference here.
assert c.P_GATE*c.LO_RUN>=F(39,100)
assert F(c.K_GATE-1,1<<128)*c.LO_RUN<F(39,100)
assert c.E_TOTAL==256+F(c.K_GATE,1<<128)*sum(c.ROWS)
assert len(c.ROWS)==27 and all(row>=0 for row in c.ROWS)
# Card69: independent population arithmetic and measured n=1,2,3 census.
table_counts={'A':[5013,5758,6503],'B':[5013,5758,6503],'EB':[5013,5758,6503]}
for cl,counts in table_counts.items():
 assert [c.tbldk_count(n,cl) for n in (1,2,3)]==counts
 assert all(4268+n*(486+6*c.EISR3[cl])<count for n,count in enumerate(counts,1))
rows_per_node=2*((1<<52)+4*(1<<51))+3*(1<<45)+2*(1<<44)
assert rows_per_node==193*(1<<47)==2*c.ROWS_52+c.ROWS_F
table_extra=c.NODES*rows_per_node
assert c.ROWS[7]==c.NODES*c.XTBLD_ALL
tablecheck=dict(fixed_intercept=4268,row_base=487,measured_counts=table_counts,
 rows_per_node=str(rows_per_node),additional_shared_instructions=str(table_extra),
 old_486_price_rejected=True,scope='Independent table-population arithmetic and comparison with measured finite assembly counts.')
# Independent B byte census: disjoint populated intervals, plus shared code/stack reserve.
bank=(1<<62)+(1<<31)+3*(1<<24)+(1<<20)+(2<<20)
assert bank==c.B_BYTES
assert c.PEAK_BYTES==(c.NODES+1)*c.G_BYTES+(1<<40)*c.R_BYTES+(1<<82)*c.H_BYTES+(1<<58)*bank+c.FABRIC_BYTES
assert c.PEAK_BYTES<(1<<138)
A=0xd800000000000000
assert A%2048==0
assert ("REQ_CURSOR_MONITOR",A-2048,A) in c.G_MEM_REGIONS
assert any(lo<=A-8 and A<=hi for _,lo,hi in c.G_MEM_REGIONS)
assert not any(lo<=A-8<hi for name,lo,hi in c.G_MEM_REGIONS if name!="REQ_CURSOR_MONITOR")
assert c.PEAK_BYTES-2048*(c.NODES+1)-(3<<20)*c.N_REL-((2<<20)+65536)*c.N_HOME==336726598245340192209911994070362352402976
# Service addresses independently transcribed from emitted source and callback ABI.
G=34*(1<<40);R=192*(1<<40);H=1<<55;B=1<<62
expect={"G":G+0x7800000,"R":R+(20<<20),"H":H-(32<<20),"B":B+(1<<34)}
assert expect==c.GENERIC_BASES
for role,base in expect.items():
 regs=c.ROLE_MAPS[role]
 for offset,length in [(0,5*65536),(3*65536+8192,16*2048)]:
  assert any(lo<=base+offset and base+offset+length<=hi for _,lo,hi in regs)
 for i in range(16):
  a=base+3*65536+8192+2048*i;b=a+2048
  assert any(lo<=a<b<=hi for _,lo,hi in regs)
assert ("R_SLOT_TABLES",R+(32<<20),R+(32<<20)+1024*8192) in c.R_MEM_REGIONS
assert ("R_SLOTS_1024x64K",196*(1<<40),196*(1<<40)+1024*65536) in c.R_MEM_REGIONS
assert ("H1_CONTIG",64<<50,(64<<50)+(1<<56)) in c.G_MEM_REGIONS and ("H2U_CONTIG",8<<60,(8<<60)+(1<<57)) in c.R_MEM_REGIONS
assert ("H_BULK_CTRL",H-(64<<20),H-(64<<20)+65536) in c.H_MEM_REGIONS
assert any(lo<=expect["H"]+4*65536+24576 and expect["H"]+5*65536<=hi for _,lo,hi in c.H_MEM_REGIONS)
for role,base in c.PROGRAM_BASES.items():
 assert any(lo<=base and base+(1<<20)<=hi for _,lo,hi in c.ROLE_MAPS[role])
assert c.G_BYTES==sum(b-a for _,a,b in c.G_MEM_REGIONS)+(3<<20)
assert c.R_BYTES==sum(b-a for _,a,b in c.R_MEM_REGIONS)+(3<<20)
assert c.H_BYTES==sum(b-a for _,a,b in c.H_MEM_REGIONS)+(3<<20)
# Broken old maps fail containment, not just a region-name comparison.
for role,base,removed in [("R",expect["R"],"R_GENERIC"),("H",expect["H"],"H_GENERIC"),("R",R+(32<<20),"R_SLOT_TABLES"),("H",H-(64<<20),"H_BULK_CTRL")]:
 old=[v for v in c.ROLE_MAPS[role] if v[0]!=removed]
 assert not any(lo<=base and base+64<=hi for _,lo,hi in old)

for m in (c.G_MEM_REGIONS,c.R_MEM_REGIONS,c.H_MEM_REGIONS,c.B_MEM_REGIONS):
 rows=sorted(m,key=lambda v:v[1]);assert all(0<=a<b<=1<<64 for _,a,b in rows)
 assert all(b<=d for (_,a,b),(_,d,e) in zip(rows,rows[1:]))
# Source-measured service bounds: begin4*n+9, commandbegin4*n+62,
# three sweep phases64*nbanks+21 perowner. 4096 includes all constants/call glue.
for q in qq:
 nb=(q+(1<<27)-1)//(1<<27)
 anchored_bound=4*q+3*nb*(64+62+128)+63
 charged=3*4096*nb+8*q+4096
 assert charged>=anchored_bound
 assert c.bank_control(q)==(1<<48)*charged+256
assert c.COLD_SCHEDULE==22001*46*16*65536 and c.COLD_SCHEDULE<1<<41
assert c.COLD_EXPLICIT_BOUND<=c.COLD_PROTOCOL_ALLOWANCE
assert c.PROGRAM_BYTES==1<<20
# Value controls, not merely assertions on constants.
controls={'old_bank_capacity_fails':max(qq)>1024*(1<<19),
          'score_one_grid_step_low_fails':loglo>c.CLAIM-F(1,10**15),
          'gate_one_unit_low_fails':F(c.K_GATE-1,1<<128)*c.LO_RUN<F(39,100),
          'omit_B_changes_memory':c.PEAK_BYTES-(1<<58)*bank!=c.PEAK_BYTES}
assert all(controls.values())
for f,expect in [(F(-123451,10000),'-12.3451'),(F(-123451,10000)+F(1,10**6),'-12.3450')]:
 assert c.dec_up(f,4)==expect
# Exhaustive finite fibers for the independently specified pair-product map.
from collections import Counter
cases=0;largest=0;over6=None
for prime in (43,61,67,79,97,127):
 pts=[(x,y) for x in range(prime) for y in range(prime) if (y*y-x*x*x-7)%prime==0]
 counts=Counter(y*y%prime for x,y in pts)
 assert max(counts.values())<=6
 for a in counts:
  for b in counts:
   if a>=b: continue
   fibers=Counter()
   for squared_y,multiplicity in counts.items():
    product=(squared_y-a)*(squared_y-b)%prime
    if product: fibers[pow(product,-1,prime)]+=multiplicity
   assert max(fibers.values(),default=0)<=12
   assert counts[a]+counts[b]<=12
   cases+=1
   for z,multiplicity in fibers.items():
    largest=max(largest,multiplicity)
    if multiplicity>6 and over6 is None:over6=dict(prime=prime,A0=a,A1=b,inverse=z,fiber=multiplicity)
assert over6 is not None
assert c.D_PAIR_INVERSE==F(36*(c.C0+1)*c.CH,c.r-1)
assert c.D_PAIR_INVERSE_CHAIN==F(12*(c.C0+1)*(c.Q2//2),c.r-1)
assert c.PAIR_INVERSE_ALL==c.NODES*c.Q1GD*c.CHK*c.D_PAIR_INVERSE
assert c.PAIR_INVERSE_ALL>=c.NODES*c.Q1GD*c.D_PAIR_INVERSE_CHAIN
paircheck=dict(cases=cases,max_fiber=largest,invalid_six_fiber_control=over6,
 chain_budget=str(c.D_PAIR_INVERSE_CHAIN),nominal_chunk_budget=str(c.D_PAIR_INVERSE),global_budget=str(c.PAIR_INVERSE_ALL),
 explicit_batch_addition=True,scope='Independent finite fiber and exact-price check; finite fields do not prove production success.')
summary=c.summary()
print(json.dumps(dict(all_checks_pass=True,source_sha256=hashlib.sha256(pathlib.Path(c.__file__).read_bytes()).hexdigest(),
 independent_schedule_boundaries=len(c.MQ1),independent_log_interval_bits=320,
 definition_checks=['exact_schedule','atanh_log_interval','minimal_gate','27_row_sum','B_map_and_peak','bank_control_envelope','cold_schedule','request_cursor_monitor_backing','ancillary_probability_anchors','signed_cap_gaps','all_role_generic_and_program_backing','source_tables_and_bulk_control'],
 controls=controls,paircheck=paircheck,tablecheck=tablecheck,summary=summary)))
