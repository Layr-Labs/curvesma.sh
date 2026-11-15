import random;rn,en,dc=range,enumerate,dict;M64=2**64-1;M41,M44=2**41-1,2**44-1;V63=1<<63;T34=34<<40;SP=T34+0x7000000;CLB=35<<40;MS=[]
# One processor with bank-selected memory (proof.md 2, premise BM). Kernel calls below run on a64.Bank machines
# (mach): near memory plus the four windows of certificate.WINS; the va55 probe runs on bare near memory. Selects
# are executed instructions: the driver's own go through bsel (dsb, stp, dsb, ret); the insert's, ystore's and
# hzero's are inside those kernels. The driver between kernel calls is this Python code, not counted instructions:
# its writes are set-up stand-ins (X entries, planted fractions, producer frames, list header) and go through the
# same device, with two exceptions: case(cur=...) presets a bank's cell count by a direct write (last_cell and
# overflow_full_drop), and selchk's err() clears the pending state by m.dsb() after the fault it expects.
class PZ(dict):
 # near memory: a never-written word reads as nonzero garbage (proof.md 5.1); nv counts such loads
 def __missing__(z,a):z.nv=getattr(z,"nv",0)+1;return pzf(a)
def pzf(a):return(a*0x9E3779B97F4A7C15^0x5A5A5A5A5A5A5A5A)%2**64|1
def mach(E,Q):m=E.a64.Bank(PZ(),Q.WINS);MS.append(m);return m
def sel(E,rt,Q,m,w,k0,k1=0):return E.arun(rt[0],rt[1],"bsel",dc(x0=Q.SEL[w],x1=k0,x2=k1),m)
def gsel(E,rt,Q,m,g):return sel(E,rt,Q,m,"S",g%2**64,g>>64)
def cl0(m,cap):m[SP]=CLB;m[SP+8]=CLB+32*cap
def cld(m):n=(m[SP]-CLB)//32;o=[tuple(m[CLB+32*i+8*w]for w in rn(4))for i in rn(n)];m[SP]=CLB;return o
class HT:
 # model of the hash banks (proof.md 4): bank (d1, K2 mod 2^44), slot K3 mod 2^bits, a slot's records oldest first;
 # ins returns the instructions the insert executes beyond the empty-slot path (keys site: 9+4-11=2 more than hocc)
 def __init__(z,Q,bits,site=2):z.Q,z.bits,z.site,z.b,z.n,z.cand,z.lcap,z.used,z.ex,z.bd,z.st=Q,bits,site,{},{},[],1<<30,set(),0,0,dc(occ=0,walk=0,k3=0,listed=0,full=0,drop=0,deep=0)
 def key(z,d1,K2,K3):return(d1&M41,K2&M44),K3&((1<<z.bits)-1)
 def ins(z,d1,K2,K3,w2,w3):
  bank,slot=z.key(d1,K2,K3);ch=z.b.setdefault((bank,slot),[]);r=(K2,K3,w2,w3);z.used.add(bank)
  if not ch:ch.append(r);return 0
  h,m,st,Q,k0,v=5,None,z.st,z.Q,z.st["k3"],1;st["occ"]+=1;e=ch[0]
  if e[1]==K3:
   h+=3;st["k3"]+=1
   if(e[0]^K2)>>41==0:m=e;h+=4
  if m is None:
   h+=2
   if len(ch)>1:
    h+=2;st["walk"]+=1;st["deep"]=max(st["deep"],len(ch))
    for e in ch[:0:-1]:
     h+=6;v+=1
     if e[1]==K3:
      h+=3;st["k3"]+=1
      if(e[0]^K2)>>41==0:m=e;h+=1;break
     h+=2
    else:h+=2
  if m is not None:
   full=len(z.cand)>=z.lcap;h+=(3 if full else 6)+2;st["full" if full else"listed"]+=1
   if not full:z.cand+=[(m[2],m[3],w2,w3)]
  n=z.n.get(bank,0)+1
  if n>>z.bits:h+=4;st["drop"]+=1
  else:h+=12;z.n[bank]=n;ch.append(r)
  h+=2+z.site;z.ex+=h;z.bd+=Q.OCC_PAIR*v+Q.OCC_VER*(st["k3"]-k0)+Q.OCC_LS1*(m is not None);return h
def chain(m,Q,bank,slot,bits):
 # a slot's records as stored: primary, then overflow cells by the numbers in the stored K3 words
 b=m.bank("H",bank)[0];a=Q.H_BASE+32*slot;k=(1<<bits)-1;o=[]
 while True:
  o+=[(b[a],b[a+8]&~k|slot,b[a+16],b[a+24])];n=(b[a+8]^slot)&k
  if not n:return o
  a=Q.OV_BASE+32*n
def hzr(E,rt,Q,m,keys,done,cq,ht=None):
 # executes the zeroing pass (hzero: select, DC ZVA, cell count) on each 128-KiB region of the primary area that a
 # record's slot falls in; the other passes are the same loop (proof.md 5.1)
 for bank,slot in keys:
  k=(bank,slot>>12)
  if k not in done:assert ht is None or bank not in ht.used,"zeroing a bank that already holds records";done.add(k);cq+=[("hzero",E.arun(*rt,"hzero",dc(x9=bank[0],x10=bank[1],x11=1,x6=Q.SEL["H"],x7=Q.H_BASE+(slot>>12<<17),x8=Q.ZW,x5=Q.OV_BASE),m),Q.HZ_BANK+Q.ZERO_PASS+1)]
def htok(m,Q,ht):
 # table contents equal the model: each slot's chain and each bank's cell count
 return all(chain(m,Q,b,s,ht.bits)==[ch[0]]+ch[:0:-1]for(b,s),ch in ht.b.items())and all(m.bank("H",b)[0][Q.OV_BASE]==ht.n.get(b,0)for b,_ in ht.b)
def kref(E,mem,I,q2,W2,Q):
 # model of keys on the fractions in memory: per-leaf path extras and the records in emission order
 meter,A,R=0,[],[]
 for l in rn(q2):
  e=I+E.ENTRY*l;w3=mem[I+E.ENTRY*(l-l%6)+288]
  for s in(0,1):
   J=E.fracs(mem,e,s);J[2],J[3]=E.ws(J[2],E.get(mem,e+256,2)),E.ws(J[3],E.get(mem,e+272,2));dev=0
   for t in E.KO:f=E.record_fields(E.leaf_w(J,t));meter+=Q.JX[t][f[0]];dev|=bool(f[0]);R+=[(f[1],f[3],f[4],W2+l%6,w3)];A+=[f]
   meter+=Q.SJOIN if dev else 0
 return meter,A,R
def ez(m,e,w=0):
 for o in(256,264,272,280):m[e+o]=0
 m[e+288]=w
def plant(E,mem,sb,I,kn=0):
 for l,J in en(I):
  ez(mem,sb+E.ENTRY*l,kn+(l//6<<42))
  for s in(0,1):
   c=E.inverse_dft(J[3*s:3*s+3]);v=[f(E.t128(E.BT[m+k]*c[i]%E.p))for i,k,f in((0,0,int),(1,0,E.wneg),(2,-1,int))for m in(0,1)]
   for o,w in zip(E.FRAC_OFF[s],v):E.put(mem,sb+E.ENTRY*l+o,w,2)
def krun(E,Q,rt,lk,m,ht,done,g,sb,q2,W2,cq,lcap=None,zero=True):
 # one chunk's keys call for virtual generator g on entries already in W_S(g): zeroing passes for its records'
 # regions, list header, keys, then the drained candidate list against the model
 meter,fl,R=kref(E,m,sb,q2,W2,Q)
 if zero:hzr(E,rt,Q,m,[ht.key(*r[:3])for r in R],done,cq,ht)
 c0=len(ht.cand);ht.lcap=c0+(Q.LIST_CAP if lcap is None else lcap);cl0(m,ht.lcap-c0);o={};n=E.arun(lk[0],lk[1],"keys",dc(x2=W2+6,x3=sb+E.KEYS0,sp=SP),m,out_regs=o,ls64=([(sb,sb+E.ENTRY*q2)],[]));x=sum(ht.ins(*r)for r in R);cq+=[("keys",n,Q.KEY_FIX+Q.KEY_INC*q2+meter+x)];return fl,cld(m)==ht.cand[c0:]and o["sp"]==SP
def part_b(E,Q):
 # Sequential mini-run: four virtual generators' planted leaves through keys into the hash banks, one after the
 # other on one machine, against the model (contents, counts, listings); then the directed cases.
 rg=random.Random(E.SEED+32452843);rt=E.kload("bank");lk=E.kload("leaf");m=mach(E,Q);ht=HT(Q,Q.SLOT_BITS);cq,done,G,q2=[],set(),4,18;D1s,D2s=(0,1,1<<40,M41),(0,1,M41-1,M41);pool=[E.planted_x(rg,rg.choice(D1s),rg.choice(D2s))[1]for _ in rn(10)];batches,shared,ok,lab={},0,True,True
 for g in rn(G):
  ar=[]
  for _ in rn(6*q2):
   if rg.random()<0.15:ar+=[E.planted_x(rg,0,0,rg.choice(pool))[0]];shared+=1
   else:ar+=[E.planted_x(rg,rg.choice(D1s),rg.choice(D2s),fx=(rg.rr(2),rg.rr(4),None))[0]]
  batches[g]=[ar[6*l:6*l+6]for l in rn(q2)]
 gs=[rg.rr(3<<64,1<<66),rg.rr(1<<64),5,(1<<66)-1];sb=Q.EB_END-E.ENTRY*q2;KN=[V63|rg.gb(42)for _ in gs]
 for i,g in en(gs):cq+=[("bsel",gsel(E,rt,Q,m,g),Q.BSEL)];plant(E,m,sb,batches[i],KN[i])
 for i,g in en(gs):cq+=[("bsel",gsel(E,rt,Q,m,g),Q.BSEL)];hzr(E,rt,Q,m,[ht.key(*r[:3])for r in kref(E,m,sb,q2,8*i,Q)[2]],done,cq)
 for i,g in en(gs):
  cq+=[("bsel",gsel(E,rt,Q,m,g),Q.BSEL)];fl,o=krun(E,Q,rt,lk,m,ht,done,g,sb,q2,8*i,cq,zero=False);ok&=o;lab&=[f[1:]for f in fl]==[E.class_label(x[3*s+t])for x in batches[i]for s in(0,1)for t in E.KO]
 st=ht.st;tab=htok(m,Q,ht);ids=all(w3>>63 and(w3-KN[w2>>3])%2**42==0 and(w3-KN[w2>>3])>>42 in(0,1,2)and w2%8<6 for ch in ht.b.values()for _,_,w2,w3 in ch);cnt,av={},{}
 for e,c,f in cq:cnt[e],av[e]=cnt.get(e,0)+c,av.get(e,0)+f
 two=sum(len(ch)>1 and len({r[:2]for r in ch})>1 for ch in ht.b.values());A=dc(leaves=6*q2*G,slots=len(ht.b),occupied_inserts=st["occ"],two_labels_one_slot=two,listed=st["listed"],k3_equal=st["k3"],walks=st["walk"],longest_chain=max(map(len,ht.b.values())),counts=cnt,formulas_ok=all(c==f for _,c,f in cq)and ok and lab and tab and ids and two>0 and st["listed"]>0 and st["walk"]>0 and shared>0 and ht.ex<=ht.bd,ledger_bound=[ht.ex,ht.bd])
 A["occupied"]=occ(E,Q,rt,lk,rg);A["select"]=selchk(E,Q,rt);A["zero"]=zchk(E,Q,rt,rg);A["y_bank"]=ybank(E,Q,rt,rg);A["phase0"]=p0x(E,Q,rt,lk,rg);A["gate"]=gtc(E,rt,Q);A["ls64"]=l6c(E,Q,rg);A["va55"]=va55_checks(E);A["encoding"]=enc(E);A["keys_edges"]=keys_edges(E,Q,rt,lk);A["meter"]=meter_chk(E,Q,rt,lk,rg);A["registry_trace"]=reg_trace(Q);A["poison_loads"]=plc(E);return A
def occ(E,Q,rt,lk,rg):
 # Directed occupied-path cases, each a keys call on planted leaves (bank (0,0,d4lo), slot klo): a slot with three
 # different K3 (walk of two cells); equal K3 with different d4 (not listed); equal labels (listed once, early
 # exit); the list full (not listed, still stored); a match in an overflow cell behind a slot record of another K3,
 # then a K3-equal record with another d4 that walks past both; the overflow area at its last cell and full
 # (dropped); after_phase0: the first case after Phase-0 inserts (xemit) filled two cells of its bank and a zeroing pass.
 B=Q.SLOT_BITS;px=lambda fx,v=None:E.planted_x(rg,0,0,v,fx=fx);A={};res=True
 def pd(fx,a):
  # a leaf like px(fx) whose d4 differs from leaf a's
  while True:
   x=px(fx)
   if E.label(x[1])[2]!=E.label(a[1])[2]:return x
 def case(name,xs,lcap=None,cur=None,exp=None,p0=0):
  nonlocal res
  m=mach(E,Q);ht=HT(Q,B);cq,done=[],set();q2=6*-(-len(xs)//36);ls=xs+[E.planted_x(rg,2+k,5)[0]for k in rn(6*q2-len(xs))];sb=Q.EB_END-E.ENTRY*q2;gsel(E,rt,Q,m,9);plant(E,m,sb,[ls[6*l:6*l+6]for l in rn(q2)],V63)
  R=kref(E,m,sb,q2,0,Q)[2]
  if p0:
   k=ht.key(*R[0][:3]);hzr(E,rt,Q,m,[k],set(),[])
   for u in rn(3):E.put(m,Q.YR_BASE+512*u+288,u<<241|k[1]<<207|k[0][1]<<143|1)
   E.arun(*lk,"xemit",dc(x1=Q.YR_BASE,x3=3,x2=9,x26=V63,x6=0,x7=Q.SEL["H"],x4=Q.H_BASE,sp=SP),m);p0-=m.bank("H",k[0])[0][Q.OV_BASE]
  hzr(E,rt,Q,m,[ht.key(*r[:3])for r in R],done,cq)
  if cur is not None:b=ht.key(*R[0][:3])[0];m.bank("H",b)[0][Q.OV_BASE]=cur;ht.n[b]=cur
  fl,o=krun(E,Q,rt,lk,m,ht,done,9,sb,q2,0,cq,lcap=lcap,zero=False);st=dc(ht.st);st.pop("deep");ok=o and not p0 and htok(m,Q,ht)and ht.ex<=ht.bd and all(c==f for _,c,f in cq)and(exp is None or all(st[k]==v for k,v in exp.items()));A[name]=dc(st,ok=ok);res&=ok;return m,ht
 a,b,c=px((1,3,5))[0],px((1,3,6))[0],px((1,3,7))[0];case("chain_of_three",[a,b,c],exp=dc(occ=2,walk=1,k3=0,listed=0,drop=0))
 g=px((2,4,9));case("k3_equal_d4_differs",[g[0],pd((2,4,9),g)[0]],exp=dc(occ=1,k3=1,listed=0))
 f=px((3,5,2))[1];case("equal_labels_listed",[px(None,f)[0]for _ in rn(4)],exp=dc(occ=3,listed=3,k3=3,walk=0));case("list_full",[px(None,f)[0]for _ in rn(4)],lcap=1,exp=dc(occ=3,listed=1,full=2));g=px((5,7,2));case("match_in_overflow_cell",[g[0],px((5,7,1))[0],px(None,g[1])[0],pd((5,7,2),g)[0]],exp=dc(occ=3,walk=2,k3=3,listed=1))
 xs=[px((4,6,k))[0]for k in(1,2,3)];m,ht=case("last_cell",xs[:2],cur=2**B-2,exp=dc(occ=1,drop=0));top=m.bank("H",next(iter(ht.b))[0])[0];A["last_cell"]["cell_at_window_top"]=Q.OV_BASE+32*(2**B-1)==Q.H_BASE+2**53-32 and Q.OV_BASE+32*(2**B-1)in top
 m,ht=case("overflow_full_drop",xs,cur=2**B-1,exp=dc(occ=2,drop=2));A["overflow_full_drop"]["only_primary_stored"]=[len(ch)for ch in ht.b.values()if len(ch)]==[1]*len(ht.b);case("after_phase0",[a,b,c],p0=2)
 pc=[Q.OCC_FIX,Q.OCC_FIX+Q.OCC_WALK+Q.OCC_VISIT,Q.OCC_FIX+Q.OCC_VER,Q.OCC_FIX+Q.OCC_VER+Q.OCC_LS1,Q.OCC_FIX+Q.OCC_VER+Q.OCC_LS0,Q.OCC_FIX-Q.OCC_DROP];h=HT(Q,B);K=1<<60;r=[h.ins(0,0,5,0,V63),h.ins(0,0,5|K,0,V63),h.ins(0,0,5|2*K,0,V63)];h2=HT(Q,B);h2.ins(0,0,5,0,V63);r+=[h2.ins(0,1<<50,5,0,V63),h2.ins(0,0,5,0,V63)];h2.lcap=1;h2.cand=[0];r+=[h2.ins(0,0,5,0,V63)];h3=HT(Q,B);h3.ins(0,0,5,0,V63);h3.n[(0,0)]=2**B-1;r+=[h3.ins(0,0,5|K,0,V63)]
 A["path_counts"]=dc(empty=r[0],occupied=r[1],walk_one=r[2],k3_equal=r[3],listed=r[4],list_full=r[5],dropped=r[6]);A["ok"]=res and r==[0]+pc and pc[3]==36 and A["last_cell"]["cell_at_window_top"]and A["overflow_full_drop"]["only_primary_stored"]and all(x.ex<=x.bd for x in(h,h2));return A
def selchk(E,Q,rt):
 # Select discipline on the device of premise BM, each case by name: a load with no bank selected, in another bank
 # or of a never-written word is poison; a select store with a window access since the last DSB (first DSB missing)
 # is an error; a window load after a select store and before a DSB (second DSB missing) is poison, a store there
 # or a RET there is an error; bank numbers are decoded modulo the stated widths; DC ZVA zeroes; windows are
 # independent; a select is three instructions.
 a64=E.a64;P,L=a64.a64_load("rd:\nldr x2,[x1]\nret\nwr:\nstr x2,[x1]\nret\nnof:\nldr x3,[x4]\nstp x1,x2,[x0]\ndsb sy\nret\nnos:\ndsb sy\nstp x1,x2,[x0]\nldr x3,[x4]\ndsb sy\nldr x5,[x4]\nret\nnosw:\ndsb sy\nstp x1,x2,[x0]\nstr x3,[x4]\ndsb sy\nret\nnor:\ndsb sy\nstp x1,x2,[x0]\nret\nzv:\ndc zva,x1\nret\nrreg:\nldr x2,[x0]\nret\n");m=mach(E,Q);H,Y,SH=Q.H_BASE+64,Q.YQ_BASE+64,Q.SEL["H"];A={};pz=a64.pz
 def rd(a):o={};E.arun(P,L,"rd",dc(x1=a),m,out_regs=o);return o["x2"]
 def err(e,r,s):
  try:E.arun(P,L,e,r,m);return False
  except AssertionError as x:m.dsb();return s in str(x)
 A["no_bank_selected"]=rd(H)==pz(H)and err("wr",dc(x1=H,x2=7),"no bank selected");c=sel(E,rt,Q,m,"H",5,9);A["never_written_word"]=rd(H)==pz(H);E.arun(P,L,"wr",dc(x1=H,x2=77),m);A["written_word"]=rd(H)==77 and c==Q.BSEL==Q.SELECT+1
 sel(E,rt,Q,m,"H",6,9);A["other_bank"]=rd(H)==pz(H);sel(E,rt,Q,m,"H",5|1<<41,9|1<<44);A["decode_mod_widths"]=rd(H)==77 and m.key[m.names.index("H")]==(5,9)
 A["missing_first_dsb"]=err("nof",dc(x0=SH,x1=5,x2=9,x4=H),"since the last DSB");o={};E.arun(P,L,"nos",dc(x0=SH,x1=5,x2=9,x4=H),m,out_regs=o);A["missing_second_dsb"]=o["x3"]==pz(H)and o["x5"]==77 and err("nosw",dc(x0=SH,x1=5,x2=9,x3=1,x4=H),"no bank selected")and err("nor",dc(x0=SH,x1=5,x2=9),"without a DSB")
 sel(E,rt,Q,m,"Y",3);E.arun(P,L,"wr",dc(x1=Y,x2=88),m);A["windows_independent"]=rd(H)==77 and rd(Y)==88;b=Q.H_BASE+4096;E.arun(P,L,"wr",dc(x1=b+8,x2=5),m);E.arun(P,L,"zv",dc(x1=b),m);A["zva_zeroes"]=rd(b+8)==0 and rd(b+2040)==0 and rd(b+2048)==pz(b+2048)and rd(H)==77
 A["zva_select_register"]=err("zv",dc(x1=SH),"Device");A["load_select_register"]=err("rreg",dc(x0=SH),"bank-select register");A["ok"]=all(A.values());return A
def zchk(E,Q,rt,rg):
 # hzero on poisoned banks: two consecutive banks, two passes each; every word of the range reads zero afterwards,
 # the next block is still poison, the cell counts are zero, and the count is the charged one
 m=mach(E,Q);n=E.arun(*rt,"hzero",dc(x9=7,x10=M44-0,x11=2,x6=Q.SEL["H"],x7=Q.H_BASE,x8=2*Q.ZW,x5=Q.OV_BASE),m);bs=[m.bank("H",(7,M44)),m.bank("H",(7,0))];P,L=E.a64.a64_load("rd:\nldr x2,[x1]\nret\n");ok=n==2*(Q.HZ_BANK+2*Q.ZERO_PASS)+1 and all(len(b[1])==2*Q.ZW//256 and b[0]=={Q.OV_BASE:0}for b in bs)
 for k,w in((M44,0),(0,1)):
  sel(E,rt,Q,m,"H",7+w,k);v=[]
  for a in(Q.H_BASE,Q.H_BASE+16*Q.ZW-8,Q.H_BASE+16*Q.ZW,Q.OV_BASE):o={};E.arun(P,L,"rd",dc(x1=a),m,out_regs=o);v+=[o["x2"]]
  ok&=v==([0,0,E.a64.pz(Q.H_BASE+16*Q.ZW),0]if not w else[E.a64.pz(a)for a in(Q.H_BASE,Q.H_BASE+16*Q.ZW-8,Q.H_BASE+16*Q.ZW,Q.OV_BASE)])
 return dc(ok=ok,insns=n,banks=2)
def ybank(E,Q,rt,rg):
 # ystore: producer records (W_S of the producer, 448 B at frame bytes 64-511, slice word at byte 56) go to
 # W_Y(sigma) at the bank's cursor; set-up writes a bank's six header words; id words; clamp at the limit with the
 # sink frame; cursors; counts
 m=mach(E,Q);g,cls,KN=rg.rr(1<<66),rg.rr(1<<15),V63|rg.gb(41)<<1;A0=Q.YR_BASE;sl=[rg.rr(1<<41)for _ in"012"];Y=[];c=gsel(E,rt,Q,m,g)==Q.BSEL;lim=Q.YQ_BASE+64+512*2
 for u,s in en((0,1,0,2,0,0,1)):
  Y+=[(sl[s],[rg.gb(64)for _ in rn(56)])];m[A0+512*u+56]=sl[s]
  for w,v in en(Y[u][1]):m[A0+512*u+64+8*w]=v
 for s in sl:c&=sel(E,rt,Q,m,"Y",s<<15|cls)==Q.BSEL and E.arun(*rt,"tinit",dc(x0=6,x1=Q.YQ_BASE-64,x2=Q.YQ_BASE+64),m)==6*Q.TINIT_ENTRY+1
 gsel(E,rt,Q,m,g);o={};n=E.arun(*rt,"ystore",dc(x15=A0+64,x11=7,x26=cls,x24=KN,x25=1<<42,x14=Q.SEL["Y"],x13=Q.YQ_BASE-64,x12=lim),m,out_regs=o,ls64=([(A0,A0+512*7)],[(Q.YQ_BASE,lim+448)]));W=lambda u:Y[u][1][:37]+[KN+(u<<42)]+Y[u][1][38:];ok=n==7*Q.YSTORE+Q.YSTORE_FIX and o["x24"]==KN+(7<<42)and c
 for i,s in en(sl):
  us=[u for u in rn(7)if Y[u][0]==s];b=m.bank("Y",(s<<15|cls,0))[0];fr=lambda j:[b.get(Q.YQ_BASE+512*j+64+8*w)for w in rn(56)];ok&=[fr(j)for j in rn(min(2,len(us)))]==[W(u)for u in us[:2]]and b[Q.YQ_BASE-64]==min(Q.YQ_BASE+64+512*len(us),lim)and all(b[Q.YQ_BASE-64+8*w]==Q.YQ_BASE+64 for w in rn(1,6))and(fr(2)==W(us[-1])if len(us)>2 else Q.YQ_BASE+512*2+64 not in b)and all(Q.YQ_BASE+512*j+8*w not in b for j in rn(3)for w in rn(8))
 src=all(m.bank("S",(g%2**64,g>>64))[0][A0+512*u+64+8*w]==v for u in rn(7)for w,v in en(Y[u][1]));return dc(ok=ok and src,insns=n,records=7,sink_used=sum(y[0]==sl[0]for y in Y)-2,counts=[Q.YSTORE,Q.YSTORE_FIX])
def keys_edges(E,Q,rt,lk):
 rg=random.Random(E.SEED+104729);m,q2,M,cq=mach(E,Q),12,E.M64,[];sb=Q.EB_END-E.ENTRY*q2;gsel(E,rt,Q,m,3)
 for l in rn(q2):
  for s in(0,1):
   w=[rg.gb(128)for _ in rn(6)];a=rg.gb(62);f=[(a,a),(a,M-2*a),(a,-2*a&M),(a,-a&M),(M-2*a,a),(-2*a&M,a)][(2*l+s)%6]
   w[4]=(f[0]+(w[0]>>64)+(w[1]>>64)-(w[2]>>64)-(w[3]>>64))%2**64<<64|w[4]&M;w[5]=(f[1]-(w[0]>>64)+(w[2]>>64))%2**64<<64|w[5]&M
   for o,v in zip(E.FRAC_OFF[s],w):E.put(m,sb+E.ENTRY*l+o,v,2)
  ez(m,sb+E.ENTRY*l,V63)
 ht=HT(Q,Q.SLOT_BITS);fl,o=krun(E,Q,rt,lk,m,ht,set(),3,sb,q2,7,cq);j=[fl[6*l+3*s+1][0]for l in rn(q2)for s in(0,1)];return dc(ok=o and all(c==f for _,c,f in cq)and htok(m,Q,ht)and j[:6]==[1,0,2,0,1,2],leaf0_j=j[:6])
def meter_chk(E,Q,rt,lk,rg):
 # the driver's skip rule (proof.md 3.4) as a Python predicate, not executed instructions: a chunk's leaf kernels
 # run only while meter <= DELTA_F - FIX_BATCH; a skipped chunk calls no kernel, so the table keeps its 36 records
 run=lambda meter:meter<=Q.DELTA_F-Q.FIX_BATCH;m=mach(E,Q);ht=HT(Q,Q.SLOT_BITS);cq=[];sb=Q.EB_END-E.ENTRY*6;gsel(E,rt,Q,m,1);plant(E,m,sb,[[E.planted_x(rg,1,1)[0]for _ in rn(6)]for _ in rn(6)],V63);n=0
 for meter in(Q.DELTA_F-Q.FIX_BATCH+1,Q.DELTA_F-Q.FIX_BATCH):
  if run(meter):krun(E,Q,rt,lk,m,ht,set(),1,sb,6,0,cq);n+=1
 return dc(ok=n==1 and len(cq)>1 and sum(map(len,ht.b.values()))==36 and all(c==f for _,c,f in cq),chunks_run=n)
def p0x(E,Q,rt,lk,rg):
 # Phase 0 through the same insert (xemit, slot width SLOT_M): A values built to share a slot with different K3,
 # equal K3 with different d4, equal A[255:102] in one pool (listed) and in two pools (not listed); then psel
 m=mach(E,Q);B=Q.SLOT_M;ht=HT(Q,B,0);cq,done=[],set();YRB=Q.YR_BASE;g0,g1=rg.rr(1<<64,1<<66),rg.rr(1<<64);p0,p1=777,778;aa=rg.gb(153)<<102;lo=lambda:rg.gb(62)|1;bit=lambda k:1<<k;ok=True
 calls=[(g0,p0,0,[aa|lo(),(aa^bit(250))|lo(),(aa^bit(200))|lo(),(aa^bit(250)^bit(245))|lo()]+[rg.gb(255)for _ in rn(28)]),(g1,p0,0,[aa|lo(),aa|lo()]+[rg.gb(255)for _ in rn(30)]),(g1,p1,32,[aa|lo()])]
 rec=lambda A,pl,g,w3:((A>>102)&M41,(A>>143)&M64,pl<<49|A>>207,g%2**64,w3);W3=lambda g,i:V63|g>>64|i<<2;RS=[[rec(A,pl,g,W3(g,i0+u))for u,A in en(As)]for g,pl,i0,As in calls]
 for R in RS:hzr(E,rt,Q,m,[ht.key(*r[:3])for r in R],done,cq,ht)
 for(g,pl,i0,As),R in zip(calls,RS):
  cq+=[("bsel",gsel(E,rt,Q,m,g),Q.BSEL)]
  for u,A in en(As):E.put(m,YRB+512*(i0+u)+288,A)
  c0=len(ht.cand);ht.lcap=c0+len(As);cl0(m,len(As));n=E.arun(lk[0],lk[1],"xemit",dc(x1=YRB+512*i0,x3=len(As),x2=g%2**64,x26=W3(g,i0),x6=pl<<49,x7=Q.SEL["H"],x4=Q.H_BASE,sp=SP),m);cq+=[("xemit",n,Q.XEMIT_FIX+Q.XEMIT_HOT*len(As)+sum(ht.ins(*r)for r in R))];ok&=cld(m)==ht.cand[c0:]
 st=ht.st;A0=calls[0][3][0];ok&=htok(m,Q,ht)and ht.ex<=ht.bd and all(c==f for _,c,f in cq)and(st["occ"],st["k3"],st["listed"],st["walk"],st["deep"])==(6,3,2,3,6)and ht.cand==[(g0%2**64,W3(g0,0),g1%2**64,W3(g1,u))for u in(0,1)]and len(chain(m,Q,*ht.key(*rec(A0,p0,0,0)[:3]),B))==7
 au,ag,J,B0=T34+0x20000,T34+0x40000,(0x987654321<<192)+0x77700000000,(0x987654321<<192);nd=mach(E,Q)
 def wrep(mm,i,a0,a1,c0,c1):
  o=au+128*i;E.put(mm,o,a0);E.put(mm,o+32,a1);E.put(mm,o+64,c0,2);E.put(mm,o+80,c1,2)
  for w in rn(4):mm[o+96+8*w]=64*i+w+1
 def rps(p_,l_,n,mm):o={};s=E.arun(p_,l_,"psel",dc(x0=au,x1=n,x2=ag,sp=SP),mm,out_regs=o);return o["x0"],s
 ac=1<<80
 for i,rv in en([(0,9,1,2),(J+50,J+10,1,2),(J+50,J+50,1,2),(J+50,J+50+(1<<63),1,2),(J+50,J+50+2**64,1,2),(J+50,J+90,ac|7,ac|7),(J+100,J+200,10,11),(J+300,J+400,11,12),(J+100,J+500,12,13),(J+600,J+200,12,13),(J+100+(1<<192),J+200+(1<<192),ac|10,ac|11),(J+200+(1<<192),J+250+(1<<192),14,15),(J+(2<<192),J+(2<<192)+((1<<63)-1),14,15),(B0+(3<<192)+(1<<63)-5,B0+(3<<192)+(1<<63)+5,16,17),(J+(4<<192),J+(4<<192)+9,18,19)]):wrep(nd,i,*rv)
 rb,rq=rt;az,av=rps(rb,rq,15,nd);ax,_=rps(rb,rq,13,nd);z0,_=rps(rb,rq,0,nd);sa=[E.get(nd,ag+128*q+32*s)for q in rn(4)for s in(0,1)];sc=[E.get(nd,ag+128*q+64+16*s,2)for q in rn(4)for s in(0,1)];ad=mach(E,Q)
 for q in rn(3):wrep(ad,q,B0+((q+1)<<192)+10,B0+((q+1)<<192)+20,2*q+1,2*q+2)
 for i in rn(3,1030):wrep(ad,i,B0+(3<<192)+15,B0+(3<<192)+20,999+2*i,1000+2*i)
 aw,al=rps(rb,rq,1030,ad);mp,ml=E.a64.a64_load(E.ksrc("bank").replace("tbnz x12,#63,5f\n",""));gm=mach(E,Q)
 for q in rn(4):wrep(gm,q,J+(q<<192)+10,J+(q<<192)+10+(1<<63),2*q+1,2*q+2)
 ok&=(az==1 and ax==0 and z0==0 and aw==0 and rps(mp,ml,4,gm)[0]==1 and rps(rb,rq,4,gm)[0]==0 and al<=12+220*1024<=Q.PSEL_POOL and((sa[4]>>63)!=(sa[5]>>63))and((sa[6]>>63)!=(sa[7]>>63))and len(set(sa))==8 and sc==[10,11,ac|10,ac|11,14,15,16,17]and all(sa[2*q]>0 and 0<sa[2*q+1]-sa[2*q]<(1<<63)for q in rn(4)));return dc(ok=ok,points=sum(len(c[3])for c in calls),listed=st["listed"],k3_equal=st["k3"],occupied=st["occ"],xemit_32_insns=cq[-3][1],psel_ok_steps=av,psel_wc_steps=al)
def gtc(E,rt,Q):
 rb,rq=rt;m=mach(E,Q);K=Q.K_GATE;a5,J=K&(2**64-1),K>>64;cases=[(0,True),(K-1,True),(K,False),(K+1,False),(((J-1)<<64)|(2**64-1),True),((J<<64)|(2**64-1),False),((1<<128)-1,False),];ok=True;I=set()
 for u,a4 in cases:o={};n=E.arun(rb,rq,"gate",dc(x0=u&(2**64-1),x1=u>>64,x2=a5,x3=J),m,out_regs=o);I.add(n);ok&=(o["x0"]==(1 if a4 else 0))
 ok&=(I=={4}and 2*Q.COIN_SCALAR+4+8<=Q.GATE_INSNS);return dc(ok=ok,counts=sorted(I),cases=len(cases))
def enc(E):
 a64=E.a64;good="ldp x1,x2,[x3,#-512]|ldp x1,x2,[x3,#504]|stp x1,x2,[x3],#504|stp x1,x2,[x3,#-512]!|ldr x1,[x3,#32760]|ldr x1,[x3,#-256]|ldr x1,[x3,#255]|str x1,[x3],#-256|add x1,x1,#4095|add x1,x1,#0x1000|tst x3,#0x3fffe0000|mov sp,x1|add sp,sp,#512|cset w9,hs|csel x1,x2,x3,lo|mov x1,#0xfff8000000000000|movk x4,#1,lsl #32|ldr x22,[x3,#-152]|and x1,x2,#0x5555555555555555|cbhi x1,x2,f|cbhs xzr,x30,f|mvn x0,x0,lsr #23|stp x16,xzr,[x14]|tbnz x0,#47,f|tst x5,#0xffff800000000000|orr x6,x4,#0xfffffffffffe0|bfi x4,x30,#5,#47|mov x4,#0x60000000000000|mov x1,#0x220400000000|dsb sy|dc zva,x1".split("|");bad="stp x25,x26,[x3],#768|ldp x21,x22,[x3,#-768]|ldp x1,x2,[x3,#508]|ldp x1,x2,[x3,#-520]|ldp x1,x2,[x3,#12]|ldr x22,[x3,#-664]|ldr x1,[x3,#32768]|ldr x1,[x3,#260]|str x1,[x3],#256|ldr x1,[x3,#-257]!|add x1,x1,#4097|tst x3,#0x1c001|and x1,x2,#0|mov x1,#0x12345|movk x4,#1,lsl #8|ldp x1,x1,[x3]|ldp x3,x2,[x3],#16|stadd x4,[sp]|dmb ish|dsb ish|ldaxr x1,[x0]|ldar x1,[x0]|stlr x1,[x0]|ldadd x1,x2,[x3]|swp x1,x2,[x3]|blr x1|wfe|ubfx x1,x2,#60,#8|bfi x1,x2,#40,#25|ldr x1,[x3,#-264]|add x1,xzr,#5|tst sp,#0x1c000|ldp xzr,xzr,[x3]|ldr sp,[x3]|mov xzr,#0x3ffe0|movk sp,#1,lsl #48|tbz x3,#64,1f|frob x1,x2|mul x1,x2,#1|adc x1,x2,#3|csel x1,x2,#3,eq|madd x1,x2,x3,#4|cset x1,x2|mvn x1,#2|mul w1,x2,x3|dmb sy|umulh w1,w2,w3|mvn x1,x2,lsl #64|movk x1,#1,lsr #16|cbhi sp,x2,f|cbhs x1,#2,f|cbhi x1,x2|ld64b x0,[x2,#0]|ldp x1,x2,[x3]!|ld64b x21,[x1]|ld64b x24,[x1]|add x1,x2,#1_0|mov x1,#0x220400010000|dc zva,sp".split("|")
 def verdict(o):
  try:a64.a64_load("f:\n\t"+o+"\n\tret\n");return True
  except AssertionError as ex:assert"not encodable"in str(ex),ex;return False
 def reach(n):
  try:a64.a64_load("f:\n\tcbhi x1,x2,g\n"+"\tadd x1,x1,#1\n"*n+"g:\n\tret\n");return 1
  except AssertionError as ex:assert"branch reach"in str(ex);return 0
 A=dc(good_accepted=sum(map(verdict,good)),bad_rejected=sum(not verdict(i)for i in bad),good=len(good),bad=len(bad),cb_reach_1020_1024=(reach(254),reach(255)));A["ok"]=A["good_accepted"]==len(good)and A["bad_rejected"]==len(bad)and A["cb_reach_1020_1024"]==(1,0);return A
def va55_checks(E):
 a64=E.a64;I,J=a64.a64_load("ld:\n\tldr x2,[x1]\n\tret\n");A={};U=2**64-2**55
 for a in(0,2**55-8,U,2**64-8,2**55-4,2**55,2**56-8,2**60,U-8,2**64-2**56,2**63,2**64-4):
  try:E.arun(I,J,"ld",dc(x1=a),PZ());A[hex(a)]="ran"
  except AssertionError as ex:A[hex(a)]="rejected"if"55-bit"in str(ex)else"other: "+str(ex)
 v=list(A.values());A["ok"]=v==["ran"]*4+["rejected"]*8;return A
def l6c(E,Q,rg):
 # LS64: 64-byte alignment, register rules, the permitted regions of a call, and no LS64 on a select register
 a64=E.a64;I,J=a64.a64_load("ld:\n\tld64b x20,[x1]\n\tret\nst:\n\tst64b x16,[x2]\n\tret\novl:\n\tld64b x12,[x15]\n\tret\n");m=mach(E,Q);sel(E,E.kload("bank"),Q,m,"Y",1);Y=Q.YQ_BASE;R=([(Y,Y+4096),(Q.SEL["Y"],Q.SEL["Y"]+64)],[(Y,Y+4096),(Q.SEL["Y"],Q.SEL["Y"]+64)]);v=[rg.gb(64)for _ in rn(8)]
 def run(e,regs,lsr=R):
  try:E.arun(I,J,e,dc(regs),m,ls64=lsr);return"ran"
  except AssertionError as ex:return"LS64 rejected"if"LS64"in str(ex)else"other: "+str(ex)
 A={};A["st64b_window"]=run("st",dc(x2=Y+64,**{"x%d"%(16+k):v[k]for k in rn(8)}));o={};E.arun(I,J,"ld",dc(x1=Y+64),m,out_regs=o,ls64=R);A["ld64b_window_reads_back"]=[o["x%d"%(20+k)]for k in rn(8)]==v;A["ld64b_odd_slot"]=run("ld",dc(x1=Y+32));A["ld_out"]=run("ld",dc(x1=Y+8192));A["st_mis"]=run("st",dc(x2=Y+96));A["st_select_register"]=run("st",dc(x2=Q.SEL["Y"]));A["ld_base_ov"]=run("ovl",dc(x15=Y+64))
 def load(text):
  try:a64.a64_load(text);return"loaded"
  except AssertionError as ex:return"LS64 rejected"if"LS64"in str(ex)else"other: "+str(ex)
 A["ld64b_x24"]=load("hi:\n\tld64b x24,[x1]\n\tret\n");A["ld64b_x21"]=load("odd:\n\tld64b x21,[x1]\n\tret\n");rj=("ld64b_odd_slot","ld_out","st_mis","st_select_register","ld_base_ov","ld64b_x24","ld64b_x21");return dc(ok=A["st64b_window"]=="ran"and A["ld64b_window_reads_back"]and all(A[k]=="LS64 rejected"for k in rj),rejected=len(rj))
def plc(E):
 # loads that returned poison (a window word that no store or DC ZVA has written, or no bank selected) or
 # never-written near-memory garbage, by the entry that executed them, over every kernel call of this run: only
 # the negative probes have any (rd and nos of select and zero, ld of va55); no kernel has one
 return dc(by_entry=dc(E.PL),ok=E.PL==dc(rd=9,nos=1,ld=4))
def reg_trace(Q):
 # every address any machine of this run wrote lies in a registry interval of its kind: near memory in N, a bank's
 # words and zeroed blocks in its window's intervals; select registers hold no memory
 iv={k:[(lo,hi)for r,_,lo,hi,_ in Q.REGISTRY if r==k]for k in Q.KINDS};n=bad=0;inside=lambda k,a:any(lo<=a<hi for lo,hi in iv[k])
 for m in MS:
  for a in m.near:n+=1;bad+=not inside("N",a)
  for(w,key),b in m.banks.items():
   k=m.names[w];bad+=not(0<=key[0]<=m.mask[w][0]and 0<=key[1]<=m.mask[w][1])
   for a in b[0]:n+=1;bad+=not inside(k,a)
   for z in b[1]:n+=1;bad+=not(inside(k,z<<11)and inside(k,(z<<11)+2040))
 return dc(machines=len(MS),banks=sum(len(m.banks)for m in MS),addresses=n,outside=bad,ok=bad==0 and n>0)
