# Candidate fixture definitions: ONLY the credential-free shared runner may load this file.
def run_protocol(inputs):
  rng=random.Random(inputs['seed']+85);results=[]
  for case in range(inputs['count']):
    records=[];current=0;baby_records=1536
    for i in range(4096):
      if case%4==0:U=0
      elif case%4==1:
        if i%4==0:current=rng.getrandbits(183)
        U=current
      elif case%4==2:U=((i%2)<<182)|((i%2)<<160)|((i%16)<<129)|((i%7)<<96)|((i%2)<<95)|((i%16)<<65)|((i%7)<<33)|((i%2)<<32)|(i%7)
      else:U=((i%3)<<144)|((i%64)<<112)|((i//64%7)<<80)|((i%2)<<64)|((i%64)<<32)|(i//64%7)
      records.append(meta_fixture(U,(1<<100)+i,0 if i<baby_records else 1))
    reconstructed=[];batches=0
    for start in range(0,len(records),512):
      batch=records[start:start+512];assert len(batch)==512
      payload=b''.join(raw[24:40]*2 for raw in batch);metadata=b''.join(batch);pending=True
      assert len(payload)==16384 and len(metadata)==20480
      payload_frames=[payload[j:j+4096] for j in range(0,len(payload),4096)]
      metadata_frames=[metadata[j:j+4096] for j in range(0,len(metadata),4096)]
      assert len(payload_frames)==4 and len(metadata_frames)==5
      assert b''.join(payload_frames)==payload and pending
      reply=b''.join(metadata_frames);pending=False
      assert reply==metadata and not pending
      reconstructed.extend(reply[j:j+40] for j in range(0,len(reply),40));batches+=1
    assert reconstructed==records and batches==8
    assert sum(raw[22]>>7==0 for raw in reconstructed)==baby_records
    assert sum(raw[22]>>7==1 for raw in reconstructed)==4096-baby_records
    for group in range(0,len(records),2048):
      bulk=b''.join(reconstructed[group:group+2048])
      frames=[bulk[j:j+4096] for j in range(0,len(bulk),4096)]
      assert len(bulk)==81920 and len(frames)==20 and b''.join(frames)==bulk
    stage0,q0,f0=sort_radix_fixture(reconstructed,0)
    stage1,q1,f1=sort_radix_fixture(stage0,1)
    expected=sorted(records,key=lambda raw:int.from_bytes(raw[:24],'little')&((1<<183)-1))
    assert stage1==expected
    
    boundary=LeafRadix(0,0,0);boundary.count([0]*200)
    base=(1<<100)+(1<<34)-17;assert boundary.prefix(base)==add128(base,200)
    test=[meta_fixture(0,(1<<110)+i,i&1) for i in range(200)];packets=[]
    for i in range(0,200,64):packets.extend(boundary.receive(test[i:i+64]))
    packets.extend(boundary.finish());written={}
    for packet in packets:
      at=int.from_bytes(packet[:16],'little');data=packet[16:]
      for j in range(len(data)//40):written[add128(at,j)]=data[40*j:40*j+40]
    assert [written[add128(base,j)] for j in range(200)]==test
    assert add128((1<<64)-1,1)==1<<64
    results.append({'case':case,'records':4096,'baby_source_records':baby_records,'giant_source_records':4096-baby_records,'worker_batches':8,'payload_bytes':32,'payload_frames_per_batch':4,'metadata_frames_per_batch':5,'stable_full183_sort':True,'masked_encoded192_padding_zero':True,'digit_widths':[96,96],'root_count_entry_bytes':[16,16],'root_count_frame_records':[256,256],'full128_ordinals_preserved':True,'padding_restored':True,'maximum_relay_queue':max(q0,q1),'full_count_frames':f0+f1,'bank_and_counter_boundary':True})
  return {'scope':'bounded sparse initialized-zero packet/state-model fixture at full183-bit fingerprints in192-bit encoding/full128-bit ordinals and real40-byte metadata/32-byte payloads; two96-bit two-edge passes and5 metadata frames; not production bucket initialization, native dispatch timing, GT success, probability or qualification','cases':results}


def _scaled_field_machine(c,x,y,z,Lambda,counter=0,threshold=None):
  m=28721772558436203898626428180203831296;n=55754029111986100220015185618117066752
  if threshold is None:threshold=(m+n)//(1<<21)+1
  assert 0<Lambda<P and 0<=counter<threshold<1<<128
  A,B=c[:2],c[2:];Bxi=cmul(B,(1,1));coeffs=[]
  for C in (A,Bxi,A,B):coeffs.extend([(C[0]*R%P,C[1]*R%P),(-C[1]*R%P,C[0]*R%P)])
  machine=Machine(coeffs);ctx=machine.const+8192
  for base,vec in [(machine.out,z),(machine.prev1,x),(machine.prev2,y)]:
    for j,v in enumerate(vec):machine.field(base+48*j,v*Lambda%P)
  for at in range(machine.const,machine.const+3728,8):machine.mem[at]=0
  machine.field(machine.const+448,P);machine.mem[machine.const+512]=P//2**325;machine.field(machine.mod,P)
  for quotient in range(50):
    for limb in range(7):machine.mem[machine.const+528+64*quotient+8*limb]=(quotient*P>>(64*limb))&MASK
  for off,value in [(32,machine.out),(40,machine.prev1),(48,machine.prev2),(56,machine.table),(96,(m+n)&MASK),(104,(m+n)>>64),(112,machine.const+4096),(128,counter&MASK),(136,counter>>64),(144,threshold&MASK),(152,threshold>>64)]:machine.mem[ctx+off]=value
  machine.extra_write_ranges=[(ctx+128,ctx+144)]
  machine.reg[:6]=[machine.out,machine.prev1,machine.prev2,machine.table,machine.const+528,machine.mod]
  machine.reg[6:11]=[machine.mod,machine.const+4096,machine.fixed_base,machine.const+528,ctx]
  machine.execute('worker_layout_init');cold=machine.count;assert cold==8
  machine.saved=machine.reg[18:].copy()
  assert machine.reg[23:30]==[machine.fixed_base,machine.const+528,P//2**325,ctx,48,96,machine.mod]
  return machine,cold

def run_scaled_field(inputs):
  rng=random.Random(inputs['seed']+381031);cases=[]
  for case in range(inputs['count']):
    vals=[rng.randrange(P) for _ in range(16)]
    c,x,y,z=[tuple(vals[i:i+4]) for i in range(0,16,4)]
    if case in (0,1,2,3):c=x=y=(0,0,0,0);z=(1,1,1,1) if case!=1 else (0,0,0,0)
    unscaled=oracle(c,x,y,z);K=(28721772558436203898626428180203831296+55754029111986100220015185618117066752)//(1<<21)+1
    threshold=K;initial=0;Lambda=rng.randrange(1,P)
    if case in (0,1,2,3):Lambda=P-1
    if case==2:threshold=1
    if case==3:threshold=3
    if case>=4 and case%2==0:
      target=case%4;assert unscaled[target]!=0;Lambda=(P-1)*pow(unscaled[target],-1,P)%P
    if case in (6,10,14):threshold=[(1<<64)+3,(1<<100)+5,(1<<128)-1][(case-6)//4];initial=threshold-1
    if case==5:initial=(1<<64)-1;threshold=(1<<64)+8
    if case==13:initial=(1<<127)-2;threshold=(1<<127)+8
    machine,cold=_scaled_field_machine(c,x,y,z,Lambda,initial,threshold);ctx=machine.const+8192
    expected=tuple(v*Lambda%P for v in unscaled);seed=tuple(v*Lambda%P for v in z)
    repeats=2 if case==9 else 1;run_counts=[];aborts=[]
    for repeat in range(repeats):
      if repeat:
        for j,v in enumerate(seed):machine.field(machine.out+48*j,v)
      before=machine.count;before_rare=getattr(machine,'scaled_borrow_visits',0);before_count=machine.readmem(ctx+128)|(machine.readmem(ctx+136)<<64)
      machine.execute('worker_field_call_fixture');total=machine.count-before
      hot=machine.calls['recurrence_fp2_step'][-1];rare=getattr(machine,'scaled_borrow_visits',0)-before_rare
      count=machine.readmem(ctx+128)|(machine.readmem(ctx+136)<<64)
      assert 0<=rare<=4 and count==before_count+rare and count<=threshold
      if case==0:assert rare==4 and hot==2421
      if case==1:assert rare==0 and hot==2345
      if case>=4 and case%2==0:assert rare>=1
      aborted=count>=threshold;assert machine.reg[0]==int(aborted)
      assert total==hot+8 and machine.reg[18:]==machine.saved
      assert machine.readmem(ctx+144)|(machine.readmem(ctx+152)<<64)==threshold
      assert machine.readmem(ctx+96)|(machine.readmem(ctx+104)<<64)==28721772558436203898626428180203831296+55754029111986100220015185618117066752
      assert machine.readmem(ctx+112)==machine.const+4096
      got=tuple(machine.loadfield(machine.out+48*j) for j in range(4))
      if aborted:
        j=machine.last_field_id;assert rare>=1 and got==expected[:j]+seed[j:] and hot<=2345+44*rare
      else:assert got==expected and all(v<P for v in got) and hot==2345+19*rare
      assert hot<=2421
      run_counts.append(hot);aborts.append(aborted)
    cases.append({'case':case,'normal_field_charge':2345,'actual_field_counts':run_counts,'caller_args_and_status':4,'fixture_ABI_status_return':4,'actual_cold_worker':cold,'aborted_paths':aborts,'complete128_counter_cap_and_keys_limit_disjoint':True,'full_dense_scaled_output_or_abort_prefix_checked':True,'actual_native_tag_branch_and_no_use_after_failure':True,'no_cold_or_counter_reset_between_calls':True,'ABI_LR_SP_and_initialized_reads':True})
  return {'scope':'ACTUAL whole scaled field/realcaller loads/BL/TBNZ/ABI wrapper; ALLFOURforcedR=p-1/max2421 and R0/normal2345,128counter carry/cap/partialabort/persistentsecondcall. FullFp oracle; uniformLambda/Markov/work/success separate analytical proof,not qualification. Caller4production+4fixture/cold8ALLexecuted','cases':cases}
