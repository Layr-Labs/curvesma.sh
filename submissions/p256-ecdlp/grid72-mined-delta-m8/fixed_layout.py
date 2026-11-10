def check(ledger):
 """Simultaneous per-node allocation (proof.md section 7): every region inside [0, 2^64), nonoverlapping, aligned."""
 mx=ledger['mx'];yb=ledger['yb'];side=1<<ledger['side_bits'];capg=ledger['cap_pg'];caph=ledger['cap_h'];cand=ledger['candidate_cap'];req=ledger['request_cap'];hb=ledger['hash_bits'];G=ledger['groups'];mp=mx//2
 ranges=[];cursor=0
 def alloc(name,size,alignment=64):
  nonlocal cursor
  base=(cursor+alignment-1)//alignment*alignment;ranges.append((name,base,base+size));cursor=base+size;return base
 fold=alloc('fold_table_and_constants',4096,4096)# F[c]=c(2^256-p) at 0 (F[1] is the Montgomery one); yprep constants at 512
 code=alloc('code_stacks',(1<<20)-4096,4096);output=alloc('output_mailboxes',16384,4096);control=alloc('contexts_ports_barriers',1<<16,4096)
 sub_ranges=[(0,5568),(6144,14336),(14336,22528),(22528,23296),(24576,40960),(40960,45056)]
 assert all(a[1]<=b[0]for a,b in zip(sub_ranges,sub_ranges[1:]))and sub_ranges[-1][1]<=1<<16
 assert(control+24576)%2048==0 and(control+40960)%4096==0
 alloc('x_points_scalars',96*mx);xa=alloc('x_areas',256*mx,512)# X area of pair J at xa+512J
 alloc('chain_coins',64*G*yb);cs=alloc('chain_state',128*yb,128);sc=alloc('chain_scratch',64*yb+64,64)
 ye=alloc('y_entries',1024*yb+1024,1024)+1024;batch=alloc('batch_entries',64*yb*mp+64,64)
 alloc('global_b_Rg_rng',4096);alloc('inverse_scratch',4096,64)
 alloc('candidate_pairs',16*cand);alloc('request_output',128*cand);alloc('source_requests',64*req);alloc('relay',64*(1<<34));alloc('home_replies_flags',16*cand*64)
 port=alloc('hop1_port_header',64);cursors=alloc('hop1_port_cursors',8*side,1<<39);assert cursors%(1<<39)==0;alloc('hop1_final_cursors',32*side)
 windows=alloc('fixed_block_payload_windows',8192*side+32,8192);assert windows%8192==0;alloc('block_remote_cursors',8*side);alloc('block_endpoint_metadata',64*side);alloc('block_delivery_context',4096,2048)
 assert cursor<1<<51
 assert(1<<51)+(1<<56)+32+32*capg*(side-1)+32*2*mx*yb*ledger['L']<1<<64# attempted remote hop1 cursor (BLOCK_PROTOCOL.md)
 tb=alloc('pair_tables',1<<56,1<<51)# pair pair P at 2^51+2^48P: interleaved sigma slots0,1; standard slots2..7 (2^45 bytes each)
 assert tb==1<<51 and(4*mp)<<45==1<<56 and cursor==(1<<51)+(1<<56)and tb%(1<<48)==0
 slots=alloc('hop1_group_slots',32*capg*side+32)+32;records=alloc('home_record_array',32*caph,8)
 assert cursor<=7<<61
 alloc('record_to_head_padding',(7<<61)-cursor,8)
 heads=alloc('hash_heads',8<<hb,1<<61)
 assert heads==7<<61 and slots%64==32
 assert all(a[2]<=b[1]for a,b in zip(ranges,ranges[1:]))and cursor==1<<64==ledger['memory_end']
 assert records%8==0 and records+32*caph<=heads and fold==0 and output==1<<20
 assert xa%512==0 and ye%1024==0 and batch%64==0 and cs%128==0
 ls64=dict(ld64b=[('pair_tables',tb,tb+(1<<56)),('batch_entries',batch,batch+64*yb*mp),('y_entries',ye,ye+1024*yb)],st64b=[('batch_entries',batch,batch+64*yb*mp)])
 return dict(nonoverlapping=True,fold_base=fold,x_area_base=xa,y_entry_base=ye,batch_base=batch,table_base=tb,chain_state_base=cs,output_base=output,port_base=port,cursor_table_base=cursors,window_base=windows,record_base=records,heads_base=heads,slots_base=slots,node_end=cursor,peak_exp=136,leaf_ls64=ls64,ranges=ranges,scope='explicit simultaneous allocation')
