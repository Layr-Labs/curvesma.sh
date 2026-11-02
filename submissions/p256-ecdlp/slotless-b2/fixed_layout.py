def check(ledger):
 q=ledger['q1'];u=1<<28;side=1<<35;cap3=ledger['cap3'];cand=1<<24;req=1<<29
 ranges=[];cursor=0
 def alloc(name,size,alignment=64):
  nonlocal cursor
  base=(cursor+alignment-1)//alignment*alignment;ranges.append((name,base,base+size));cursor=base+size;return base
 alloc('code_stacks',1<<20,4096);output=alloc('output_mailboxes',16384,4096);control=alloc('contexts_ports_barriers',1<<16,4096)
 sub_ranges=[(0,5568),(6144,14336),(14336,22528),(22528,23296),(24576,40960)]
 assert all(a[1]<=b[0] for a,b in zip(sub_ranges,sub_ranges[1:])) and sub_ranges[-1][1]<=1<<16
 assert (control+24576)%2048==0
 alloc('source_X',64*q);alloc('source_Y',64*u);alloc('a_labels',32*q);alloc('row_b_labels',32*q);alloc('c_labels',32*u);alloc('global_b_Rg_rng',4096)
 batch=alloc('prefix_and_pair_batch',64*u+64,64)+64
 alloc('range_tree_records',72*u*29);alloc('range_tree_nodes',64*(2*u-1));alloc('sort_buffers',144*u)
 alloc('candidate_pairs',16*cand);alloc('request_output',128*cand);alloc('source_requests',64*req);alloc('relay',64*(1<<34));alloc('home_replies_flags',16*cand*64)
 alloc('home_port_window_state',2048,64);alloc('d3_cursor_end',16*(1<<26))
 assert cursor<1<<44
 ranges.append(('d4_counters',1<<51,1<<52));cursor=1<<52
 alloc('full_record_buckets',32*cap3*(1<<26));alloc('pointer_array',8*cap3*(1<<26))
 assert all(a[2]<=b[1] for a,b in zip(ranges,ranges[1:]));assert cursor<1<<64
 assert output==(1<<20) and output%4096==0 and batch%64==0
 assert ranges[-1][2]==ledger['memory_end']
 return dict(nonoverlapping=True,fixed_end_before_counters=ranges[-4][2],per_node_end=cursor,total_peak_upper_exponent=134,output_base=output,batch_base=batch,ranges=ranges,scope='explicit simultaneous allocation; all contexts must use these bases in final construction')
