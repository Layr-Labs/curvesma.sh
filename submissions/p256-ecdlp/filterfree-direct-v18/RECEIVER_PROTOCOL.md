# Two-phase direct-bucket receiver machine

This specifies the changed transport component. It is an abstract scalar AArch64 instruction-work construction, not a claim about measured network throughput. All processors' instructions and the stated delivery/in-flight work are charged. No direct-home shortcut or uncharged key oracle is used.

## 1. Contexts, channels and bounded state

Each algorithmic node has its main scalar context and one dedicated scalar receiver context sharing that node's coherent private memory. The receiver is resident, not invoked by a per-packet software dispatcher. Its code and register parameters are fixed before hop2 begins and persist across both phases. It executes independently while the main context stalls on a send or waits for completion. Both contexts' work is included.

Hop1, final-cursor, barrier and phase4 messages retain the original addressed delivery handler. Hop2 payloads have a separate fixed row-payload endpoint/port bound to the receiver window; control messages cannot enter it. This is a fixed channel distinction, not a content-dependent classification operation. The extra context, stack,48-byte packet window, valid/ack state and payload-port state require less than2KiB per node, included in the fixed-region reservation. The existing sender ports and their in-flight state remain charged separately.

Use the source blocking-port convention: a sender's completing STP waits if its bounded port/window cannot accept another message; the receiver's first LDP waits until a complete packet is atomically published. A waiting instruction counts once, not repeatedly as a poll. Arbitration is fair, backpressured and the same transport mechanism as the source delivery construction. There is one bounded receiver window, no hidden unbounded software queue, and no software context-switch loop. The published packet remains immutable until acknowledgment. Packet publication and the blocking device load have acquire ordering; the window is an explicitly ordered device endpoint, not an ordinary unsynchronized shared-memory buffer.

## 2. Receiver body and acknowledgment

Resident inputs are the packet address, d3 mask, cursor/end-table base, d4 mask, counter base, acknowledgment endpoint, and constants1,32,-32. The body reads the four payload words in two LDPs. It computes d3, addresses the interleaved cursor/end pair, reserves32bytes with LDADD, loads the end and compares. On a full bucket it undoes the increment. Otherwise it stores the two payload pairs, computes d4 and increments its counter.

Both exits execute STLR xzr,[ack] after those writes, then branch to the blocking first load. STLR is a base scalar release store with an unindexed address; its encoding is restricted in the interpreter. The finite callable test uses RET in place of that final branch, with the same one-instruction price. The accepted path costs15 instructions; a full-bucket drop costs11. Charge15 for every potential record, including ones that in fact never reach hop2 because of earlier drops.

The explicit release acknowledgment is charged in addition to the source's two per-message in-flight updates. In-flight decrement occurs only after the acknowledgment, after insertion/drop and the d4 update are visible. Delivery completion means processing completion, not inbox arrival. The main context's completion wait observes that acknowledgment with acquire semantics; its DSB/barrier cannot pass merely because the packet reached the window. This release/completion contract is a material part of the machine definition.

There is no per-message call, argument setup, register save/restore or polling: parameters stay in the dedicated context. Startup and public parameter/address setup are bounded by the separate2^14 per-node phase/context-control allowance. A final issued blocking receive after the last packet costs at most one additional instruction and is included there; with no further messages it executes no further work. Program, context and port initialization also remain within the charged2^20 node initialization row. These bounds do not excuse an implementation that adds a software scheduler; such an implementation would require its own extra ledger.

## 3. Phases and progress

Partition absolute X indices into two contiguous phases, each containing at mostceil(q1/2) indices. Keep absolute k in w2. For each phase initialize only the hop1 cursor table; emit records and final-cursor messages; drain and barrier. Then h2s streams those local slots through the dedicated row-payload port. Extract destination d2 before BFI replaces w3bits0..34 with original source column j. Sign/l and w2 are untouched.

The receiver never sends another payload or waits on main-context work while processing a packet: it only performs bounded local operations and acknowledges. Fair execution therefore drains finite sender workloads even if main contexts block on their ports. After h2s, each sender waits for every acknowledgment and all nodes barrier. Only then may hop1 slot storage be reused. Buckets and d4 counts persist through both phases. Prefix4, full-record pass4 and scanning start only after the second phase fully drains.

The two first-hop cursor initializations and finishing sweeps are charged. h2s costs7n+2floor(n/8)+11*SIDE+2 per phase. Bucket cursor/end initialization costs5*2^26+1. After the phases, full-record pass4 costs6n+2floor(n/8)+11*2^26+1. Prefix and scanner costs are unchanged. Sixty-four global barriers at1024instructions/node and2^14 additional phase/context-control instructions per node conservatively cover the finite public phase transitions, arguments and completion setup. Receiver record work is not hidden inside those fixed allowances.

## 4. Addresses and resource obligations

Reserve[0,2^52) for fixed memory; add the d3 end table and2KiB extra receiver state explicitly. Place phase receive slots, tight32*CAP_3 full-record buckets, then the8-byte pointer array in disjoint increasing regions. success_probe.py checks the exact endpoint below2^64, fixed reservation, positive cap gaps and attempted hop1 cursor bound. The bucket increment/rollback cannot wrap because its end is strictly below the subsequent pointer-array region and the final address limit. No second-hop R2 receive array or final-cursor table is needed. The bound does not rely on reusing dead storage.

Finite tests run actual emit, first-hop delivery, final cursors, sender, receiver acknowledgment, prefix, scatter and scanner over sixteen nodes and two phases, with filtered-free high raw labels, slot/bucket caps and independent reference contents/counts. They establish the stated finite behavior and instruction counts, not a physical asynchronous implementation. The protocol and all ideal-machine assumptions remain explicit review obligations.
