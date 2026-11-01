# Two-phase direct full-record receiver: machine, work and memory

## Scalar contexts and ordered endpoints

Each of N=2^72 nodes has a main scalar context and a dedicated resident receiver context, both included in the summed computation. Private memory is below2^64 bytes. Base scalar AArch64 plus RNG/LSE is used. A separate ordinary delivery context executes the unchanged6-instruction dlv work for fixed-address messages. There are no free interrupt entries, callbacks, queues, polling loops or uncharged dispatch. The machine is an explicit instruction-equivalent construction, not a throughput claim for an existing system.

A packet endpoint exposes a fixed48-byte window: destination index/address and four payload words. Stores complete a message on the final word. A store stalls while its endpoint cannot accept the next message. For the hop2 payload endpoint only, the resident receiver's first LDP blocks until a complete packet is present. Its payload window remains unchanged until the receiver's STLR xzr,[x13] release acknowledgement. Only then may a new packet become visible. This ordered blocking register behavior is part of the machine contract; stalled instructions count once. The link's in-flight decrement occurs only after the release, which orders payload/cursor/counter writes before completion. The receiver's branch back to its entry replaces the test harness RET at the same single-instruction cost. Its x0 window, x1/x3 masks, x2/x4 table bases, x13 ACK address, x17=1,x18=32,x19=-32 and loop return address are initialized once, remain resident and do not share the main context's clobbered registers. Startup/state space is explicitly charged.

Hop1 payload and fin messages retain fixed-address dlv delivery. Their ordinary endpoint releases its window and decrements in-flight state after dlv completes its stores, as in the inherited6+2instruction delivery contract; it neither calls receive nor requires an extra software STLR. Hop2 payload uses a separate endpoint executing receive; a fin packet cannot be misclassified as payload. No hop2 fin protocol exists. The receiving main context reads d4/buckets only after all sender ports drain, release ACKs complete and the barrier ends. Buckets remain stable until pointer scatter and candidate resolution finish.

DSB waits for a node's own stores and in-flight counts to complete. Every message charges two scalar-equivalent count updates, one on acceptance and one after completed delivery/ACK. A combining barrier tree uses separate epoch/reservation words. At most five waits per node each use LDAXR,CMP,B.EQ and WFE with a single store wakeup; at most five14-instruction messages, DSB and20setup instructions give131, below the charged1024/node/barrier. Unrelated wakeups are excluded by the explicit machine. Sixty-four barriers suffice for program/Q distribution, preprocessing, both first-hop/drain phases, grouping,12resolution transitions and output.

## Exact kernel charges

- key emit:10/record, or5when the first-hop slot is full. Add ordinary delivery6 and in-flight2. Charge18for every attempted leaf.
- h2s:7/retained record (twoLDP, destination AND, source-column BFI and threeSTP),2floor(n/8) per source slot,11per slot,2fixed per invocation. There are2invocations per node. Destination d2 is extracted before BFI replaces low36 bits of w3 with source j; type36..38 and leaf index39..63 are untouched.
- receive:15on acceptance (twoLDP, mask/address, LDADD, end load, CMP/branch, twoSTP, d4mask/address/STADD, STLR, loop jump);11on drop with rollback and release. Charge15for every delivered record, plus in-flight2. No ordinary dlv is additionally assumed for these same payloads.
- pass4:6/record,2floor(n/8) per bucket,11per bucket and1fixed. It reads full records and places their addresses into d4 ranges.
- init3:5per bucket plus1, writing both cursor and end for32*CAP3tight stride.
- zero4:4per four counters plus1. prefix4:10per four counters plus1. All four grouping fixed returns (init3,zero4,prefix4,pass4) are separate ledger rows.
- scan_ranges is unchanged:195/64per record plus26per range,8per candidate and one extra full-list append attempt, plus2fixed.
- tinit and fin are unchanged:3SIDE+1 and9SIDE+3. Fin delivery/in-flight contributes8SIDE, once per first-hop phase. No second-hop cursor initialization or finishing charge exists.

Thus hop2 charges97S/4 plus2N(11SIDE+2), and grouping record work595S/64. Every failed, dropped or interrupted prefix is bounded by the full charge. Additional setup2^14/node covers receiver registers/window and twice-phase control; candidate decoding has an additional64/candidate. Arithmetic-fixup and j-path meter caps and admission reserves persist across both phases; they are never reset with hop1 cursors.

## Phase-4 exchange conservation

Relay cap2^34 cannot be justified by the old per-node2^32-message assertion. Instead bound total messages across nodes. Each request subphase has at most N*2*CAND_CAP records; each reply subphase has at most N*8*REQ_CAP. Relays forward without duplication and dropped records only reduce this total. The sum for each subphase is at most12*N*SIDE+88*global_records. Since CAND_CAP=2^28,REQ_CAP=2^29,SIDE=2^36, this is below N*2^41. Charge four such complete subphases. Per-source/request/relay success losses still use their separate caps, not an average cap. This conservation does not amortize across targets.

## Full fixed allocation and disjoint layout

Let q=145289002,u=2^25,s=2^36,b=14755795394. Reserve

F=2^20+5*2048+1024+288q+160u+32(q+u+1)+256u
 +2*64+8s+32s+16*2^24+8*2^48+16*2^28
 +64*2^29+16*2^28*64+64*2^34+96s+2048+2^34.

Terms cover code/stack, barriers/header, X/Y/scalars/batch, ports and hop1 cursor/E1 tables, interleaved bucket cursor/end and d4/candidate arrays, request/reply buffers, relay buffer, phase4tables, receiver state and explicit alignment padding. F<2^51+2^44<2^52. The2^33-byte batch remains in the lower half of a2^34-aligned block, preserving the inherited leaf-pointer encoding. No extra per-leaf metadata is needed.

Fixed reservation is[0,2^52). Next is hop1 R1 of7*2^60bytes. Then32*2^24*b full-record bucket bytes and8*2^24*b pointer bytes. Their endpoint E=17977400794953023488<2^64. A simple bound is b<14*2^30, giving E<2^52+63*2^58<2^64. Total memory N*E<2^136. All regions remain allocated; no live-storage overlay is needed. There is no R2 or E2 and no old2^38bucket-mask constraint.

A phase has qphi=72644501<2^27 batches and at most6qphi*2^25attempted records per node. Even if every attempt selects one destination, the cursor is below2^52+32*CAP1(s-1)+32*6qphi*2^25<2^52+7*2^60+3*2^58<2^63. Only in-cap addresses are written; fin clamps attempted cursors. A full bucket's temporary reservation is at most bucket_end+32<=pointer_base+32<E and is rolled back before ACK. Total d4counts are at most2^24*b<2^58. No address or counter wraps.

## Proof transfer boundary

Only the first-hop phase-slot mean/variance and one rather than two slot-loss terms change. Full buckets see the whole run. Original row order, allgenerator coins, arithmetic/label computations and persistent meter guards remain identical. The two-fraction label-error row is retained. No transfer of an incompatible orbit-min key or unsupported clean-run conditioning is made. The finite emulator checks instruction semantics and synchronous ordering; it does not certify real asynchronous hardware. That limitation and this abstract endpoint contract are explicit parts of review.
