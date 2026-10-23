# Unchanged production mechanism and hard caps

## Production group and algorithm: unchanged mechanism

Use protected secp256k1 y^2=x^3+7, p and prime r in `ledger.py` and
`certificate-values.json`. The matching beta and lambda satisfy cube-root
identities. Units H={+/-1,+/-lambda,+/-lambda^2} act freely outside infinity;
7 is a nonresidue, excluding x=0 points. Hence the exact quotient size is

    N = 19298681539552699237261830834781317975472927379845817397100860523586360249056.

Canonicalisation chooses least X,beta X,beta^2 X and even Y, retaining its unit u.
Infinity input returns scalar zero. Intermediate aP+bQ=O gives -a/b if b!=0,
then full target verification. Equal-x additions use a complete branch.

Generate R=2^20 independent full scalar pairs (s_j,t_j), table points
T_j=s_j P+t_j Q, and 1024 full starting pairs. Uniform scalars use rejection of
256-bit draws >=r, at most four draws per scalar. The random pairs are stored
in full, not replaced by a short-seed independence assertion. Fixed branch
h(C) is bits 64..83 of canonical X. Initial DP depth d=56, deepening only to
57 and 58; h never changes and never overlaps the DP bits.

An ordinary transition is canon(C+T_h(C)). Keep the last eight canonical classes
in an ordered ring. Test DP before escaping. On a repeat in the ring, isolate
the actual closed-cycle suffix, rewind along f to its minimum, replace that
minimum by canon(2C), clear the ring, and test DP. Ordinary/replay-of-f moves
debit one main tick, doubling TWO. Rewind never selects a minimum in the tail.
Brent compares the complete expanded (point,ordered ring) state at ordinary
boundaries, aborting exact repeats; detector checkpoints never choose escapes.
A deterministic endpoint survives coalescence unless a stated cap/detector
aborts. This does not imply independence of the doubled minimum and its basin.

The new integer budget is

    T = ceil(sqrt((98859/100000)*2^2*N))
      = 276249768022971580840080560454142166451.
    L = 40*2^58 ticks per trail.
    K = T+1024*L = 276249768022971592645996767628255200691.

Stop NEW trails after T; drain only in-flight trails, all subject to L. Global
work includes wasted escapes, rewinds, failed trails and detector work. Per-trail
length fits 64 bits, global K fits 128 bits. Labels are replayed only for key
matches: (a,b)<-u(a+s_j,b+t_j), and (a,b)<-2u(a,b) on doubling. Reconstruct
full endpoints, solve x=(a1-a2)/(b2-b1) mod r if invertible, and verify [x]P=Q.
An unequal full endpoint is a spurious key, never an answer.

## Unchanged hard memory and capped recovery

Allocate M=2^73 nodes, exactly 192 bytes each: independent 192-bit universal
fingerprint 24 B; full starting a0,b0 64 B; tick length 8 B; parent/left/right
96-bit byte pointers 36 B; color 1 B; next-two-DP-bit flags 1 B; metadata/padding
58 B. Total node pool is 3*2^79 B. Adding entries (160 B each), worker states,
code and routing scratch share a separate hard 2^30 B auxiliary cap. Total
<2^81 B, independent of every statistical assumption.

Order records by (fingerprint, initial scalars, length) in a red-black tree,
height <=148. Charge multiword addresses, transfers, comparisons and balancing
per visited level. The independent Toeplitz 192-by-256 matrix uses 447 random
bits and hashes full canonical X; DP bits do not consume fingerprint entropy.

At full pool at d=56 or 57, filter on the stored next DP bit and rebuild in
place under d+1. At fullness at d=58 fail. Bound DP handling and starts by
S=2^74 each, key replay pairs by A=2^24, complete exceptional repairs by 2^16.
Every excess fails, not reallocates. Three pool scans/two rebuilds are charged.
Depth-58 records survive both filters, so final fullness requires M distinct
depth-58 DPs even if initial fullness was common. Initial fullness is no longer
conservatively declared failure; the already specified deepening is credited.

