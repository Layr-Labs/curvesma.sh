# Static interpolation coefficient certificate

This is a mathematical integer-DAG certificate and a text-identity check,NOT execution of candidate assembly/Python or a production experiment. The actual submitted archive was fetched from evaluated commit832a8782b89a1754831bbcc98e662f9ca28563c3. Its project_toom.S SHA256 is6fba18a21e7602bcf9c04803ce6f5003df1300ffc3d9f132bc5ba0b411093d2e,identical to the local submittedpackage. No archive byte was changed to obtain this evidence.

The v38 rejection said temp1920=2Q1+Q1=3Q1. The archived lines actually are:

    FADD x21,1824,x20,384,x20,384      // 2Q1
    FADD x21,1872,x21,1824,x21,1824   // 4Q1
    FADD x21,1920,x21,1872,x20,384    // 4Q1+Q1=5Q1
    ...
    FADD x19,0,x19,0,x21,1920        // adds5Q1 to a0r

The middle address1872 is material. This contradicts the stated3Q1 rejection reason;it does not itself confer acceptance. This newrevision leaves all arithmetic instructions unchanged and adds annotations/certificate for renewed complete review. No previous qualification/score is inherited.

x20 holds P0..P6,Q0..Q6 at offsets48i;Q group was formed from the two shifted product groups by seven adds. x21 is the disjoint multiple area. Each rowbelow is exactlyONE nativeFADD and declaresits independently checked integer multiple. All61 source rows were compared as plain TEXT with the corresponding integer-DAG design;no instructions interpreted or executed. The integer output vectors were summed separately and compared with Lagrange columns derived from N_i,V modulo t4−2t2+2,not merely with cached desiredlabels.

| Node | x21 offset | Input | Multiple | Native source addresses added |
| --- | ---: | --- | ---: | --- |
| 0 | 0 | P0 | 2 | x20+0 + x20+0 |
| 1 | 48 | P0 | 3 | x21+0 + x20+0 |
| 2 | 96 | P0 | 4 | x21+0 + x21+0 |
| 3 | 144 | P0 | 8 | x21+96 + x21+96 |
| 4 | 192 | P0 | 9 | x21+144 + x20+0 |
| 5 | 240 | P1 | 2 | x20+48 + x20+48 |
| 6 | 288 | P1 | 4 | x21+240 + x21+240 |
| 7 | 336 | P1 | 8 | x21+288 + x21+288 |
| 8 | 384 | P1 | 5 | x21+288 + x20+48 |
| 9 | 432 | P2 | 2 | x20+96 + x20+96 |
| 10 | 480 | P2 | 4 | x21+432 + x21+432 |
| 11 | 528 | P2 | 8 | x21+480 + x21+480 |
| 12 | 576 | P2 | 16 | x21+528 + x21+528 |
| 13 | 624 | P2 | 3 | x21+432 + x20+96 |
| 14 | 672 | P2 | 6 | x21+624 + x21+624 |
| 15 | 720 | P2 | 12 | x21+672 + x21+672 |
| 16 | 768 | P2 | 13 | x21+720 + x20+96 |
| 17 | 816 | P3 | 2 | x20+144 + x20+144 |
| 18 | 864 | P3 | 4 | x21+816 + x21+816 |
| 19 | 912 | P3 | 5 | x21+864 + x20+144 |
| 20 | 960 | P4 | 2 | x20+192 + x20+192 |
| 21 | 1008 | P4 | 4 | x21+960 + x21+960 |
| 22 | 1056 | P4 | 5 | x21+1008 + x20+192 |
| 23 | 1104 | P4 | 3 | x21+960 + x20+192 |
| 24 | 1152 | P4 | 6 | x21+1104 + x21+1104 |
| 25 | 1200 | P4 | 7 | x21+1152 + x20+192 |
| 26 | 1248 | P5 | 2 | x20+240 + x20+240 |
| 27 | 1296 | P5 | 3 | x21+1248 + x20+240 |
| 28 | 1344 | P6 | 2 | x20+288 + x20+288 |
| 29 | 1392 | P6 | 4 | x21+1344 + x21+1344 |
| 30 | 1440 | P6 | 3 | x21+1344 + x20+288 |
| 31 | 1488 | P6 | 8 | x21+1392 + x21+1392 |
| 32 | 1536 | P6 | 9 | x21+1488 + x20+288 |
| 33 | 1584 | Q0 | 2 | x20+336 + x20+336 |
| 34 | 1632 | Q0 | 4 | x21+1584 + x21+1584 |
| 35 | 1680 | Q0 | 3 | x21+1584 + x20+336 |
| 36 | 1728 | Q0 | 8 | x21+1632 + x21+1632 |
| 37 | 1776 | Q0 | 9 | x21+1728 + x20+336 |
| 38 | 1824 | Q1 | 2 | x20+384 + x20+384 |
| 39 | 1872 | Q1 | 4 | x21+1824 + x21+1824 |
| 40 | 1920 | Q1 | 5 | x21+1872 + x20+384 |
| 41 | 1968 | Q2 | 2 | x20+432 + x20+432 |
| 42 | 2016 | Q2 | 3 | x21+1968 + x20+432 |
| 43 | 2064 | Q2 | 6 | x21+2016 + x21+2016 |
| 44 | 2112 | Q2 | 7 | x21+2064 + x20+432 |
| 45 | 2160 | Q2 | 4 | x21+1968 + x21+1968 |
| 46 | 2208 | Q2 | 8 | x21+2160 + x21+2160 |
| 47 | 2256 | Q3 | 2 | x20+480 + x20+480 |
| 48 | 2304 | Q3 | 4 | x21+2256 + x21+2256 |
| 49 | 2352 | Q3 | 3 | x21+2256 + x20+480 |
| 50 | 2400 | Q3 | 6 | x21+2352 + x21+2352 |
| 51 | 2448 | Q4 | 2 | x20+528 + x20+528 |
| 52 | 2496 | Q4 | 4 | x21+2448 + x21+2448 |
| 53 | 2544 | Q4 | 8 | x21+2496 + x21+2496 |
| 54 | 2592 | Q4 | 3 | x21+2448 + x20+528 |
| 55 | 2640 | Q4 | 6 | x21+2592 + x21+2592 |
| 56 | 2688 | Q4 | 5 | x21+2496 + x20+528 |
| 57 | 2736 | Q5 | 2 | x20+576 + x20+576 |
| 58 | 2784 | Q6 | 2 | x20+624 + x20+624 |
| 59 | 2832 | Q6 | 3 | x21+2784 + x20+624 |
| 60 | 2880 | Q6 | 6 | x21+2832 + x21+2832 |

Output coefficient rows in P0..P6,Q0..Q6 order:
[3, 8, 16, 1, 5, 0, 2, 2, 5, -7, 4, -8, 1, -3]
[9, 4, 8, -1, -5, 0, -4, -4, 0, -6, -6, 6, -2, 6]
[-1, 5, -13, -1, -1, -1, 3, 3, 4, 8, 1, 5, 0, 1]
[-3, -5, 1, -5, 7, -3, 9, 9, 2, 4, -1, -5, 0, -2]

All21 FMUL argument triples exactlymatch U0V0,U1V2,U2V1;all7 Q-sum FADDs andall4FCOPY/47outputFADD/FSUB rows exactlymatch the checked design text. Canonical full-p field add/sub realizes these integer linear identities modulo p. Source core/frame/ABI/canonical carry proofs remain in proof.md. This certificate certifies the stated symbolic linear circuit only,not H/R/storage/production inversion or an observed runtimepass. Shared llm-verifier still owns review/scoring.

This revision uses the same arithmetic macro/instruction token sequence. Its annotated project_toom.S SHA256 is58ec4dfffdba3a89da0605a0e8cb394bbfaacc95a9d5546206465fd750dc41e1. Added comments name every positive-multiple slot;they execute no instruction and change no frame/count.

For this revision,project_toom.S is a retained reference circuit,not called by the algorithm. Its added REFERENCE comment changes the file hash but no coefficient/instruction tokens. Actual combine_products.S implements the identical61-node interpolation after the twenty paid lookup products and one pinned Montgomery product. The old annotated-file SHA above identifies the earlier file version only.

Currentfile scope: the oldproject_toom.S/combine_products.S are archivedsource references,NOT included orcalled here. The actual lazy_interpolate.S carries the SAME115integerarithmetic/4copy coefficientcalls under an explicit invertible48-to64slot mapping. Oldfilehashes above identify earlierpublic/referencearchives,not missingcurrentfiles. Allnativehotwork andfinalcorrections are in the actual newsource;quotient_certificate.md/lazy_interpolation.md prove its semantics/counts.
