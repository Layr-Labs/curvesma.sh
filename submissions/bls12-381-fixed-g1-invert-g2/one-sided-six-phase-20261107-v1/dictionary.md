# One ordinary hash: joint generator label and capped dictionary

The global dictionary stores every visited orbit and terminates at the first repeat. It retains the initialized head array and hard traversal cap. This revision removes the separate ordinary PHAS hash. It uses standard BLAKE3-256 only; H is explicitly applied to its full256-bit outputs. No wider hash, changed round count or additional heuristic is introduced. Actual labels and trajectories change from the predecessor.

## Exact bit assignment and ordinary label

Compute the ordinary JUMP hash of the same104-byte message JUMP16||seed32||zero16||U40. Prefix64 len64 flags1 is cached cold; final40 len40 flags10 counter0 uses full seven-round BLAKE3-256. Read the32-byte digest as unsigned256-bit LITTLE-ENDIAN integer Y. Set
 j=Y mod2^76, V=floor(Y/2^76), e=floor(6V/2^180), bucket=V mod2^128.
No k correction or inverse-previous-edge rejection is used. V has180 bits. The phase and bucket share bits; they are NOT claimed independent. j is independent of V under H.

For little-endian64-bit digest words D0..D3, j uses D0 andlow12 D1. V has limbs V0=EXTR(D2,D1,12), V1=EXTR(D3,D2,12), V2=D3>>12. Bucket is(V0,V1). Phase uses an ordinary three-limb multiplication by6;the upper52-bit boundary of the product gives e. V2<2^52 and6V<2^183,so no256-bit unit primitive or hidden fourth wide operation is required. MUL/UMULH,ADDS/ADC andshifts fit the retained20 phase instructions. Arguments,loads andbucket copies fit the explicit96-instructioncontroller. The only ordinary compression costs776,not1594.

## Coupling to independent six-phase increments

For uniform V in[0,2^180),six phase bins I_e={V:floor(6V/2^180)=e} have sizes floor(2^180/6) orceil. Choose ideal e uniformly then V uniformly in its bin; conditional V-given-e law is unchanged. Joint total variation per query<=6/2^180, so overK+1 inputs the loss<2^-50. The independent low76 jump bits are unchanged. The ideal outgoing(j,e) is uniform independent of the group past. There are no rejected labels, initial fictitious generator or retry domains. Initial JUMP/key/bucket work is cold. Group-only orientation k remains solely for terminal twist recovery.

## Conditional bucket distributions and collision variance

Expose the ideal group trajectory and all generator-choice information (ordinary phases and indices),but not the residual V inside each phase bin. Extend the ideal chain with fresh labels after first repeat only for analysis. Group transitions inspect only j,e,not residual V. Deferred exposure therefore leaves independent residual V values,each uniform in its exposed bin. This statement is about the unaborted ideal chain,not conditioned on dictionary-cap survival. Actual execution agrees with its prefix until an explicit failure or first repeat.

Each bin has length ell>=2^177. For m=2^128,the number of integers of that bin with any fixed residue modm differs from ell/m by at most1. Hence each conditional bucket law p_i satisfies
 (1-epsilon)/m <= p_i(b) <= (1+epsilon)/m,
 epsilon=2^-47,
 since m/ell<=2^-49<epsilon. Bucket laws may differ with the exposed phases. They need not be exactly uniform or identical.

Except for the separately subtracted fingerprint failure,the distinct orbit keys have distinct U inputs. Take the first-repeat-free prefix and extend to n=K+1 independent auxiliary bucket draws with allowed phase-bin laws. Let C=sum_(i<j)1[H_i=H_j]. Conditional on all the exposed group information,
 E C <= mu=(1+epsilon)n(n-1)/(2m).
Disjoint pair indicators are independent. For overlapping pairs (i,j),(i,k),their joint expectation is at most(1+epsilon)^2/m²<=4/m² by conditioning on H_i. There are at most n³ ordered overlapping pairs. Dropping negative covariance terms only enlarges the bound. Consequently
 Var C <= mu+4n³/m² < K.
The exact protected values give mu/K<0.101 and a=floor(K/8) satisfying(a-mu)/K>0.024. Uniformly in the exposed group information,Chebyshev gives
 Pr[C>a] <= (mu+4n³/m²)/(a-mu)² <2^-100.
Average this bound over that information. This handles the phase/bucket correlation explicitly;it does not assert pairwise-independent collision indicators or use a union bound on every bucket-law perturbation.

## Dictionary implementation and hard work cap

Allocate m=2^128 heads,each16 bytes,ALL initialized to zero during charged cold work. Heads and next links are128-bit logical payload indices;0=nil,allocation starts1. Each512-byte payload contains canonicalS0..191,zero padding192..303,A304..335,B336..367,k368..375,zero376..383,U384..423,zero424..431,next432..447,zero448..511. All paired logical pointers are aligned16. Every allocated payload isinitialized before publication.

Query the bucket head,then follow next pointers. Before each old-payload visit charge the128-bit remaining traversal counter initializedfloor(K/4)+2;fail before reading beyond it. Select its bank andcompare all five U words. On different U continue;on equal U compare all24 canonicalSwords andterminate. DifferentS means fingerprint failure,equalS the exact first orbit repeat andterminal relation recovery. On a miss allocate/write payload andprepend it. Copies of live fields,head andnext areheld insharedscratch across bank switches. No unchecked hash match is accepted.

Distinct insertions traverse at most C records total. The final repeated query traverses at most its distinct co-bucket members plusits matching record,at most C+1. Thus total visits<=2C+1. On C<=floor(K/8),this is<=floor(K/4)+1,within the hard cap. No final repeated input is treated as a fresh hash. Earlier failure/cap termination only reduces work. Every execution performs at mostfloor(K/4)+2 visits;at128instructions each,this is<=32K+256. Fixed256 ispaidcold. The extra dictionary-cap failure remains<2^-100.

Per ordinary transition,head read/publication192 includes bucket/sharedcopies16,two selectors64,head loads/stores4,logicaladdress/allocation/context/control64 and44 slack. Each old-payload visit128 includes selector32,Uloads6,XOR/ORR/branch10,nextload1,two-limb traversalcounter/check10,cursor/context/control32 and37 slack. FullSconfirmation occurs at mostonce andispaidterminalcold. Payloadzero/write224 includeszero51233,selector32,Stransfers24,coefficients/phase/U/next transfers32,allocation/context/interfaces64,and39 slack. No freeallocator,memcpy ormemorysearch isused.

Every bucket head is initialized:zero8192 is512 STP zero pairs+RET. Each of2^119 blocks additionally pays24 argument,address/loop/context instructions,537*2^119<2^129. Selectors across2^92 banks cost<=32*2^92=2^97;bank-loop/setup<2^100. Add all priorcompleteconstruction/finish<2^120:totalcold<(540*2^119). Initial key/bucket/payload,all failure finishes,pairings andG2 output areincluded. No pre-zeroed memory,advice orfuturetargetamortization isassumed.

Theexactkeyisnow192-bytecanonicalFp4traceS;fp4-trace.mdprovesittestsTHE SAMEfullrawGTtwelve-orbitrelation. Thereforedictionaryfailure/Cvariance/firstrepeat/dummyextensionproofsaboveareunchanged. Payload512keepsA304,B336,k368,U384,next432;192..303ispaddinginitializedzero. FullUconfirmationcompares24Swords,thenpaidcoldrawpointreconstructionresolvesoneofthreeactualtwistsbeforetheoriginalverifiedFr/G2extraction. No directtrace-onlyanswerisaccepted. Head192,payload224,chain128remainconservative.
