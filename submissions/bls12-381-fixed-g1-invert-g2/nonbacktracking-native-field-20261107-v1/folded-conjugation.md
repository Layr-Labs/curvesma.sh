# Fold conjugation into signed structured constants

This is an exact scalar implementation refinement of the same canonical function. It changes no fieldvalue,canonicalorientation,hashlabel,trajectory,probabilitymodel orworkcap. Allconstants aretarget-independent publicfieldvalues constructedcold;theirpositive andnegative Montgomeryencodings arebothstoredandcharged.

Let the rawinputpairbe(a,b),andletq indicate whether conjugation isrequired. The following formulas compute c*(a+(-1)^q*b*u),with u²=-1:

- c=(s,0):even(sa,sb),odd(sa,(-s)b).
- c=(0,s):even((-s)b,sa),odd(sb,sa).
- c=(s,s):even(s(a-b),s(a+b)),odd(s(a+b),s(a-b)).
- c=(s,-s):even(s(a+b),s(b-a)),odd(s(a-b),(-s)(a+b)).

Each uses two native Fp multiplications. Real/pureimaginarycases need no otherfieldoperation;diagonalcases needexactlytwo canonicaladd/suboperations. Selectcoldconstant s or-s asindicated;do not materialize -b ornegate the output afterward. Allsign/conjugationeffectsareinthemultiplicationconstants oradd/sub choice. Input/outputpathsuseordinaryscalarcallsandpaidconstant/pointerselection,notfreefieldautomorphismprimitives. Intermediatebuffersandcanonicalwritesarechargedinthe unchangedkey/control256 andwrapper256 allowances. AllnativeFp384 interfacesand42-instructionfield-operationpricesremainunchanged.

Onlythree such productsoccur:T(x),andtheother twob-pairtransforms. Theythereforeneedatmostsixfieldoperations,replacingthepreviousnine. For thecheapk0/6branches,theother twopairsinsteadneedfourcopy/negationoperations;T(x) needsatmosttwo,so6stillboundsallbranches. Otherevenk needsatmosttwo b1signoperationsplusT(x)'stwo,andtheotherrealconstantproduct needszeroadd/sub. Candidate-generationrotationsstilluse16operations,withtheexplicit96-byteinitialcopyinside384minimum/copy asbefore. Noothercopyornegationisomitted.

The exactnonmultiplicationledgerbecomes:
 6*42 +16*42 +384minimum/copy +256key/control +256wrappers
 =252+672+384+256+256=1820.
Worstcanonicalhelper10*384+1820=5660. Its orientation-dependenthardbound is
 f(k)=4124+768c(k), c=0fork0/6,1forothereven,2forodd.
WiththeunchangedB=floor1027K/768 counter,
 4124K+768B<=5151K.
Theexistingstationaryvariance/Chebyshevproofstillgivescapfailure<2^-100. Everyfailedpartialhelperispaidcold. No expected-only work ornewheuristic enters.

Completeordinaryledger:
 5151canonical+807tabulation+797jointlabelhash+46scalar
 +848controller+192head+224payload+32cappedtraversals+1retry
 =8098K.
Retain8128K,30slack. CompletecoldD=540*2^119 alreadyincludesallfixedconstants/code,pairings,hashes,generatorconstruction,rejections,headinitialization,failedfinishes,canonicalG2outputandexactverification. Additional-sconstants andfixedbranchcodearefarbelowthepreviouscompletecold2^120 term. Peakmemoryremains<2^136bytes,zeroadvice. ProtectedC20000 andone denseFp12referenceperstepgive
 W<=(28128K+540*2^119)/20000.
Exactintegerinequality(28128K+540*2^119)^1000<20000^1000*2^126173 provesproposed126173/1000=126.173,abovediagnostic126.17238483423418. Existingjointphase/dictionary/canonicalcaplossesareunchangedandsuccessremains>0.39 underH andtheexplicitfinitebankedstoragecondition. Noofficialscore,qualification,promotion orproductionexecutionisclaimed.
