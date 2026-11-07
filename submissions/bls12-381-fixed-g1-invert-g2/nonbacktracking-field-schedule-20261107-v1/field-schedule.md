# Explicit canonical field-operation schedule

This refines only deterministic implementation bounds. Group values,canonical order,labels,hash bridge,workcap andall success calculations remain the same. The prior pointer-minimum leaf remains327 internal instructions;there is no candidate execution on the host.

## Twelve candidate values:exactly sixteen field operations beyond six products

Let x be the selected nonzero Fp2 pair,y=T(x),b=beta*x,c=beta*y. The six native Fp multiplications are two fory,two forb,two forc. Write outputs y,b,c directly into candidate slots k1,k2,k3. Copy x toslot k0 withsixLDP/sixSTP,12 scalar data instructions plusatmosttwo pointer instructions. These14 instructions fit the retained384minimum/copy allowance:327minimum+16interfaces+14=357<384. No additional Fp copy price is needed or omitted.

For the even slots from x,b,write:
 k4=b-x, k6=-x, k8=-b, k10=x-b.
For the odd slots from y,c,write:
 k5=c-y, k7=-y, k9=-c, k11=y-c.
Each assignment is two canonical Fp additions/subtractions/negations,one perreal/imaginary coordinate. Eight assignments therefore use exactly16 field operations,not the old19-operation allowance. Each costs<=48 scalar instructions including its concrete interface,as infield-ops.S/copy-negate.S. All results arewritten directly into the contiguous1152-byte candidate arena before canonical-min.S reads them. Everyarena byte used isinitialized bythese writes;nozero initializationorfreecopyisassumed.

## Structured Fp2 constants:at most nine further field operations

The only three structured products are T(x) andtheother two b-pair transforms. All constants c[k,d] are real,purely imaginary,or(s,s)/(s,-s). Their signedreal constants areprecomputedcold,including-s. For one such product on(a,b),firstperform optionalconjugation b'=-b (oneFpnegation),oruseb'=b for even k.

For c=(s,0),outputs(sa,sb') use twoFp products andnoadditionalfield operation. For c=(0,s),outputs((-s)b',sa) use thecold-s constant,twoFp products,andnoadditionalfield operation. For c=(s,s),outputs(s(a-b'),s(a+b')) use twoadd/sub operations andtwoFp products. For c=(s,-s),outputs(s(a+b'),s(b'-a)) likewise use twoadd/sub andtwoFp products. Including optionalconjugation,every structured product uses atmostthreefieldoperations,notfour. Its copies/argument/control areinthe unchangedkey/control768 andwrapper512 allowances andnativeFp416 interfaces. AlloutputsarecanonicalMontgomeryvalues;nocanonicalizationinverseorresidueconversionisinserted.

Thus allthree products need atmost9 fieldoperations,not theold12. For cheapk0/6 classes,theother twoFp2pairs arecopied/negated inatmostfour48-instructionfield-operation allowances;T(x) stillusesatmostthree,so7<=9. For otherevenk,b1 isjustsigned andtheotherpair isreal-scaled,also<=9. Both ordinary andcheapbranches therefore fit the same fixed allowance. No field-operation omission ishiddenby theweightedclass cap.

## Complete revised certificate

The nonmultiplication allowance is now2864:
 9*48 structured-productfieldops +16*48 rotations
 +384minimum/copy +768otherkey/control +512wrappers
 =432+768+384+768+512=2864.
It replaces3152 (oldsixextrafield/rotationoperations included288 more). Native Fp products remain416 each,withreferenceinternalsnotdoublecounted. Theworsthelper is10*416+2864=7024. Its orientation-dependentbound is
 f(k)=5360+832c(k), c=0fork0/6,1forothereven,2forodd.
The existing128-bit hardclasscounter B=floor1027K/768 isunchanged,and
 5360K+832B<=6474K.
Its proven capfailure remains<2^-100. Allpartialfailedhelpersarepaidfixedcold. Noexpected-onlycostornewprobabilityhypothesisisused.

Completeordinaryledger:
 6474canonical+807tabulation+797jointlabelhash+46scalar
 +848controller+192head+224payload+32cappedtraversals+1retry
 =9421K.
Retain9440K,19slack. CompletecoldD=540*2^119 stillpaysallheadzeroing,generator/table/scalar construction,hashes,rejections,code/constants,validation,pairings,failurefinishes,verification andcanonicalG2output. Peakstorage<2^136bytes,zeroadvice. Therefore
 W<=(29440K+540*2^119)/20000.
Exactintegerinequality(29440K+540*2^119)^1000<20000^1000*2^126239 supportsproposed126.239,abovediagnostic126.23814349539336. H andtheexplicitfinitebankedRAMmodelremainmaterial;thisisnotanofficialscore,qualification orproductionexecution.
