# FILTER20 correction

In FILTER20.md replace exactly:

m* x (1 + 0.0037): the upper 2-SE bound of filtfast's paired understatement of exact N.

with:

The quoted excess-1.0 paired differences give a two-SE upper understatement of max(0, -1.18% + 2*0.51%, 0.11% + 2*0.57%) = 1.25%, not 0.37%. Apply 1.0125 to the filtfast scale 2^26 estimate only; scale 2^25 uses exact filt10. The two mean-plus-two-SE values are 0.009309177803458114 and 1.0125*0.008984999557547182 = 0.009097312052016521. Their maximum is below 0.0093092. The variance maximum also stays at scale 2^25: max(4.10363742714952e-08, 1.0125^2*2.034198985394874e-08) < 4.1037e-08. The retained m* = 0.0093092*1.0037 uses 1.0037 solely as additional safety margin; it is not an estimate of sampler understatement. Thus sqrt(m*^2 + 4.1037e-08) < 0.009346, and the existing NU, size cap, ledger and claim are unchanged. This remains H3's extrapolation from the two largest simulated scales, not a measured production bound.
