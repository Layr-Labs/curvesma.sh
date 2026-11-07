# Finite fourth-moment graph certificate

This is a finite combinatorial calculation, not a production attack or extra heuristic. The enumeration has no challenge-secret or hash input. It determines only graph coefficients in the analytic probability bound. No candidate program is executed locally.

For m=4,...,8, enumerate each set G of four distinct pairs from the m ordered vertices0,...,m-1. Retain only graphs incident to every vertex. Compute equality components using union-find; their number is q, rank=m-q. For every subset F of the m-1 adjacent time gaps:

1. Contract endpoints using G's component labels.
2. For each selected gap i, add an edge between the labels of i and i+1. Count a loop twice in degrees. If any nonzero degree is one, discard this zero Fourier term.
3. Set f=|F|. For nonempty F, let v be the number of incident contracted vertices and c their number of connected components; b=f-v+c. For empty F use b=0.
4. Compute twice the exponent: e2=-2*rank+2*b+m-f.
5. Mark `forced` exactly when some(a,b) in G has every gap a,...,b-1 selected in F. This requires that original equality's total selected gap length is at least T. Otherwise mark `free`.
6. Count by m,f,e2,forced. Empty F is recorded separately.

This fully specifies every finite loop and test; union-find uses ordinary undirected connectivity, while Fourier orientations can be arbitrary. No interval permutations or extra ordering multiplicity is inserted: G is an unordered set of four equality events on already ordered time endpoints.

| m | covering four-edge graphs | free E=0 counts by f=1,2,3,4,5,6 |
|---|---:|---|
|4|15|0,1,0,0,0,0|
|5|135|0,0,26,0,0,0|
|6|330|0,48,0,73,0,0|
|7|315|90,0,96,0,72,0|
|8|105|0,90,0,48,0,24|
|total|900|90,139,122,121,72,24|

There are22246 retained G/F terms including empty terms. No nonempty free term has e2>0. Every other free term has e2<=-1. Every forced term has e2<=4. Empty m=8 terms number105 and have e2=0; all other empty terms have e2<=-1. These are the only enumeration facts used by probability.md.

As a reproducibility aid, the complete enumeration in language-independent pseudocode is:

```
for m in 4..8:
  for G in four-element subsets of pairs(0..m-1):
    if endpoints(G) != all m vertices: continue
    label = connected_component_labels(G)
    q = number_of_labels(label)
    rank = m-q
    for F in subsets(0..m-2):
      if F empty: record(empty,0,2*q-m); continue
      edges = [(label[i],label[i+1]) for i in F]
      degree = endpoint_multiplicity(edges) # loop counts twice
      if any positive degree equals 1: continue
      v = number_of_incident_vertices(edges)
      c = number_of_nonempty_connected_components(edges)
      b = len(F)-v+c
      e2 = -2*rank+2*b+m-len(F)
      forced = any(all(i in F for i in range(a,b)) for (a,b) in G)
      record(forced,len(F),e2)
```

The index b in the final comprehension is the original edge endpoint; the cycle dimension named b above can equivalently be named cycle_rank to avoid shadowing. This pseudocode supplies no executable candidate experiment and makes no request to run it on the host. The bounded shared verifier may independently check the finite certificate in its credential-free execution environment. The proof is conditional only on the already disclosed H and good-table bound, not on this enumeration being a probabilistic model.
