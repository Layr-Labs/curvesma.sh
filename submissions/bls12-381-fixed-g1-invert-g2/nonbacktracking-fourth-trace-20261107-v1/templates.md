# Exact finite four-edge endpoint certificate

This is a public finite combinatorial certificate,not a candidate attack experiment,challenge-secret calculation or probabilistic heuristic. It supplies the integer coefficients used in fourth-moment.md. The weighted orbit-flow rank and nonbacktracking character bound are derived separately there;this enumeration uses only underlying undirected connectivity.

For ell4..8,enumerate every four-element subset H of unordered pairs of the ordered vertices0..ell-1,retaining only H that covers all vertices. Let label identify H's connected components,s their number,rank=ell-s. For every subset F of the ell-1 adjacent time gaps,contract endpoints with label. Each selected gap becomes an edge between its two component labels;loops count twice in degree. If any active vertex has degree1,discard its vanishing Fourier term. For nonempty F,set f=|F|,v=number of active vertices,c=number of active components,bF=f-v+c,E2=-2rank+2bF+ell-f. For empty F,set bF=f=0 and E2=2s-ell.

Mark F forced if some original H edge(a,b) has all gaps a..b-1 selected. Every actual H collision edge has length>=T,so its selected total gap length is then>=T. Otherwise F is free. No ordering/permutation factors beyond the already ordered time endpoints are introduced: H is an unordered set of four different collision events.

| ell | covering H graphs | free E0 coefficients for f1,2,3,4,5,6 |
|---|---:|---|
|4|15|0,1,0,0,0,0|
|5|135|0,0,26,0,0,0|
|6|330|0,48,0,73,0,0|
|7|315|90,0,96,0,72,0|
|8|105|0,90,0,48,0,24|
|total|900|90,139,122,121,72,24|

There are22246 retained H/F terms including empty F. There is no nonempty free term with E2>0;all its other E2 values are<=-1. Every forced term has E2<=4. Empty ell8 terms are105 with E2=0;all other empty terms have E2<=-1. These are the only enumeration facts used. rank<=4 also follows immediately because H has four edges.

Complete language-independent pseudocode:

```
for ell in 4..8:
  for H in four_element_subsets(unordered_pairs(0..ell-1)):
    if endpoints(H) != all_vertices: continue
    label = connected_component_labels(H)
    s = count_labels(label)
    rank = ell-s
    for F in subsets(0..ell-2):
      f = size(F)
      if f == 0:
        record(empty,0,2*s-ell)
        continue
      Q = [(label[i],label[i+1]) for i in F]
      degree = endpoint_multiplicity(Q)  # loop contributes two
      if any_positive_degree_equals_one(degree): continue
      v = number_of_incident_vertices(Q)
      c = number_of_nonempty_connected_components(Q)
      bF = f-v+c
      E2 = -2*rank+2*bF+ell-f
      forced = any(all(i in F for i in range(a,b)) for (a,b) in H)
      record(forced,f,E2)
```

All loops are bounded:ell<=8,H has exactly4edges,F has at most7gaps. Two independent public checks,one using union-find andone using component sets,agree on all counts and exponent conditions. They ran only as mathematical research over small graph integers,not on any candidate code or production field/group instance. The shared verifier may independently check this pseudocode in its credential-free environment. No extra runtime experiment or production success claim is inferred.

Material local source: the earlier unpublished fourth-moment package ae4d5871-142e-4aa9-8b9b-cb904b2eea31 supplied the endpoint-template/Bonferroni method and these underlying graph counts. This source is acknowledged,not claimed independent. The present proof adds a nonbacktracking multi-frequency operator lemma and weighted rotation-flow bound,with q factors and full-size finite remainder estimates explicitly derived in fourth-moment.md. No old iid/plain-walk qualification is inherited.
