# Required pairing preimage

After the internally recovered scalar x has passed g^x=z, compute P=[x]P0 in
the protected G1 and set Q=Q0, the protected G2 generator. Return {"P":P,"Q":Q}.
P is null for infinity or an affine {x,y}; Q0 uses its two-component Fp2 affine
coordinates. All coefficients are exactly 64 lowercase hex digits in the pinned
bases. When x=0, P=null and Q=Q0, so the identity target is covered.

Use left-to-right binary scalar multiplication with at most 256 doublings and
256 additions, starting at infinity. Affine addition over Fp handles infinity
first. Equal x and opposite y returns infinity (including doubling with y=0).
For other equal points use slope 3*x^2/(2*y); otherwise slope (y2-y1)/(x2-x1).
Set x3=slope^2-x1-x2 and y3=slope*(x1-x3)-y1. Every actual denominator is nonzero
after the cases above. Inverse is denominator^(p-2) by binary exponentiation,
at most 512 base-field multiplies. Field primeness is protected. Point closure
and the binary invariant prove P=[x]P0 in the order-r subgroup. Optionally check
the curve equation and [r]P=infinity by the same formulas; those checks are
included. Q0's subgroup membership is a protected public constant.

Each point operation uses fewer than 1024 base-field products and 1024 adds,
including inversion, loop control and copies, under 2^32 scalar instructions.
At most 1024 point operations for multiplication and its order check cost below
2^42; curve tests, canonical serialization, root handling and GT verification
fit below 2^46, within machine.md's 2^50 fixed budget. Do not pretend a BN254
point operation costs one reference unit. All actual arithmetic is charged.

By the protected pairing's bilinearity, e([x]P0,Q0)=e(P0,Q0)^x=g^x=z. Thus a
second full pairing evaluation is unnecessary: the already performed scalar
check plus this exact algebraic reduction verifies the output. This argument
does not assert equivalence of complexity or scores between tracks.
