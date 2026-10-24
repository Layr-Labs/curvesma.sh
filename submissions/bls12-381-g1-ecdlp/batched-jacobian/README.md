# Width-5 wNAF scalar recurrence with blocked Montgomery normalization

Recovers the G1 discrete logarithm from the supplied d-th power with
d = 2^20. Two-stage Cheon search over j in Z_m, m = (r-1)/2^20. Baby
records hold a 96-byte uncompressed affine point and a two-limb index
at offset 96 (112-byte records; mergesort moves them whole). Each
recurrence step recodes its scalar as a width-5 wNAF: at most 255
doublings and 64 nonzero-digit additions plus a 16-operation
odd-multiple table. Lists are normalized to affine form by blocked
Montgomery inversion: 12 field multiplications per point and one
Fermat inversion per 2^15-point block.
