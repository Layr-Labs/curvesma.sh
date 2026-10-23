# Capped linear-hash BSGS

Fixed non-identity `Q` in G2, target `z` in GT. Build `g = e(P0, Q)`, run a square baby/giant search for `a` with `z = g^a`, and join the lists by a 256-bit linear hash. Every reported match is an equality of the canonical 72-limb GT encoding. A hard cap on those exact comparisons makes the instruction total independent of collision frequency. The score rational is `128137645/1000000`.
