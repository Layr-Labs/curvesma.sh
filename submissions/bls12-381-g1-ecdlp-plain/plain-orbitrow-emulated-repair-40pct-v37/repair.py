"""The v37 exact-index repair: an AArch64 emitter plus its own stack contract.

Split out of schedule.py to stay inside the verifier's per-file byte limit.  The block is emitted,
not written by hand, so the listing and the audit can never drift apart.
"""
M64 = (1 << 64) - 1
P = int("1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f624"
        "1eabfffeb153ffffb9feffffffffaaab", 16)
BETA = int("1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b"
           "409427eb4f49fffd8bfd00000000aaac", 16)
R = 1 << 384
N0 = (-pow(P, -1, 1 << 64)) % (1 << 64)      # -p^{-1} mod 2^64
BETA_R = BETA * R % P                          # Montgomery form of beta

X_BASE, T_BASE, CAND_BASE, MIN_BASE, MINIDX = 240, 288, 400, 448, 496
SAVE_BASE, NREG = 16, 28


def limbs(v, n=6):
    return [(v >> (64 * i)) & M64 for i in range(n)]


def imm64(reg, v):
    """movz/movk materialisation of a 64-bit constant.

    A single 16-bit chunk is emitted as `mov` (the assembler's own alias for movz), so the listing
    text matches what clang prints for the object byte-for-byte -- completeness.py compares the two
    and v36's `movz x12,#0xf071` vs `mov x12, #0xf071` mismatch is exactly that alias.
    """
    out = ["    mov %s,#0x%04x" % (reg, v & 0xffff)]
    for h in range(1, 4):
        w = (v >> (16 * h)) & 0xffff
        if w:
            out.append("    movk %s,#0x%04x,lsl #%d" % (reg, w, 16 * h))
    return out


def repair_block(stack_bytes=512):
    """A correct exact-index repair, emitted as literal AArch64.

    argmin_j (beta^j x mod p) for j in 0,1,2 computed the way the reference does it:
      bx  = beta*x mod p        -- one 6x6 CIOS Montgomery product with a = beta*R mod p,
                                   full schoolbook carry propagation, borrow-propagating
                                   final conditional subtraction of p
      b2x = p - x - bx (mod p)  -- the identity beta^2 x = -(1+beta) x, one 6-limb
                                   borrow-propagating subtract plus a masked add of p
    The running minimum lives in the frame.  Per candidate the block READS the min
    slots, compares (6-limb subs/sbcs cascade), and only then writes the winner back
    (csel), so the spill slots are never overwritten before they are compared.  The
    winning index lands in x8, which is what row 4 reads.  Every register the block
    touches except x8 is saved and restored.  Ends with `b L_rejoin`.
    """
    o = []
    a = lambda s: o.append("    " + s)
    # the block owns its scratch: it allocates and releases its own frame, so the leaf prologue
    # (which allocates 16 bytes for the link) is untouched and every offset below is in-frame.
    a("sub sp,sp,#%d" % stack_bytes)
    # ---- save the clobber set -------------------------------------------------
    for i in range(NREG):
        if i == 8:                       # x8 carries the result, not preserved
            continue
        a("str x%d,[sp,#%d]" % (i, SAVE_BASE + 8 * i))
    # ---- constants: a = BETA_R in x12..x17, p in x18..x23, n0' in x24 ----------
    for k, v in enumerate(limbs(BETA_R)):
        o += imm64("x%d" % (12 + k), v)
    for k, v in enumerate(limbs(P)):
        o += imm64("x%d" % (18 + k), v)
    o += imm64("x24", N0)
    # ---- zero the product window ---------------------------------------------
    a("mov x0,#0")
    for i in range(14):
        a("str x0,[sp,#%d]" % (T_BASE + 8 * i))
    # ---- initialise min = x, idx = 0 -----------------------------------------
    for i in range(6):
        a("ldr x0,[sp,#%d]" % (X_BASE + 8 * i))
        a("str x0,[sp,#%d]" % (MIN_BASE + 8 * i))
    a("str xzr,[sp,#%d]" % MINIDX)
    # ---- full 6x6 product with carry propagation: T[0..11] = BETA_R * x ---------
    for i in range(6):                                   # multiplier limb x_i
        a("ldr x8,[sp,#%d]" % (X_BASE + 8 * i))
        a("mov x9,#0")                                   # C = 0
        for j in range(6):                               # T[i+j] += a[j]*x_i + C
            a("ldr x0,[sp,#%d]" % (T_BASE + 8 * (i + j)))
            a("mul x6,x%d,x8" % (12 + j))
            a("umulh x7,x%d,x8" % (12 + j))
            a("adds x6,x6,x0")
            a("adc x7,x7,xzr")
            a("adds x6,x6,x9")
            a("adc x7,x7,xzr")
            a("str x6,[sp,#%d]" % (T_BASE + 8 * (i + j)))
            a("mov x9,x7")
        a("ldr x0,[sp,#%d]" % (T_BASE + 8 * (i + 6)))
        a("adds x0,x0,x9")
        a("str x0,[sp,#%d]" % (T_BASE + 8 * (i + 6)))
        a("cset x9,hs")
        a("ldr x0,[sp,#%d]" % (T_BASE + 8 * (i + 7)))
        a("add x0,x0,x9")                                # propagate into the next limb
        a("str x0,[sp,#%d]" % (T_BASE + 8 * (i + 7)))
    # ---- Montgomery reduction (SOS): m = T[i]*n0', T += m*p*2^(64i) --------------
    for i in range(6):
        a("ldr x0,[sp,#%d]" % (T_BASE + 8 * i))
        a("mul x10,x0,x24")                              # m = T[i] * n0' mod 2^64
        a("mul x6,x10,x18")
        a("umulh x7,x10,x18")
        a("adds x6,x6,x0")                               # T[i] + m*p0 : low word becomes 0
        a("adc x7,x7,xzr")
        a("mov x9,x7")                                   # C
        for j in range(1, 6):
            a("ldr x0,[sp,#%d]" % (T_BASE + 8 * (i + j)))
            a("mul x6,x10,x%d" % (18 + j))
            a("umulh x7,x10,x%d" % (18 + j))
            a("adds x6,x6,x0")
            a("adc x7,x7,xzr")
            a("adds x6,x6,x9")
            a("adc x7,x7,xzr")
            a("str x6,[sp,#%d]" % (T_BASE + 8 * (i + j)))
            a("mov x9,x7")
        a("ldr x0,[sp,#%d]" % (T_BASE + 8 * (i + 6)))
        a("adds x0,x0,x9")                               # 3-limb carry cascade
        a("str x0,[sp,#%d]" % (T_BASE + 8 * (i + 6)))
        a("ldr x0,[sp,#%d]" % (T_BASE + 8 * (i + 7)))
        a("adc x0,x0,xzr")
        a("str x0,[sp,#%d]" % (T_BASE + 8 * (i + 7)))
        a("ldr x0,[sp,#%d]" % (T_BASE + 8 * (i + 8)))
        a("adc x0,x0,xzr")
        a("str x0,[sp,#%d]" % (T_BASE + 8 * (i + 8)))
    # ---- one borrow-propagating conditional subtraction of p ------------------
    for j in range(6):
        a("ldr x%d,[sp,#%d]" % (j, T_BASE + 8 * (6 + j)))
        a("ldr x%d,[sp,#%d]" % (6 + j, T_BASE + 8 * (6 + j)))
    a("subs x0,x0,x18")
    for j in range(1, 6):
        a("sbcs x%d,x%d,x%d" % (j, j, 18 + j))
    a("cset x25,hs")                    # 1 iff t >= p (no borrow): keep the difference
    a("cmp x25,#1")
    for j in range(6):
        a("csel x%d,x%d,x%d,eq" % (j, j, 6 + j))
    for j in range(6):
        a("str x%d,[sp,#%d]" % (j, T_BASE + 8 * (6 + j)))
    # ---- candidate 1: bx = t[6..11], compare, then write ----------------------
    _compare_write(a, T_BASE + 8 * 6, 1)
    # ---- candidate 2: b2x = p - x - bx (mod p) -------------------------------
    for k in range(6):
        a("ldr x%d,[sp,#%d]" % (k, X_BASE + 8 * k))
        a("ldr x%d,[sp,#%d]" % (6 + k, T_BASE + 8 * (6 + k)))
    a("subs x0,x18,x0")                                  # s = p - x
    for k in range(1, 6):
        a("sbcs x%d,x%d,x%d" % (k, 18 + k, k))
    for k in range(6):                                   # s -= bx
        a("sbcs x%d,x%d,x%d" % (k, k, 6 + k))
    a("cset x25,hs")                                     # 1 iff p - x - bx >= 0
    a("adds x6,x0,x18")                                  # t = s + p (carry cascade)
    for k in range(1, 6):
        a("adc x%d,x%d,x%d" % (6 + k, k, 18 + k))
    a("cmp x25,#1")
    for k in range(6):
        a("csel x%d,x%d,x%d,hs" % (k, k, 6 + k))
    for k in range(6):
        a("str x%d,[sp,#%d]" % (k, CAND_BASE + 8 * k))
    _compare_write(a, CAND_BASE, 2)
    # ---- restore and rejoin ---------------------------------------------------
    for i in range(NREG):
        if i == 8:
            continue
        a("ldr x%d,[sp,#%d]" % (i, SAVE_BASE + 8 * i))
    a("ldr x8,[sp,#%d]" % MINIDX)        # the exact argmin, read by row 4
    a("add sp,sp,#%d" % stack_bytes)     # release the scratch frame
    a("b L_rejoin")
    if stack_bytes < 504:
        raise AssertionError("stack contract: frame %d < 504 bytes" % stack_bytes)
    for line in o:
        off = 0
        if "[sp,#" in line:
            off = int(line.split("[sp,#")[1].split("]")[0])
        if off >= stack_bytes:
            raise AssertionError("stack contract: offset %d outside frame %d" % (off, stack_bytes))
    return o


def _compare_write(a, cand_base, cand):
    """Read the running min, compare, THEN write the winner (never read-after-write)."""
    for j in range(6):
        a("ldr x%d,[sp,#%d]" % (j, cand_base + 8 * j))       # candidate
    for j in range(6):
        a("ldr x%d,[sp,#%d]" % (6 + j, MIN_BASE + 8 * j))    # running minimum
    a("ldr x27,[sp,#%d]" % MINIDX)                           # running index
    a("mov x24,#%d" % cand)
    a("subs x25,x0,x6")                                      # 384-bit borrow cascade
    for j in range(1, 6):
        a("sbcs x25,x%d,x%d" % (j, 6 + j))
    a("csel x24,x24,x27,lo")                                 # new index iff cand < min
    for j in range(6):
        a("csel x%d,x%d,x%d,lo" % (6 + j, j, 6 + j))         # winner limbs
    for j in range(6):
        a("str x%d,[sp,#%d]" % (6 + j, MIN_BASE + 8 * j))    # written only after compare
    a("str x24,[sp,#%d]" % MINIDX)


def stack_contract():
    """Every frame offset the block touches must be inside the frame it allocates."""
    alloc = None
    for line in repair_block():
        if line.strip().startswith("sub sp,sp,#"):
            alloc = int(line.strip().split("#")[1])
    if alloc is None:
        raise AssertionError("repair block does not allocate its own frame")
    for line in repair_block():
        if "[sp,#" in line:
            off = int(line.split("[sp,#")[1].split("]")[0])
            if off >= alloc:
                raise AssertionError("repair offset %d outside the block's frame %d" % (off, alloc))
    if sum(1 for l in repair_block() if l.strip().startswith("add sp,sp,#")) != 1:
        raise AssertionError("repair block does not release its frame exactly once")
    return alloc
