"""band3_a64vm.py — Pure-stdlib baseline scalar ARMv8.0-A parser, budget verifier, and live VM for KERNELS.md §5."""
import math, re, struct
from fractions import Fraction as Frac
from pathlib import Path
import model_bn254_v6 as M

ALLOWED_A64_OPS = frozenset({
    "mov", "add", "sub", "subs", "neg", "cneg", "mul", "umulh", "madd", "msub",
    "lsl", "lsr", "asr", "and", "ubfx", "bfi", "cmp", "tst", "csel",
    "ldr", "ldrh", "ldp", "str", "strh", "stp", "cbz", "cbnz", "tbnz",
    "b", "bl", "b.eq", "b.ne", "b.lt", "b.le", "b.gt", "b.hs", "ret",
})
_MASK64 = (1 << 64) - 1
_SIGN64 = 1 << 63

def _s64(v):
    v &= _MASK64
    return v - (1 << 64) if (v & _SIGN64) else v

def _reg_id(name):
    n = name.strip().lower()
    if n in ("xzr", "wzr"): return (32, False)
    if n == "sp": return (31, False)
    return (int(n[1:]), n[0] == "w")

class A64VM:
    """Fast pre-compiled baseline scalar ARMv8.0-A integer VM executing KERNELS.md §5 live."""
    def __init__(self, insns, labels, mem_bytes):
        self.labels = labels
        self.mem = bytearray(mem_bytes)
        self.r = [0] * 33  # 0..30=x0..x30, 31=sp, 32=xzr
        self.N = self.Z = self.C = self.V = 0
        self.exec_cnt = 0
        self.hit_labels = set()
        self.pc_lbl = {}
        for k, v in labels.items():
            self.pc_lbl.setdefault(v, []).append(k)
        self.prog = self._compile(insns, labels)

    def _parse_op2(self, parts):
        if len(parts) == 1:
            t = parts[0].strip()
            return ("imm", int(t[1:], 0) & _MASK64) if t.startswith("#") else ("reg", _reg_id(t))
        return ("sh", _reg_id(parts[0]), parts[1].strip().lower().split()[0], int(parts[1].split("#")[1], 0))

    def _compile(self, insns, labels):
        prog = []
        for op, argstr in insns:
            if op == "ret":
                prog.append(("ret",))
            elif op in ("b", "bl"):
                prog.append((op, labels[argstr.strip()]))
            elif op.startswith("b."):
                prog.append(("bcc", op[2:], labels[argstr.strip()]))
            elif "[" in argstr:
                i0 = argstr.index("[")
                head = [_reg_id(x) for x in argstr[:i0].rstrip(", ").split(",") if x.strip()]
                ms = argstr[i0:].strip()
                pre = False; post = 0
                if ms.endswith("!"): pre = True; inner = ms[1:-2]
                elif "]," in ms: p0, p1 = ms.split("],", 1); inner = p0[1:]; post = int(p1.strip().lstrip("#"), 0)
                else: inner = ms[1:-1]
                ps = [x.strip() for x in inner.split(",")]
                br = _reg_id(ps[0])
                if len(ps) == 1: mode = ("imm", 0)
                elif len(ps) == 2: mode = ("imm", int(ps[1][1:], 0)) if ps[1].startswith("#") else ("reg", _reg_id(ps[1]), 0)
                else: mode = ("reg", _reg_id(ps[1]), int(ps[2].split("#")[1], 0))
                prog.append(("mem", op, head, br, mode, pre, post))
            else:
                a = [x.strip() for x in argstr.split(",")]
                if op == "mov":
                    prog.append(("mov", _reg_id(a[0]), ("imm", int(a[1][1:], 0) & _MASK64) if a[1].startswith("#") else ("reg", _reg_id(a[1]))))
                elif op == "neg": prog.append(("neg", _reg_id(a[0]), _reg_id(a[1])))
                elif op == "cneg": prog.append(("cneg", _reg_id(a[0]), _reg_id(a[1]), a[2]))
                elif op == "csel": prog.append(("csel", _reg_id(a[0]), _reg_id(a[1]), _reg_id(a[2]), a[3]))
                elif op in ("add", "sub", "subs", "and"): prog.append((op, _reg_id(a[0]), _reg_id(a[1]), self._parse_op2(a[2:])))
                elif op in ("cmp", "tst"): prog.append((op, _reg_id(a[0]), self._parse_op2(a[1:])))
                elif op in ("lsl", "lsr", "asr"): prog.append((op, _reg_id(a[0]), _reg_id(a[1]), int(a[2][1:], 0)))
                elif op in ("ubfx", "bfi"): prog.append((op, _reg_id(a[0]), _reg_id(a[1]), int(a[2][1:], 0), int(a[3][1:], 0)))
                elif op in ("mul", "umulh"): prog.append((op, _reg_id(a[0]), _reg_id(a[1]), _reg_id(a[2])))
                elif op in ("madd", "msub"): prog.append((op, _reg_id(a[0]), _reg_id(a[1]), _reg_id(a[2]), _reg_id(a[3])))
                elif op in ("cbz", "cbnz"): prog.append((op, _reg_id(a[0]), labels[a[1]]))
                elif op == "tbnz": prog.append(("tbnz", _reg_id(a[0]), int(a[1][1:], 0), labels[a[2]]))
                else: raise ValueError(f"Unsupported op: {op} {argstr}")
        return prog

    def load64(self, a): return struct.unpack_from("<Q", self.mem, a)[0]
    def store64(self, a, v): struct.pack_into("<Q", self.mem, a, v & _MASK64)
    def load16(self, a): return struct.unpack_from("<H", self.mem, a)[0]
    def store16(self, a, v): struct.pack_into("<H", self.mem, a, v & 0xFFFF)

    def _cond(self, c):
        if c == "eq": return self.Z == 1
        if c == "ne": return self.Z == 0
        if c == "lt": return self.N != self.V
        if c == "ge": return self.N == self.V
        if c == "gt": return self.Z == 0 and self.N == self.V
        if c == "le": return self.Z == 1 or self.N != self.V
        if c == "hs": return self.C == 1
        return self.C == 0

    def call(self, sym):
        SENT = _MASK64
        r = self.r; mem = self.mem; prog = self.prog; pc_lbl = self.pc_lbl; hl = self.hit_labels
        u64 = struct.unpack_from; p64 = struct.pack_into
        r[30] = SENT; pc = self.labels[sym]
        def gr(ri):
            v = r[ri[0]]
            return (v & 0xFFFFFFFF) if ri[1] else (v & _MASK64)
        def wr(ri, v):
            if ri[0] < 32: r[ri[0]] = (v & 0xFFFFFFFF) if ri[1] else (v & _MASK64)
        def e2(x):
            if x[0] == "imm": return x[1]
            if x[0] == "reg": return gr(x[1])
            v = gr(x[1]); sh = x[3]
            if x[2] == "lsl": return (v << sh) & _MASK64
            if x[2] == "lsr": return v >> sh
            return (_s64(v) >> sh) & _MASK64
        def fsub(a, b):
            res = (a - b) & _MASK64
            self.N = (res >> 63) & 1; self.Z = 1 if res == 0 else 0; self.C = 1 if a >= b else 0
            sa = a - (1 << 64) if (a & _SIGN64) else a
            sb = b - (1 << 64) if (b & _SIGN64) else b
            sr = res - (1 << 64) if (res & _SIGN64) else res
            self.V = 1 if ((sa >= 0 and sb < 0 and sr < 0) or (sa < 0 and sb >= 0 and sr >= 0)) else 0
            return res
        while pc != SENT:
            if pc in pc_lbl:
                for lb in pc_lbl[pc]: hl.add(lb)
            ins = prog[pc]; k = ins[0]; self.exec_cnt += 1
            if k == "ret": pc = r[30]; continue
            if k == "b": pc = ins[1]; continue
            if k == "bl": r[30] = pc + 1; pc = ins[1]; continue
            if k == "bcc": pc = ins[2] if self._cond(ins[1]) else pc + 1; continue
            if k == "mem":
                _, op, head, br, mode, pre, post = ins
                base = gr(br)
                off = mode[1] if mode[0] == "imm" else (_s64(gr(mode[1])) << mode[2])
                addr = (base + off) & _MASK64
                if pre: wr(br, addr)
                if op == "ldr": wr(head[0], u64("<Q", mem, addr)[0])
                elif op == "ldrh": wr(head[0], u64("<H", mem, addr)[0])
                elif op == "ldp": wr(head[0], u64("<Q", mem, addr)[0]); wr(head[1], u64("<Q", mem, addr + 8)[0])
                elif op == "str": p64("<Q", mem, addr, gr(head[0]))
                elif op == "strh": p64("<H", mem, addr, gr(head[0]) & 0xFFFF)
                else: p64("<Q", mem, addr, gr(head[0])); p64("<Q", mem, addr + 8, gr(head[1]))
                if post: wr(br, gr(br) + post)
            elif k == "mov": wr(ins[1], e2(ins[2]))
            elif k == "neg": wr(ins[1], (-_s64(gr(ins[2]))) & _MASK64)
            elif k == "cneg": v = gr(ins[2]); wr(ins[1], ((-_s64(v)) & _MASK64) if self._cond(ins[3]) else v)
            elif k == "csel": wr(ins[1], gr(ins[2]) if self._cond(ins[4]) else gr(ins[3]))
            elif k == "add": wr(ins[1], gr(ins[2]) + e2(ins[3]))
            elif k == "sub": wr(ins[1], gr(ins[2]) - e2(ins[3]))
            elif k == "subs": wr(ins[1], fsub(gr(ins[2]), e2(ins[3])))
            elif k == "and": wr(ins[1], gr(ins[2]) & e2(ins[3]))
            elif k == "cmp": fsub(gr(ins[1]), e2(ins[2]))
            elif k == "tst":
                res = gr(ins[1]) & e2(ins[2])
                self.N = (res >> 63) & 1; self.Z = 1 if res == 0 else 0; self.C = 0; self.V = 0
            elif k == "lsl": wr(ins[1], (gr(ins[2]) << ins[3]) & _MASK64)
            elif k == "lsr": wr(ins[1], gr(ins[2]) >> ins[3])
            elif k == "asr": wr(ins[1], (_s64(gr(ins[2])) >> ins[3]) & _MASK64)
            elif k == "ubfx": wr(ins[1], (gr(ins[2]) >> ins[3]) & ((1 << ins[4]) - 1))
            elif k == "bfi":
                mk = ((1 << ins[4]) - 1) << ins[3]
                wr(ins[1], (gr(ins[1]) & ~mk) | ((gr(ins[2]) << ins[3]) & mk))
            elif k == "mul": wr(ins[1], (_s64(gr(ins[2])) * _s64(gr(ins[3]))) & _MASK64)
            elif k == "umulh": wr(ins[1], ((gr(ins[2]) * gr(ins[3])) >> 64) & _MASK64)
            elif k == "madd": wr(ins[1], (_s64(gr(ins[2])) * _s64(gr(ins[3])) + _s64(gr(ins[4]))) & _MASK64)
            elif k == "msub": wr(ins[1], (_s64(gr(ins[4])) - _s64(gr(ins[2])) * _s64(gr(ins[3]))) & _MASK64)
            elif k == "cbz": pc = ins[2] if gr(ins[1]) == 0 else pc + 1; continue
            elif k == "cbnz": pc = ins[2] if gr(ins[1]) != 0 else pc + 1; continue
            elif k == "tbnz": pc = ins[3] if ((gr(ins[1]) >> ins[2]) & 1) else pc + 1; continue
            pc += 1

def verify_band3_asm_schedules(kernels_md_path):
    """Parse KERNELS.md §5 fenced asm block live, enforce baseline scalar ARMv8.0-A allowlist + zero x18, and count all levels k=0..11."""
    txt = Path(kernels_md_path).read_text()
    m = re.search(r"```asm\n(.*?)\n```", txt, flags=re.S)
    if not m:
        return dict(ok=False, err="missing asm block"), [], {}
    blocks, labels, insns = {}, {}, []
    cur_label = None; globals_seen = set(); bad_ops = []; x18_uses = 0
    for raw_ln in m.group(1).splitlines():
        ln = raw_ln.split("//", 1)[0].strip()
        if not ln: continue
        if ln.startswith(".global "):
            globals_seen.add(ln.split()[1]); continue
        if ln.startswith(".") and not ln.endswith(":"): continue
        if ln.endswith(":"):
            cur_label = ln[:-1]; labels[cur_label] = len(insns); blocks.setdefault(cur_label, 0); continue
        toks = ln.split(None, 1)
        op = toks[0].lower(); ops_str = toks[1] if len(toks) > 1 else ""
        if op not in ALLOWED_A64_OPS: bad_ops.append(op)
        if re.search(r"\b[xw]18\b", ops_str, flags=re.I): x18_uses += 1
        insns.append((op, ops_str))
        blocks[cur_label] = blocks.get(cur_label, 0) + 1
    req_syms = {"band3_k1_node", "band3_k1_to_k0_scan", "band3_k_ge_2_node", "band3_apply16_unroll"}
    line_ab_ops = blocks.get(".Lk1_loop", 0)
    cand_ops = blocks.get(".Lloop", 0)
    line_f_ops = blocks.get(".Lk1_next", 0) - 1
    line_ops = line_ab_ops + line_f_ops
    hit_ops = blocks.get(".Lhit", 0) + blocks.get(".Ltable_b_loaded", 0) + blocks.get(".Lhit_done", 0) - 1
    alloc_ops = blocks.get(".Lalloc", 0)
    k1_node_ops = blocks.get("band3_k1_node", 0) + blocks.get(".Lk1_nz_ok", 0) + 1 + blocks.get(".Lk1_ret", 0)
    k2_loop_ops = blocks.get(".Lk2_child_loop", 0) + blocks.get(".Lk2_do_k1", 0)
    k2_node_ops = (blocks.get("band3_k_ge_2_node", 0) + blocks.get(".Lk_b", 0) + blocks.get(".Lk_done_int", 0)
                   + blocks.get(".Lk_loop_chk", 0) + blocks.get(".Lk2_init", 0)
                   + blocks.get(".Lk2_restore", 0) + blocks.get(".Lk_restore_ab", 0) + blocks.get(".Lk_ret", 0))
    k3_node_ops = (blocks.get("band3_k_ge_2_node", 0) + blocks.get(".Lk_b", 0) + blocks.get(".Lk_done_int", 0)
                   + blocks.get(".Lk_loop_chk", 0) + blocks.get(".Lk3_init_c", 0) + 4 * blocks.get(".Lk_init_c", 0)
                   + blocks.get(".Lk3_rst_prep", 0) + 4 * blocks.get(".Lk_rst_c", 0)
                   + blocks.get(".Lk_restore_ab", 0) + blocks.get(".Lk_ret", 0))
    k3_child_fixed = blocks.get(".Lk3_child_loop", 0) + blocks.get(".Lk3_chk_b", 0) + blocks.get(".Lk3_do_rec", 0) + blocks.get(".Lk3_call", 0) + blocks.get(".Lk3_step", 0) + 3
    k3_child_per_i = blocks.get(".Lk3_upd_c", 0) - 3
    k3_child_ops = k3_child_fixed + 4.0 * k3_child_per_i
    apply16_ops = blocks.get(".Lapp16", 0) - 1
    hit_eff_ir = hit_ops + 0.45 + alloc_ops / 511.0
    apply_eff_ir = apply16_ops / 16.0 + (15.0 * 1.35 + 12.0) / 511.0
    derived_ir = (
        2.58 * cand_ops + hit_eff_ir + 0.45 * line_ops
        + 0.068 * (k1_node_ops + k2_loop_ops)
        + 0.009 * k2_node_ops + 0.0013 * k3_node_ops + 0.010 * k3_child_ops + apply_eff_ir
    )
    ok = (
        req_syms.issubset(globals_seen) and not bad_ops and x18_uses == 0
        and line_ops == 40 and cand_ops == 5 and hit_ops == 31 and alloc_ops == 9
        and k1_node_ops == 111 and k2_loop_ops == 56 and k2_node_ops == 194
        and k3_node_ops == 237 and abs(k3_child_ops - 137.0) < 1e-6 and apply16_ops == 58
        and derived_ir <= M.C_LBS
    )
    return dict(
        ok=ok, total_insns=len(insns), bad_ops=bad_ops, x18_uses=x18_uses,
        line_ops=line_ops, cand_ops=cand_ops, hit_ops=hit_ops, hit_eff_ir=round(hit_eff_ir, 3),
        alloc_ops=alloc_ops, k1_node_ops=k1_node_ops, k2_loop_ops=k2_loop_ops, k2_node_ops=k2_node_ops,
        k3_node_ops=k3_node_ops, k3_child_ops=round(k3_child_ops, 1),
        apply16_ops=apply16_ops, apply_eff_ir=round(apply_eff_ir, 3), derived_ir=round(derived_ir, 2),
    ), insns, labels

def run_orbit_in_a64vm(insns, labels, canon_a, table_b, SQRT16, num_sbs, Sa_box, idx_a, sidx_b5,
                       S, S2, mv6, B, Da_f, Db_f, Dw_f, T_a, T_b, Tw_i,
                       Da_sc, Db_sc, Dw_sc, inva, invb, invw, Ua_x, Ub_x, Uw_x, ETA, logp_q):
    """Execute KERNELS.md §5 band3_k_ge_2_node LIVE on one real BN254 C6 orbit (with forced .Lalloc)."""
    S24 = 1 << 24; S30 = 1 << 30; S48 = 1 << 48
    def _q48p(dfk, ti):
        dn = dfk / Frac(ti)
        return ((2 * dn * S48).numerator // (2 * dn * S48).denominator,
                (dn * S30).numerator // (dn * S30).denominator,
                int(math.ceil(2.0 * float(dn) * (ETA / S24) * S48)) + (1 << 26))
    ddqa1, da1_q30, msha1 = _q48p(Da_f[1], T_a)
    ddqb1, db1_q30, mshb1 = _q48p(Db_f[1], T_b)
    ddqw2, dw2_q30, mshw2 = _q48p(Dw_f[2], Tw_i)
    ddqa2, da2_q30, msha2 = _q48p(Da_f[2], T_a)
    ddqb2, db2_q30, mshb2 = _q48p(Db_f[2], T_b)
    mem_sz = 8 * len(canon_a) + 16 * len(table_b) + (num_sbs + 256) * 4096 + 524288
    vm = A64VM(insns, labels, mem_sz)
    vptr = [4096]
    def _va(nb, al=4096):
        vptr[0] = (vptr[0] + al - 1) & ~(al - 1); out_a = vptr[0]; vptr[0] += nb; return out_a
    ca_ad = _va(8 * len(canon_a), 64)
    for i_a, v_a in enumerate(canon_a):
        if v_a: vm.store64(ca_ad + 8 * i_a, v_a)
    tb_ad = _va(16 * len(table_b), 64)
    for i_b, (w0_b, w1_b) in enumerate(table_b):
        vm.store64(tb_ad + 16 * i_b, w0_b); vm.store64(tb_ad + 16 * i_b + 8, w1_b)
    sq_ad = _va(2 * len(SQRT16), 64)
    for i_s, v_s in enumerate(SQRT16): vm.store16(sq_ad + 2 * i_s, v_s)
    ctx_ad = _va(256, 64); lvl_ad = _va(12 * 256, 64); cs_ad = _va(3 * 96, 64)
    bl_ad = _va(3 * 64, 64); sb_ad = _va(8 * num_sbs, 64); fp_ad = _va(8, 64)
    pg_base = _va((num_sbs + 192) * 4096, 4096)
    init_cur = []
    for sbi in range(num_sbs):
        pg0 = pg_base + sbi * 4096; vm.store64(pg0, 0)
        cur0 = pg0 + 4080; vm.store64(sb_ad + 8 * sbi, cur0); init_cur.append(cur0)
    fp0 = pg_base + num_sbs * 4096; vm.store64(fp_ad, fp0)
    da_st = [8 * sum(B[k][ci] * Sa_box[ci] for ci in range(6)) for k in range(12)]
    sb5_st = [[0] * 12 for _ in range(3)]; db0_st = [[0] * 12 for _ in range(3)]
    for k in range(12):
        bk = tuple(B[k][6:]); bks = (bk, mv6(S, bk), mv6(S2, bk))
        for kc in range(3): sb5_st[kc][k] = sidx_b5(bks[kc]); db0_st[kc][k] = bks[kc][0]
    for kc in range(3):
        b0_a = bl_ad + 64 * kc
        vm.store64(b0_a + 16, sb5_st[kc][1]); vm.store64(b0_a + 24, db0_st[kc][1])
        vm.store64(b0_a + 32, sb5_st[kc][0]); vm.store64(b0_a + 40, db0_st[kc][0])
    for k in range(12):
        lb = lvl_ad + 256 * k
        vm.store64(lb + 0, invw[k]); vm.store64(lb + 8, inva[k]); vm.store64(lb + 16, invb[k])
        vm.store64(lb + 24, Dw_sc[k]); vm.store64(lb + 32, Da_sc[k]); vm.store64(lb + 40, Db_sc[k])
        vm.store64(lb + 48, da_st[k])
        for kc in range(3):
            vm.store64(lb + 56 + 16 * kc, sb5_st[kc][k]); vm.store64(lb + 64 + 16 * kc, db0_st[kc][k])
        uwc, uac, ubc = _va(8 * max(1, k), 8), _va(8 * max(1, k), 8), _va(8 * max(1, k), 8)
        for i_u in range(k):
            vm.store64(uwc + 8 * i_u, Uw_x[i_u][k]); vm.store64(uac + 8 * i_u, Ua_x[i_u][k]); vm.store64(ubc + 8 * i_u, Ub_x[i_u][k])
        vm.store64(lb + 104, uwc); vm.store64(lb + 112, uac); vm.store64(lb + 120, ubc)
        if k == 2:
            for fi, fv in enumerate((ddqw2, ddqa2, ddqb2, dw2_q30, mshw2, da2_q30, msha2, db2_q30, mshb2,
                                     Uw_x[1][2], Ua_x[0][2], Ua_x[1][2], Ub_x[0][2], Ub_x[1][2])):
                vm.store64(lb + 128 + 8 * fi, fv)
    for ci, cv in enumerate((ddqa1, ddqb1, Ua_x[0][1], Ub_x[0][1], inva[0], invb[0], fp_ad,
                             da_st[0], da_st[1], Da_sc[1], da1_q30, msha1, Db_sc[1], db1_q30, mshb1,
                             ETA, ETA + 1, lvl_ad, cs_ad, invw[1], inva[1], invb[1])):
        vm.store64(ctx_ad + 8 * ci, cv)
    vm.r[31] = mem_sz - 64; vm.r[19] = ctx_ad; vm.r[24] = ca_ad + 8 * idx_a((0,) * 6)
    vm.r[25] = sb_ad; vm.r[26] = tb_ad; vm.r[27] = bl_ad; vm.r[28] = sq_ad
    vm.r[0] = vm.r[1] = vm.r[2] = vm.r[3] = 0; vm.r[4] = 11
    vm.call("band3_k_ge_2_node")
    vm_sc, sample = {}, []
    for sbi in range(num_sbs):
        cur = vm.load64(sb_ad + 8 * sbi); pg_c = cur & ~4095; end_c = cur
        while pg_c != 0:
            st_c = init_cur[sbi] if pg_c == (init_cur[sbi] & ~4095) else (pg_c + 8)
            for ap in range(st_c, end_c, 8):
                ei = vm.load64(ap)
                vm_sc[(sbi, ei)] = vm_sc.get((sbi, ei), 0) + logp_q
                if len(sample) < 32: sample.append(ei & 0x3FF)
            pg_c = vm.load64(pg_c); end_c = pg_c + 4096
    return vm_sc, (vm.load64(fp_ad) - fp0) // 4096, vm.exec_cnt, vm.hit_labels, sample

def verify_aux_kernels_in_a64vm(insns, labels, sample32):
    """Execute band3_apply16_unroll (32 emitted entries), explicit 3-lane .Lk2, and .Lsidx0 live."""
    vm_app = A64VM(insns, labels, 65536)
    for i_e, e_v in enumerate(sample32[:32]): vm_app.store64(4096 + 8 * i_e, e_v)
    vm_app.r[31] = 65000; vm_app.r[0] = 4096; vm_app.r[1] = 32; vm_app.r[2] = 8192; vm_app.r[3] = 147
    vm_app.call("band3_apply16_unroll")
    exp_app = {}
    for e_v in sample32[:32]: exp_app[e_v] = exp_app.get(e_v, 0) + 147
    app_errs = sum(1 for i_e in range(1024) if vm_app.load16(8192 + 2 * i_e) != exp_app.get(i_e, 0))
    vm_k2 = A64VM(insns, labels, 131072)
    vm_k2.store64(4096, (1 << 63) | 777)
    for kc in range(3):
        vm_k2.store64(4224 + 16 * (kc + 1), kc)
        vm_k2.store64(4224 + 16 * (kc + 1) + 8, ((-2 & 0x3FFF) << 50) | (10 << 39) | (kc + 1))
        vm_k2.store64(4608 + 64 * kc, kc + 1); vm_k2.store64(4608 + 64 * kc + 8, kc)
        vm_k2.store64(4864 + 8 * kc, 8192 + 4096 * kc + 8)
    vm_k2.store64(5120 + 48, 32768)
    vm_k2.r[31] = 130000; vm_k2.r[19] = 5120; vm_k2.r[25] = 4864; vm_k2.r[26] = 4224; vm_k2.r[27] = 4608
    vm_k2.r[20] = 0; vm_k2.r[23] = -1; vm_k2.r[8] = 1; vm_k2.r[12] = 1; vm_k2.r[10] = 4096; vm_k2.r[6] = 8
    vm_k2.call(".Lorig_enter")
    k2_ok = ([vm_k2.load64(8192 + 4096 * kc + 8) for kc in range(3)] == [((kc + 1) << 55) | ((kc + 2) << 45) | 777 for kc in range(3)] and ".Lk2" in vm_k2.hit_labels)
    # Explicit live .Lsidx0 test: sidx == 0, b0 = -3 -> b0_s = +3 (ch_id = 19, sb_id = 0, b0_min = 1, L_b = 10 -> b0_rel = 2)
    vm_s0 = A64VM(insns, labels, 131072)
    vm_s0.store64(4096, 555)
    vm_s0.store64(4224, 0)
    vm_s0.store64(4224 + 8, ((1 & 0x3FFF) << 50) | (10 << 39) | 19)
    vm_s0.store64(4608 + 0, 0); vm_s0.store64(4608 + 8, (-3) & _MASK64)
    vm_s0.store64(4864, 8192 + 8); vm_s0.store64(5120 + 48, 32768)
    vm_s0.r[31] = 130000; vm_s0.r[19] = 5120; vm_s0.r[25] = 4864; vm_s0.r[26] = 4224; vm_s0.r[27] = 4608
    vm_s0.r[20] = 0; vm_s0.r[23] = -1; vm_s0.r[8] = 1; vm_s0.r[12] = 1; vm_s0.r[10] = 4096; vm_s0.r[6] = 8
    vm_s0.call(".Lorig_enter")
    sidx0_ok = (vm_s0.load64(8192 + 8) == ((19 << 55) | (2 << 45) | 555) and ".Lsidx0" in vm_s0.hit_labels)
    return app_errs, k2_ok, sidx0_ok
