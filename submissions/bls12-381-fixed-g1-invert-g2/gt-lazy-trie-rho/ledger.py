"""Part D of experiment.py (proof.md section 5): the exact work ledger and the success chain.
Loaded by experiment.py with runpy and its globals (split out for the 64 KiB per-file limit)."""

# Ledger (part D)
CLAIM_TXT = "124.209"
P_MODEL_TXT = "0.4014"


def ledger(cyc):
    """Exact recomputation of the work ledger of proof.md section 5 (instruction counts / C, C = 20,000)."""
    from decimal import Decimal, getcontext, localcontext
    getcontext().prec = 60
    ln2 = Decimal(2).ln()
    lg = lambda v: float(Decimal(v).ln() / ln2)
    D = lambda fr: Decimal(fr.numerator) / Decimal(fr.denominator)
    C = C_UNIT
    Rt, Lt = R_PROD, L_PROD
    # step prices (instructions) and the per-step average over a complete block
    c = Fraction(BLOCK_I, Lt * C)
    cmax = Fraction(STEP_FLUSH_I, C)
    assert c == Fraction(41821, 120000)
    assert all(Fraction(k * STEP_LAZY_I, 1) <= k * Lt * c * C for k in range(1, Lt))
    assert STEP_ROOT_I < STEP_LAZY_I < c * C < STEP_FLUSH_I
    # partial blocks: k < L steps cost at most STEP_ROOT_I + (k-1) STEP_LAZY_I < k c C
    assert all(STEP_ROOT_I + (k - 1) * STEP_LAZY_I < k * c * C for k in range(1, Lt))
    # collision law
    N = Fraction(r - 1, 12)
    Pm = Decimal(P_MODEL_TXT)
    sumq2 = Fraction(1, Rt)
    Np = N / (1 - sumq2)
    Npd = D(Np)
    S = (2 * Npd * (1 / (1 - Pm)).ln()).sqrt()
    # fruitless closed walks (part F): eps in step equivalents, 1.25 x the H2 model
    esc_steps = Fraction(ESC_INSTR, C) / c
    P2 = Fraction(1, 12 * Rt)
    eps2 = P2 * (2 + esc_steps)
    eps_long = sum(Fraction(pl).limit_denominator(10 ** 30) * (64 + l + l * cmax / c + esc_steps)
                   for l, pl in zip(range(3, 65), cyc["P_l_by_length_3_to_64"]))
    EPS_MARGIN = Fraction(5, 4)
    assert cyc["P_2_equals_1_over_12R"]
    # trails, restarts, detection window: the stuck counter counts every step, escape traversals included, and a
    # trail restarts at A_STUCK = a / theta, a = 120 (proof.md 3.5, detection lemma)
    WALKS, RESTART_CAP, A_DET = 2 ** 20, 2 ** 20, 120
    STUCK = A_DET * Decimal(2) ** THETA_BITS
    assert STUCK == 0x1E000000 * Decimal(2) ** 64
    E = (WALKS + RESTART_CAP) * STUCK
    # closures longer than 64 (part F rate, per first-visit iteration) are caught by the doubling checkpoint within
    # 3 (mu + l) steps and traversed once (l steps), l <= mu + l <= STUCK: at most 4 STUCK steps at cmax, plus the escape
    LONG_RATE = Decimal(2) ** -108
    assert cyc["rate_useless_walks_longer_than_64_log2"] < -108
    eps_brent = Fraction(1, 2 ** 108) * (4 * Fraction(int(STUCK)) * cmax / c + esc_steps)
    eps = D(EPS_MARGIN * (eps2 + eps_long + eps_brent))
    # Hard iteration maximum: a block of k <= 6 iterations with the escape or report and fresh start ending it is
    # charged at least k c C (asserted below), so n <= B_walk / (c C) + 6 M (last partial block per walk), B_walk < W_cap C.
    for k in range(1, Lt + 1):
        blk = STEP_ROOT_I + (k - 1) * STEP_LAZY_I + (FLUSH if k == Lt else 0)
        assert blk + ESC_INSTR >= k * c * C and blk + 4000 + int(2 ** 10) * C >= k * c * C
        assert k < Lt or blk == Lt * c * C
    N_MAX = Decimal(2) ** Decimal(CLAIM_TXT) / D(c) + 6 * 2 ** 20
    THETA = Decimal(2) ** (-THETA_BITS)
    W_BAR = D(eps2 + eps_long + eps_brent)               # H2-model waste per first-visit iteration (unscaled)
    B_MAX = D(64 + 64 + 64 * cmax / c + esc_steps)
    assert B_MAX < 312
    # H1 (proof.md 3.4): a random mapping on the escape-reduced class walk; an image landing in the trail's current
    # segment is a structural return (a fruitless closure), charged in eps and escaped, never counted as a collision.
    # H1+ (disclosed, not proven): the realized waste of the run is at most T_TOL times the H2 rate W_BAR per useful
    # iteration. Under it the meter still buys S_TOL = S (1 + eps) / (1 + T_TOL W_BAR) useful iterations, and the loss
    # is the drop of the collision law from S to S_TOL (computed, not flat).
    T_TOL = 2
    # part B's toy metric uses these rounded-down constants and the same model point
    assert Decimal(str(TOY_EPS)) <= eps and Decimal(str(TOY_WBAR)) <= W_BAR and Decimal(str(TOY_PM)) == Pm
    assert eps - Decimal(str(TOY_EPS)) < Decimal("1e-6") and W_BAR - Decimal(str(TOY_WBAR)) < Decimal("1e-6")
    assert 1 + W_BAR < 1 + eps < 1 + T_TOL * W_BAR
    X_COL = S * S / Npd
    SLOPE = X_COL * (-(X_COL / 2)).exp()      # dP/dS * S, below 2/e for every S, used for the restart-cap fraction
    assert SLOPE < 1
    # Detection lemma (proof.md 3.5): outside a self-collision inside the trail, a trail of U steps visits at least U/4
    # distinct classes, so a trail is stuck (no DP in STUCK steps) with probability at most (1 - theta)^(a/4 theta)
    # <= e^-30, plus the self-collision probability STUCK^2 / (2 N') (pairwise H1).
    P_STUCK = (-Decimal(A_DET) / 4).exp() + STUCK * STUCK / (2 * Npd)
    assert P_STUCK < Decimal("9.36e-14")
    # Markov on every cap. Reports: one per trail end at a DP; a trail ending on a fresh class does so with probability
    # theta (coupling), a merge into an earlier trail adds at most one report, merges are pairs of first visits
    MU_REPORTS = N_MAX * THETA + N_MAX * N_MAX / (2 * Npd) + WALKS
    EXP_RESTARTS = P_STUCK * (MU_REPORTS + WALKS) / (1 - P_STUCK)
    assert lg(EXP_RESTARTS) < 20 - 20
    RESTART_STEPS = Decimal(RESTART_CAP) * STUCK
    RESTART_FRAC = RESTART_STEPS / S              # conservative: dP/dS * S = SLOPE < 1
    # reports and their log
    REPORT_CAP, SLOT_BYTES = 2 ** 55, 80
    assert lg(MU_REPORTS) < 55 - 15
    REPORT_TAIL = MU_REPORTS / REPORT_CAP
    RESTART_TAIL = EXP_RESTARTS / RESTART_CAP
    # false fingerprint pairs: R <= REPORT_CAP, so E[R^2] <= REPORT_CAP E[R]; the matcher stops at the first true pair,
    # so the budget of 2^4 pairs fails only with at least 2^4 - 1 false pairs
    PAIR_BUDGET = 2 ** 4
    FALSE_PAIRS = Decimal(REPORT_CAP) * MU_REPORTS / 2 / Decimal(2) ** 112
    PAIR_TAIL = FALSE_PAIRS / (PAIR_BUDGET - 1)
    assert REPORT_TAIL < Decimal(2) ** -15 and RESTART_TAIL < Decimal(2) ** -23 and PAIR_TAIL < Decimal(2) ** -21
    report_instr = 4000
    # sampler
    SAMPLE_TABLE_DRAWS = 2 * Rt
    SAMPLES_MAX = Decimal(SAMPLE_TABLE_DRAWS + 2 * (WALKS + REPORT_CAP + RESTART_CAP))
    assert lg(SAMPLES_MAX) < Decimal("56.0001") and SAMPLES_MAX * SAMPLE_R_MAX < Decimal(2) ** 64
    # setup: g, z check, table, trie (each node one unit), zero-fill of the trie pages, walk starts
    NODES_INT = sum((12 * Rt) ** d for d in range(1, Lt))
    NODES_LEAF = (12 * Rt) ** Lt
    INT_BYTES, LEAF_BYTES = 18 * 48 + 54 * 48, 54 * 48
    TRIE_BYTES = NODES_INT * INT_BYTES + NODES_LEAF * LEAF_BYTES
    node_items = {"copy the parent's 12 raw coordinates out of its mul.s entry into a contiguous operand (36 ldp, 36 stp)": 72,
                  "mul.s: Pi_child = Pi_parent * T[h][t] (parent's raw coordinates against the table entry)": MUL_CORE + 8,
                  "42 Karatsuba pre-sums of the child's mul.s entry, each a sum of two earlier vectors (6 ldp, 6 adds, 3 stp, address)": 42 * 20,
                  "node vector gamma: four xi products (8 modular add/sub at 45) and six pre-sums at 20": 8 * 45 + 6 * 20,
                  "store the 72 vectors (432 words, 216 stp) and address/loop control": 216 + 50}
    NODE_INSTR = sum(node_items.values())
    assert NODE_INSTR <= C
    ZERO_FILL_PER_BYTE = Fraction(2 * 131072 + 65536, 2 ** 21)
    table_units = (Decimal(Rt * (1021 + 11)) + Decimal(SAMPLE_TABLE_DRAWS * SAMPLE_MAX_INSTR) / C
                   + Decimal(2 * 12 * Rt * 800 + 12 * Rt * 54 * 200) / C)
    trie_units = Decimal(NODES_INT + NODES_LEAF) + D(ZERO_FILL_PER_BYTE * TRIE_BYTES) / C
    ONE_START = Decimal(2) ** Decimal("10.1")
    START_MAX = Decimal(1021) + Decimal(2 * SAMPLE_MAX_INSTR) / C + Decimal(43) / C
    assert START_MAX < ONE_START
    SETUP_ITEMS = (Decimal(2) ** Decimal("14.1") + 510 + table_units + trie_units + WALKS * ONE_START
                   + Decimal(32 + 2 * 1940) / C)
    SETUP_UNITS = Decimal(2) ** Decimal("49.06")
    assert SETUP_ITEMS <= SETUP_UNITS
    # walk-phase non-step work
    METER_UPDATE = 40
    restart_count_units = Decimal(RESTART_CAP) * 3 / C
    dp_units = Decimal(REPORT_CAP) * report_instr / C
    start_units = Decimal(WALKS) * METER_UPDATE / C
    restart_units = Decimal(RESTART_CAP) * METER_UPDATE / C
    restart_start_units = Decimal(RESTART_CAP) * ONE_START
    dp_restart_units = Decimal(REPORT_CAP) * (ONE_START + Decimal(METER_UPDATE) / C)
    WALK_NONSTEP = dp_units + start_units + restart_units + restart_start_units + restart_count_units + dp_restart_units
    CHUNK_K = 2 ** 63
    SWITCH_INSTR = 600
    SWITCHES = N_MAX / CHUNK_K + WALKS + RESTART_CAP + 1
    SWITCH_UNITS = SWITCHES * SWITCH_INSTR / C
    # stopping allowance: heapsort, dedup, pairs (two replays each), output, verification
    LEVEL_INSTR, OUTER_INSTR = 132, 40
    LEVELS = heapsort_level_bound(REPORT_CAP)
    assert LEVELS == 54 * 2 ** 55 + 57
    SORT_INSTR = LEVEL_INSTR * LEVELS + OUTER_INSTR * (REPORT_CAP // 2 + REPORT_CAP - 1)
    DEDUP_INSTR = 100 * REPORT_CAP
    A_MAX = 2 ** 25
    ESC_LONG_INSTR = ESC_INSTR + 64 * STEP_FLUSH_I + 65 * 64
    assert ESC_LONG_INSTR <= A_MAX
    REPLAY_ITERS = 2 * STUCK + 3
    REPLAY_ITER_INSTR = 2 * A_MAX
    REPLAY_UNITS = REPLAY_ITERS * REPLAY_ITER_INSTR / C
    PAIR_UNITS = 2 * REPLAY_UNITS + 2 * Decimal(4000) / C + 600
    OUTPUT_UNITS = 510 + 2 ** 11 + 2 ** 14   # g^x = z, [x]Q0 in G2 (proof.md section 5), verification pairing
    R_STOP_INSTR = SORT_INSTR + DEDUP_INSTR + 100 + int(((PAIR_BUDGET * PAIR_UNITS + OUTPUT_UNITS) * C)
                                                          .to_integral_value(rounding="ROUND_CEILING"))
    R_STOP = Decimal(R_STOP_INSTR) / C
    NONWALK = WALK_NONSTEP + SETUP_UNITS + R_STOP
    cd = D(c)
    # separate units (proof.md 5): useful iterations c S, the waste allowance X_ALLOW = eps c S in meter work (each
    # escape at its full cost once), the stuck allowance c E, setup, walk non-step work and the stopping allowance
    X_ALLOW = eps * cd * S
    W = cd * S + X_ALLOW + cd * E + NONWALK + SWITCH_UNITS
    assert abs(W - (cd * S * (1 + eps) + cd * E + NONWALK + SWITCH_UNITS)) <= W * Decimal(10) ** -55
    claim_thousandths = int((W.ln() / ln2 * 1000).to_integral_value(rounding="ROUND_CEILING"))
    CLAIM = Decimal(claim_thousandths) / 1000
    assert CLAIM == Decimal(CLAIM_TXT), (CLAIM, lg(W))
    W_up = W.quantize(Decimal(1).scaleb(int(W.adjusted()) - 40), rounding="ROUND_CEILING")
    fnum, fden = Fraction(W_up).numerator, Fraction(W_up).denominator
    assert fnum ** 1000 <= fden ** 1000 * 2 ** claim_thousandths, "ledger exceeds the claim"
    assert fnum ** 1000 > fden ** 1000 * 2 ** (claim_thousandths - 1), "claim is not the least thousandth"
    with localcontext() as ctx:
        ctx.prec = 90
        WCAP_INSTR = int((Decimal(2) ** CLAIM * C).to_integral_value(rounding="ROUND_FLOOR"))
    assert WCAP_INSTR ** 1000 <= 2 ** claim_thousandths * C ** 1000 < (WCAP_INSTR + 1) ** 1000
    B_WALK = WCAP_INSTR - int(SETUP_UNITS) * C - R_STOP_INSTR
    L_THRESH = B_WALK - A_MAX
    assert 0 < L_THRESH < B_WALK < 2 ** 192
    walk_need = (cd * S + X_ALLOW + cd * E + WALK_NONSTEP + SWITCH_UNITS) * C
    assert walk_need <= L_THRESH
    # meter lemma: while realized waste X <= X_ALLOW and restart plus in-progress steps <= E, the meter (pre-debit, never
    # above L_THRESH + A_MAX <= B_WALK) cannot refuse a charge before S useful iterations are paid
    assert walk_need + A_MAX <= B_WALK
    budget_units = Decimal(L_THRESH) / C
    slack_frac = (budget_units - WALK_NONSTEP - SWITCH_UNITS) / (cd * (S * (1 + eps) + E)) - 1
    EXP_REPORTS = Decimal(B_WALK) / C / cd * Decimal(2) ** (-THETA_BITS)
    assert REPORT_CAP / EXP_REPORTS > 2
    # success chain
    SAMPLER_LOSS = Decimal(2) ** Decimal("-71.99")
    assert SAMPLES_MAX * Decimal(2) ** -SAMPLE_R_MAX <= SAMPLER_LOSS
    TIE = Decimal(2) ** Decimal("-240")
    assert 70 * N_MAX / p < TIE
    # a colliding trail is size-biased; with at least U/4 distinct classes in U steps it reaches no DP within STUCK
    # steps with probability at most (a + 4) e^(-a/4) (proof.md 3.5)
    DETECT = (A_DET + 4) * (-Decimal(A_DET) / 4).exp()
    assert DETECT < Decimal("1.2e-11")
    # within-trail pairs: a trail of U <= 4 D steps, D geometric of mean 1/theta, E[U^2] <= 32/theta^2, over at most
    # MU_REPORTS + WALKS trails plus the restarted ones (<= RESTART_CAP of STUCK steps), each pair 1/N' under H1
    SELF_COL = ((MU_REPORTS + WALKS) * 32 / THETA ** 2 + RESTART_CAP * STUCK * STUCK) / (2 * Npd)
    assert SELF_COL < Decimal(2) ** -35 < Decimal("1e-6")
    S_TOL = S * (1 + eps) / (1 + T_TOL * W_BAR)
    H1PLUS = Pm - (1 - (-(S_TOL * S_TOL / (2 * Npd))).exp())
    assert 0 < H1PLUS < Decimal("0.0107")
    LOSS_ITEMS = {"H1+ tolerance: waste up to 2 x the H2 rate, collision-law drop from S to S_TOL": H1PLUS,
                  "report cap 2^55 exceeded (Markov)": REPORT_TAIL,
                  "restart cap 2^20 exceeded (Markov)": RESTART_TAIL,
                  "restart cap worst case, as a fraction of S": RESTART_FRAC,
                  "detection failure on the two colliding trails": 2 * DETECT,
                  "collision inside a single trail": Decimal("1e-6"),
                  "pair budget exceeded (Markov)": PAIR_TAIL,
                  "sampler exhaustion allowance": SAMPLER_LOSS, "key tie or zero key": TIE,
                  "beta^1 = beta^2": Decimal(2) ** -255}
    LOSS_SUM = sum(LOSS_ITEMS.values())
    SUCCESS_LB = (Pm - LOSS_SUM).quantize(Decimal("0.00001"), rounding="ROUND_FLOOR")
    assert SUCCESS_LB >= Decimal("0.39") and Pm - LOSS_SUM >= SUCCESS_LB
    # the 0.39 bound is conditional on H1+: success after every loss as a function of the tolerance t, its slope per
    # unit of t at t = T_TOL (from t = 2 to 3), and the largest t (in hundredths) at which it still reaches 0.39
    def success_at_t(t):
        st = S * (1 + eps) / (1 + t * W_BAR)
        return Pm - (Pm - (1 - (-(st * st / (2 * Npd))).exp())) - (LOSS_SUM - H1PLUS)
    T_MAX = max(k for k in range(200, 400) if success_at_t(Decimal(k) / 100) >= Decimal("0.39"))
    assert T_MAX >= 200
    # memory
    LOG_BYTES = REPORT_CAP * SLOT_BYTES
    MEM = LOG_BYTES + TRIE_BYTES + WALKS * 4096 + 2 ** 26
    assert LOG_BYTES < 2 ** 63 and NODES_LEAF * LEAF_BYTES < 2 ** 63 and MEM < 2 ** 62
    # bit layout of K (stored, below p/2 < 2^380)
    assert 64 + THETA_BITS + 112 == 262 <= 380
    # sensitivity
    def W_with(cc=cd, ee=eps, ss=S):
        return cc * ss * (1 + ee) + cc * E + NONWALK + SWITCH_UNITS
    return {
        "every_step_items": EVERY_ITEMS, "every_step_instr": EVERY, "flush_items": FLUSH_ITEMS, "flush_instr": FLUSH,
        "escape_items": ESC_ITEMS, "escape_instr": ESC_INSTR,
        "step_prices_instr": {"root": STEP_ROOT_I, "lazy": STEP_LAZY_I, "lazy_with_flush": STEP_FLUSH_I},
        "block_instr": BLOCK_I, "c_units_per_step": float(c),
        "R": Rt, "trie_depth_L": Lt, "sum_q2_conservative": float(sumq2),
        "log2_N": lg(D(N)), "log2_N_effective": lg(Npd), "P_model": float(Pm), "log2_S": lg(S),
        "eps_two_cycles_model": float(eps2), "eps_cycles_3_to_64_model": float(eps_long),
        "eps_charged_1.25x": float(eps), "escape_in_step_equivalents": float(esc_steps),
        "useless_walks_longer_than_64_rate_charged_log2": -108, "eps_doubling_checkpoint_closures_model": float(eps_brent),
        "stuck_constant_a": A_DET, "stuck_cap_steps_log2": lg(STUCK), "H1plus_tolerance_t": T_TOL,
        "log2_S_tolerance": lg(S_TOL), "report_tail_markov": float(REPORT_TAIL), "restart_tail_markov": float(RESTART_TAIL),
        "pair_budget": PAIR_BUDGET, "self_collision_bound": float(SELF_COL), "waste_allowance_units_log2": lg(X_ALLOW),
        "hard_max_walk_steps_log2": lg(N_MAX), "expected_restarts_log2": lg(EXP_RESTARTS),
        "restart_cap_log2": 20,
        "stuck_trail_probability_log2": lg(P_STUCK), "waste_per_first_visit_model": float(W_BAR),
        "S_squared_over_N_effective": float(X_COL), "success_slope_dP_dS_times_S": float(SLOPE),
        "restart_worst_case_step_fraction": float(RESTART_FRAC),
        "trie_nodes_internal_log2": math.log2(NODES_INT), "trie_leaves_log2": math.log2(NODES_LEAF),
        "trie_bytes_log2": math.log2(TRIE_BYTES), "node_items": node_items, "node_instr": NODE_INSTR,
        "trie_units_log2": lg(trie_units), "setup_items_units_log2": lg(SETUP_ITEMS),
        "preprocessing_bound_units_log2": 49.06,
        "report_log_bytes_log2": math.log2(LOG_BYTES), "peak_memory_bytes_log2": math.log2(MEM),
        "expected_reports_log2": lg(EXP_REPORTS), "report_cap_over_expected": float(REPORT_CAP / EXP_REPORTS),
        "walk_nonstep_units_log2": lg(WALK_NONSTEP), "switch_units_log2": lg(SWITCH_UNITS),
        "replay_units_per_trail_log2": lg(REPLAY_UNITS), "stopping_allowance_units_log2": lg(R_STOP),
        "nonwalk_units_log2": lg(NONWALK), "detection_window_steps_log2": lg(E),
        "log2_W": lg(W), "submitted_bound": str(CLAIM),
        "walk_budget_slack_fraction_over_need": float(slack_frac),
        "walk_budget_instr_log2": math.log2(B_WALK), "threshold_instr_log2": math.log2(L_THRESH),
        "success_loss_items": {k: float(v) for k, v in LOSS_ITEMS.items()},
        "success_loss_sum": float(LOSS_SUM), "success_after_every_loss": float(Pm - LOSS_SUM),
        "success_lower_bound_floor_5_digits": str(SUCCESS_LB),
        "success_drop_per_unit_t_from_2_to_3": float(success_at_t(Decimal(2)) - success_at_t(Decimal(3))),
        "largest_t_hundredths_keeping_0.39": T_MAX, "success_at_t_3": float(success_at_t(Decimal(3))),
        "expected_reports_at_iteration_max_log2": lg(MU_REPORTS), "expected_false_fingerprint_pairs_log2": lg(FALSE_PAIRS),
        "sensitivity_bits_eps_margin_1_instead_of_1.25": lg(W_with(ee=eps / Decimal("1.25"))) - lg(W),
        "sensitivity_bits_per_100_instr_on_every_step": lg(W_with(cc=cd + Decimal(100) / C)) - lg(W),
        "sensitivity_bits_if_lazy_kernel_doubled": lg(W_with(cc=cd + Decimal(4 * LAZY_CORE) / Decimal(6 * C))) - lg(W),
    }
