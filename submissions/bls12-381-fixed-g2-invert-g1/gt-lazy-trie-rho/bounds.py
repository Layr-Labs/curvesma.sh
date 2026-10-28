"""Success bounds of proof.md sections 3.4 and 4, computed without independence across steps, trails or reports:
expectations by deferred decisions (a fresh class's labels are uniform given everything revealed, H2), H1 for merges, and linearity;
cap overflows by Markov; the waste of fruitless closed walks by Freedman's martingale inequality. Loaded by
experiment.py part D; every figure the text quotes is asserted here."""

from decimal import Decimal as Dec
from fractions import Fraction


def freedman_log2(lam, v, b, split):
    """log2 of split * exp(-lam^2 / (2 (v + b lam / 3))): Freedman's bound on each of `split` interleaved
    martingales (increments <= b, predictable variance <= v, deviation lam), union over the split."""
    nats = lam * lam / (2 * (v + b * lam / 3))
    return (Dec(split).ln() - nats) / Dec(2).ln()


def bounds(D, lg, theta_bits, n_max, S, Npd, P2, P_l, cmax_c, esc, r_model, eps_charged, rep_cap, rst_cap, pair_cap,
           long_rate, walks):
    theta = Dec(2) ** -theta_bits
    Pl = [Fraction(p).limit_denominator(10 ** 30) for p in P_l]   # as experiment.py charges it
    # the weights of section 3.3, in step equivalents: 2 + esc for a 2-cycle, 64 + l + l cmax/c + esc for length l;
    # their sum weighted by the model rates is r = 0.026346, and the charged eps is 1.25 r
    w2 = 2 + esc
    wl = [64 + l + l * cmax_c + esc for l in range(3, 65)]
    assert P2 * w2 + sum(p * w for p, w in zip(Pl, wl)) == r_model, "section 3.3 weights"
    dec = lambda x: x if isinstance(x, Dec) else D(x)   # experiment.py passes eps as the Decimal it charges
    assert dec(eps_charged) == D(Fraction(5, 4) * r_model), "eps is 1.25 x the model rate"
    b = D(wl[-1])          # the largest charge one iteration can trigger, in step equivalents
    lam = (dec(eps_charged) - D(r_model)) * S
    waste_log2 = freedman_log2(lam / 128, b * D(r_model) * S, b, 128)
    assert waste_log2 < -1000
    r_it = 2 * D(P2) + sum(D(p) * (64 + l) for l, p in zip(range(3, 65), Pl))
    m = 19 / theta
    trail_log2 = freedman_log2((1 / theta - r_it * m) / 128, 128 * r_it * m, b, 128)
    assert r_it < Dec("0.01112") and r_it * 19 < Dec("0.2113") and trail_log2 < -1000
    trail_tail = Dec(2) ** -1000
    # apart from a merge (merges below) and a useless closed walk over 64 (exp_long), a trail is restarted only if
    # its first 19/theta fresh classes hold no distinguished class ((1 - theta)^m <= e^-19, one fresh label at a
    # time) or it wastes more than 1/theta of its 20/theta iterations (trail_tail)
    stuck_each = (-Dec(19)).exp() + trail_tail
    detect_each = 20 * (-Dec(19)).exp() * (1 + theta) + trail_tail
    assert stuck_each < Dec("5.61e-9") and detect_each < Dec("1.13e-7")
    merges = n_max * n_max / (2 * Npd)
    assert merges < Dec("0.546")
    exp_rep = n_max * theta + 1
    exp_stuck = (exp_rep + walks + rst_cap) * stuck_each + merges
    exp_long = n_max * long_rate
    exp_rst = exp_stuck + exp_long
    assert lg(exp_stuck) < 18.32
    false_pairs = Dec(rep_cap) * rep_cap / 2 / Dec(2) ** 112
    caps = {"reports": (rep_cap, exp_rep), "restarts": (rst_cap, exp_rst), "pairs": (pair_cap, false_pairs)}
    ratio = {k: Dec(c) / e for k, (c, e) in caps.items()}
    share = {k: e / Dec(c) for k, (c, e) in caps.items()}
    assert min(ratio.values()) >= 500
    single = 2 * 20 / theta / S
    assert single < Dec("1e-6")

    def window(spare_instr, c, C, esc_max):   # proof.md 3.4: window iterations at c plus one escape each
        iters = walks * (20 / theta + Dec(2) ** 63)
        units, spare = iters * (c + Dec(esc_max) / C), Dec(spare_instr) / C
        assert lg(iters) < 104.33 and lg(units) < 110.09 and lg(spare) > 111.05
        return {"iterations_log2": lg(iters), "units_log2": lg(units), "spare_units_log2": lg(spare)}

    return {"window": window, "b": b, "lam": lam, "waste_log2": waste_log2, "r_it": r_it, "trail_log2": trail_log2,
            "stuck_each": stuck_each, "detect_each": detect_each, "merges": merges, "exp_rep": exp_rep,
            "exp_stuck": exp_stuck, "exp_long": exp_long, "exp_rst": exp_rst, "false_pairs": false_pairs,
            "ratio": ratio, "share": share, "single": single}
