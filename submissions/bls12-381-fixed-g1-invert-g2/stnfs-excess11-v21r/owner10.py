"""Ten active vector owners: allocation and row-tail formula checks."""
import math
ACTIVE_OWNERS = 10
C_ROW = 1640.0
def install(model):
    assert model.C_ROW == C_ROW
    assert (model.BW_M, model.BW_N) == (19, 4)
    model.C_ROW = C_ROW
def check(model, result):
    owners = math.ceil(2 ** result['log2_N'] / (1 << 56))
    row = (
        9 * (15 + 18) + 9 * (62 + 2)
        + 92 + 180 + 183 + 2
        + 10 * ACTIVE_OWNERS + 40 + 168
    )
    return (
        owners <= ACTIVE_OWNERS
        and result['log2_N'] < 59.12
        and row == 1638
        and row <= model.C_ROW == C_ROW
        and (model.BW_M, model.BW_N) == (19, 4)
        and result['log2_la'] <= 117.7461
    )
