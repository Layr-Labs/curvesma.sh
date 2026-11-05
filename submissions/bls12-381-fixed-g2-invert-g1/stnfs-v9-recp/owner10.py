"""Ten active vector owners (N < 10 x 2^56 = 2^59.32 at merge weight 300, excess 1.05); blocks (28, 5) (v7)."""
import math
ACTIVE_OWNERS = 10
C_ROW = 1650.0
LA_BOUND = 117.925074783
def row_tail(send=15, combine=62, finish=92, fold9=180, reduce6=183, owners=ACTIVE_OWNERS):
    # (owners - 1) sends and combinations up the fixed tree, finish at the row owner, 10 header
    # instructions per owner, routing 40 and group overhead 168 (STORAGE.md)
    return (owners - 1) * (send + 18) + (owners - 1) * (combine + 2) + finish + fold9 + reduce6 + 2 + 10 * owners + 40 + 168
def install(model):
    assert model.C_ROW == 2710.0
    assert (model.BW_M, model.BW_N) == (28, 5)
    model.C_ROW = C_ROW
def check(model, result):
    # The production bound log2(N) < log2(10) + 56 implies N < 10*2^56; the allocation test is applied to the model N.
    owners = math.ceil(2 ** result['log2_N'] / (1 << 56))
    row = row_tail()
    return (
        owners <= ACTIVE_OWNERS
        and result['log2_N'] < math.log2(ACTIVE_OWNERS) + 56
        and row == 1638
        and row <= model.C_ROW == C_ROW
        and (model.BW_M, model.BW_N) == (28, 5)
        and result['log2_la'] <= LA_BOUND
    )
