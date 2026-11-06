"""Eight active vector owners (N < 8 x 2^56 = 2^59 at merge weight 300); blocks (28, 5)."""
import math
ACTIVE_OWNERS=8
C_ROW=1430.0
LA_BOUND=116.978857066
def row_tail(send=15,combine=62,finish=92,fold9=180,reduce6=183,owners=ACTIVE_OWNERS):
		return(owners-1) * (send+18) + (owners-1) * (combine+2) +finish+fold9+reduce6+2+10*owners+40+168
def install(model):
		assert model.C_ROW==2710.0
		assert(model.BW_M,model.BW_N) == (30,6)
		model.C_ROW=C_ROW
def check(model,result):
		owners=math.ceil(2**result['log2_N'] / (1<<56))
		row=row_tail()
		return(
		  owners<=ACTIVE_OWNERS
		  and result['log2_N'] <59.0
		  and row==1424
		  and row<=model.C_ROW==C_ROW
		  and(model.BW_M,model.BW_N) == (30,6)
		  and result['log2_la'] <=LA_BOUND
		  )
