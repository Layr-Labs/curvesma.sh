""
import math
ACTIVE_OWNERS=7
C_ROW=1320.0
LA_BOUND=113.815127389
def row_tail(send=15,combine=62,finish=92,fold9=180,reduce6=183,owners=ACTIVE_OWNERS):
	return(owners-1)*(send+18)+(owners-1)*(combine+2)+finish+fold9+reduce6+2+10*owners+40+168
def install(model):
	assert model.C_ROW==2710.0
	assert(model.BW_M,model.BW_N)==(72,5)
	model.C_ROW=C_ROW
def check(model,result):
	owners=math.ceil(2**result['log2_N']/(1<<55))
	row=row_tail()
	return(
	owners<=ACTIVE_OWNERS
	and 2**result['log2_N']<7*(1<<55)
	and row==1317
	and row<=model.C_ROW==C_ROW
	and(model.BW_M,model.BW_N)==(72,5)
	and result['log2_la']<=LA_BOUND
	)
