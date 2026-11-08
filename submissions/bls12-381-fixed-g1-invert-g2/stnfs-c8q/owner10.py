""
import math
ACTIVE_OWNERS=6
C_ROW=1020.0
LA_BOUND=113.231356737
def row_tail(send=15,combine=62,redc=260,owners=ACTIVE_OWNERS):
	return(owners-1)*(send+18)+(owners-1)*(combine+2)+redc+6+10*owners+40+168
def install(model):
	assert model.C_ROW==2710.0
	assert(model.BW_M,model.BW_N)==(144,17)
	model.C_ROW=C_ROW
def check(model,result):
	owners=math.ceil(2**result['log2_N']/(1<<55))
	row=row_tail()
	return(
	owners<=ACTIVE_OWNERS
	and 2**result['log2_N']<ACTIVE_OWNERS*(1<<55)
	and row==1019
	and row<=model.C_ROW==C_ROW
	and(model.BW_M,model.BW_N)==(144,17)
	and result['log2_la']<=LA_BOUND
	)
