""
import math
ACTIVE_OWNERS=3
C_ROW=700.0
LA_BOUND=112.821878470
def row_tail(send=15,combine=62,redc=260,part=267,owners=ACTIVE_OWNERS):
	return(owners-1)*(send+18)+(owners-1)*(combine+2)+redc+6+10*owners+40+owners*math.ceil((part+6)/5)
def install(model):
	assert model.C_ROW==2710.0
	assert(model.BW_M,model.BW_N)==(144,17)
	model.C_ROW=C_ROW
def check(model,result):
	per=2**64//320
	owners=math.ceil(2**result['log2_N']/per)
	row=row_tail()
	return(
	owners<=ACTIVE_OWNERS
	and 2**result['log2_N']<ACTIVE_OWNERS*per
	and row==695
	and row<=model.C_ROW==C_ROW
	and(model.BW_M,model.BW_N)==(144,17)
	and result['log2_la']<=LA_BOUND
	)
