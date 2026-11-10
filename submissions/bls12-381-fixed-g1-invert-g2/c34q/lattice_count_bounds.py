"";from fractions import Fraction as F;from math import isqrt;from ideal_canonical import determinant;from common12 import integer_metric;from itertools import product
def atan_interval(q,terms=100):value=sum((F((-1)**k,(2*k+1)*q**(2*k+1))for k in range(terms)),F(0));remainder=F((-1)**terms,(2*terms+1)*q**(2*terms+1));return min(value,value+remainder),max(value,value+remainder)
def pi_interval():a,b=atan_interval(5),atan_interval(239);return 16*a[0]-4*b[1],16*a[1]-4*b[0]
def sqrt_interval(value,bits=128):
	value=F(value)
	if value<0:raise ValueError('negative square root')
	scale=1<<bits;k=isqrt((value.numerator*scale*scale)//value.denominator);lower=F(k,scale);upper=lower if lower*lower==value else F(k+1,scale);assert lower*lower<=value<=upper*upper;return lower,upper
def determinant_rational(matrix):q,denominator,_=integer_metric(matrix);return F(determinant(q),denominator**len(q))
def cell_radius(matrix):
	m=[[F(v)for v in row]for row in matrix];vmax=max(sum(v[i]*m[i][j]*v[j]for i in range(6)for j in range(6))for v in product((-1,1),repeat=6))/4;k=isqrt(int(vmax*4096))
	while F(k,64)**2<vmax:k+=1
	return F(k,64),vmax
def half_count_ratios(matrix,minimum_radius=1024):
	matrix=[[F(v)for v in row]for row in matrix]
	if len(matrix)!=6 or any(len(row)!=6 for row in matrix):raise ValueError('six-dimensional metric')
	if any(matrix[i][j]!=matrix[j][i]for i in range(6)for j in range(6)):raise ValueError('symmetric metric required')
	if not all(determinant_rational([row[:k]for row in matrix[:k]])>0 for k in range(1,7)):raise ValueError('positive-definite metric required (Sylvester)')
	t,_=cell_radius(matrix)
	if minimum_radius<=t:raise ValueError('radius too small')
	p=pi_interval();v6=p[0]**3/6,p[1]**3/6;v5=F(8,15)*p[0]**2,F(8,15)*p[1]**2;d6=determinant_rational(matrix);d5=determinant_rational([row[:5]for row in matrix[:5]]);sd6,sd5=sqrt_interval(d6),sqrt_interval(d5);section_coefficient_upper=v5[1]/v6[0]*sd6[1]/sd5[0];r=F(minimum_radius);lower=(1-t/r)**6-section_coefficient_upper*(1+t/r)**5/r;upper=(1+t/r)**6;return lower,upper,t
def correction(cfg,minimum_radius=1024):
	a=half_count_ratios(cfg['matrix_a'],minimum_radius);b=half_count_ratios(cfg['matrix_b'],minimum_radius);yield_factor=a[0]*b[0];work_factor=a[1]*b[1]
	return dict(a_count_ratio=[float(x)for x in a[:2]],
	b_count_ratio_uniform=[float(x)for x in b[:2]],
	safe_yield_factor=yield_factor,safe_work_population_factor=work_factor)
