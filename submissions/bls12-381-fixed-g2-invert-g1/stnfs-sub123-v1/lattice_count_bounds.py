"""Exact rational corrections from real ellipsoid volumes to lattice counts."""
from fractions import Fraction as F
from math import isqrt
from ideal_canonical import determinant
from common12 import integer_metric
def atan_interval(q, terms=100):
    value = sum((F((-1)**k, (2*k+1)*q**(2*k+1)) for k in range(terms)), F(0))
    remainder = F((-1)**terms, (2*terms+1)*q**(2*terms+1))
    return min(value, value+remainder), max(value, value+remainder)
def pi_interval():
    a, b = atan_interval(5), atan_interval(239)
    return 16*a[0]-4*b[1], 16*a[1]-4*b[0]
def sqrt_interval(value, bits=128):
    value = F(value)
    if value < 0:
        raise ValueError('negative square root')
    scale = 1 << bits
    k = isqrt((value.numerator*scale*scale)//value.denominator)
    lower = F(k, scale)
    upper = lower if lower*lower == value else F(k+1, scale)
    assert lower*lower <= value <= upper*upper
    return lower, upper
def determinant_rational(matrix):
    q, denominator, _ = integer_metric(matrix)
    return F(determinant(q), denominator**len(q))
def half_count_ratios(matrix, minimum_radius=1024, cell_radius=2):
    """Bounds valid at EVERY R>=minimum_radius, with positive last coordinate.
    Unit cubes have Q-norm displacement <=1/2 sqrt(sum|Qij|)<=cell_radius.
    Apply covering bounds in 6D and to the x5=0 section in 5D.
    """
    matrix = [[F(v) for v in row] for row in matrix]
    if len(matrix) != 6 or any(len(row) != 6 for row in matrix):
        raise ValueError('six-dimensional metric')
    if any(matrix[i][j] != matrix[j][i] for i in range(6) for j in range(6)):
        raise ValueError('symmetric metric required')
    if min(matrix[i][i]-sum(abs(matrix[i][j]) for j in range(6) if i != j)
           for i in range(6)) <= 0:
        raise ValueError('positive diagonal-dominance certificate required')
    if sum(abs(v) for row in matrix for v in row) > 4*cell_radius**2:
        raise ValueError('uncertified cell radius')
    if minimum_radius <= cell_radius:
        raise ValueError('radius too small')
    p = pi_interval()
    v6 = p[0]**3/6, p[1]**3/6
    v5 = F(8,15)*p[0]**2, F(8,15)*p[1]**2
    d6 = determinant_rational(matrix)
    d5 = determinant_rational([row[:5] for row in matrix[:5]])
    if d6 <= 0 or d5 <= 0:
        raise ValueError('positive-definite metrics required')
    sd6, sd5 = sqrt_interval(d6), sqrt_interval(d5)
    section_coefficient_upper = v5[1]/v6[0]*sd6[1]/sd5[0]
    r, t = F(minimum_radius), F(cell_radius)
    lower = (1-t/r)**6-section_coefficient_upper*(1+t/r)**5/r
    upper = (1+t/r)**6
    return lower, upper
def correction(cfg):
    a = half_count_ratios(cfg['matrix_a'])
    b = half_count_ratios(cfg['matrix_b'])
    assert a[0] >= F(98,100) and b[0] >= F(98,100)
    assert a[1] <= F(102,100) and b[1] <= F(102,100)
    return dict(a_count_ratio=[float(x) for x in a],
                b_count_ratio_uniform=[float(x) for x in b],
                safe_yield_factor='2401/2500', safe_work_population_factor='2601/2500',
                minimum_radius=1024, cube_radius=2,
                certification='exact rational Machin pi and determinant/sqrt intervals')
