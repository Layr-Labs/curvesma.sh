"""Finite and numerical checks only; never executes the production NFS attack."""
import json, sys
from decimal import Decimal, localcontext

def main():
    # Protocol accepts bounded input but no input controls filesystem or processes.
    text=sys.stdin.read()
    if text.strip(): json.loads(text)
    counts={}
    for r in (3,5,7,11,17,31):
        image={(b,(b*a)%r) for b in range(1,r) for a in range(r)}
        assert len(image)==r*(r-1)
        assert all((b,x) in image for b in range(1,r) for x in range(r))
        counts[str(r)]=len(image)
    with localcontext() as ctx:
        ctx.prec=90
        D=Decimal
        pow2=lambda x:(D(x)*D(2).ln()).exp()
        core=D('1.01')*pow2('122.0169')
        wrapper=pow2('80')
        p=core+wrapper/D(20000)
        g=D(10)*core+wrapper/D(2000)
        success=(D(1)-D(2)**-128)*D('0.391663')-D(2)**-48
        checks={'toy_distribution_bijections':True,
          'recover_p_budget':p<pow2('122.04'),
          'g1_plain_budget':g<pow2('125.3531834'),
          'sampling_success':success>D('.39')}
        assert all(checks.values())
        print(json.dumps({'checks':checks,'toy_pair_counts':counts,
          'g1_plain_log2_upper_model':str(g.ln()/D(2).ln()),
          'recover_p_log2_upper_model':str(p.ln()/D(2).ln()),
          'success_lower_model':str(success),
          'scope':'Conditional conversion arithmetic and finite toy bijections only. Does not validate source H1-H5, wrapper implementation, production success or kernels.'}))
if __name__=='__main__': main()
