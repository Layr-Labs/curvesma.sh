"""Finite cyclic-group r-adding rho checks; not production security evidence."""
import hashlib,json,math,random,sys
params=json.load(sys.stdin)
trials=params.get('trials',64)
if not isinstance(trials,int) or not 1<=trials<=128: raise ValueError('trials')
orders=[1019,10007,65537,1000003]
results=[]
for order in orders:
 successes=fruitless=timeouts=0; lengths=[]
 for trial in range(trials):
  rng=random.Random(10000000*order+trial)
  secret=rng.randrange(order)
  recovered=None; total=0
  for attempt in range(8):
   seed=rng.getrandbits(256).to_bytes(32,'big')
   pairs=[(rng.randrange(order),rng.randrange(order)) for _ in range(32)]
   table=[(a+b*secret)%order for a,b in pairs]
   a,b=rng.randrange(order),rng.randrange(order)
   initial=((a+b*secret)%order,a,b)
   def step(state):
    x,a,b=state
    index=hashlib.sha256(seed+x.to_bytes(32,'big')).digest()[0]&31
    da,db=pairs[index]
    return ((x+table[index])%order,(a+da)%order,(b+db)%order)
   slow=fast=initial
   for iteration in range(16*math.isqrt(order)+16):
    slow=step(slow);fast=step(step(fast));total+=3
    if slow[0]==fast[0]:
     denominator=(slow[2]-fast[2])%order
     if denominator:
      x=(fast[1]-slow[1])*pow(denominator,-1,order)%order
      if x==secret: recovered=x
      else: raise AssertionError('incorrect collision recovery')
     else: fruitless+=1
     break
   else: timeouts+=1
   if recovered is not None:break
  if recovered==secret:successes+=1
  lengths.append(total)
 results.append({'order':order,'trials':trials,'successes':successes,'fruitless_restarts':fruitless,'capped_attempts':timeouts,'max_total_steps':max(lengths),'mean_total_steps':sum(lengths)/trials})
print(json.dumps({'scope':'toy additive prime-order groups; deterministic test seeds; not production validity or scalar-instruction measurements','partitions':32,'attempts':8,'results':results}))
