import json,sys,math,bisect
q=json.load(sys.stdin)['order']
method='bsgs'
m=math.isqrt(q-1)+1
checked=0
for base in range(1,q):
 if math.gcd(base,q)!=1: continue
 for secret in range(q):
  target=base*secret%q
  if method=='enumeration':
   value=0
   for answer in range(q):
    if value==target: break
    value=(value+base)%q
   else: raise AssertionError('missing')
  else:
   babies=sorted((j*base%q,j) for j in range(m))
   keys=[v for v,j in babies]
   value=target
   for i in range((q+m-1)//m):
    k=bisect.bisect_left(keys,value)
    if k<len(keys) and keys[k]==value:
     answer=i*m+babies[k][1]
     if answer<q: break
    value=(value-m*base)%q
   else: raise AssertionError('missing')
  assert answer==secret
  assert base*answer%q==target
  checked+=1
print(json.dumps({'method':method,'order':q,'checked':checked,'all_recovered':True,'scope':'abstract cyclic groups only; no production performance claim'}))
