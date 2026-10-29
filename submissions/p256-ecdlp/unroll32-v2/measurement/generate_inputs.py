import random, sys
P=2**256-2**224+2**192+2**96-1
A=P-3
B=int("5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",16)
n,seed=map(int,sys.argv[1:]);rng=random.Random(seed);count=0
while count<n:
 x=rng.randrange(P);z=(x*x*x+A*x+B)%P;y=pow(z,(P+1)//4,P)
 if y*y%P!=z:continue
 print(f"{x:064x}\n{y:064x}");count+=1
