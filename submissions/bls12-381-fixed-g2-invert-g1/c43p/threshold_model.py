"";from math import exp,log,log2,log1p,ceil,floor,sqrt;from collections import Counter;from array import array;LN2=log(2)
def primes_upto(n):
	sieve=bytearray(b'\1')*(n+1);sieve[:2]=b'\0\0'
	for p in range(2,int(sqrt(n))+1):
		if sieve[p]:sieve[p*p::p]=b'\0'*((n-p*p)//p+1)
	return[p for p in range(2,n+1)if sieve[p]]
def coefficients(weights,XN,length):
	terms=[(j,j*v)for j,v in enumerate(weights)if j and abs(v)>=1e-24];out=[exp(-XN)]+[.0]*length
	for n in range(1,length+1):
		value=.0
		for j,w in terms:
			if j>n:break
			value+=w*out[n-j]
		out[n]=value/n
		if value<-1e-13:raise ArithmeticError('negative coefficient beyond roundoff tolerance')
	return out
class Retention:
	def __init__(self,XI,XK,*,prime_min=256,step=.25,max_bits=1024,tilt=.05,explicit_bits=12,power_terms=8,odd=0):
		if not 1<prime_min<=2**XK<=2**XI:raise ValueError('invalid sieve/factor-base range')
		if not 0<step<=1 or max_bits<=XI or tilt>=LN2:raise ValueError('invalid grid/tilt')
		self.B,self.X,self.minimum,self.step=XI,XK,prime_min,step;self.tilt,self.length=tilt,ceil(max_bits/step);self.names='small','sieve_geo','large_geo','sieve_sf','large_sf';weights={k:[.0]*(self.length+1)for k in self.names};XG=Counter();XL=min(float(explicit_bits),XI)
		def XQ(p):
			b=log2(p);XF=round(b/step);q=exp(tilt*XF*step-LN2*b)
			if not 0<q<1:raise ValueError('tilt makes prime series divergent')
			if p<prime_min:XS=[('small',False)]
			else:XH='sieve'if b<XK else 'large';XS=[(XH+'_geo',False),(XH+'_sf',True)]
			for name,squarefree in XS:
				XG[name]+=log1p(q)if squarefree else-log1p(-q);qk=q
				for k in range(1,self.length+1):
					index=k*XF
					if index>self.length:break
					sign=-1 if squarefree and k%2==0 else 1;weights[name][index]+=sign*qk/k;qk*=q
					if qk<1e-25:break
		for p in primes_upto(int(2**XL)-1)[odd:]:XQ(p)
		def XJ(XT,XR,XH):
			if XT>=XR:return
			width=step/4;cells=ceil((XR-XT)/width);width=(XR-XT)/cells
			for i in range(cells):
				b=XT+(i+.5)*width;XF=round(b/step)
				for k in range(1,power_terms+1):
					XO=exp(k*tilt*XF*step-(k-1)*LN2*b)*width/(k*b);index=k*XF;sign=-1 if k%2==0 else 1
					for g,x in((XH,XO),)if XH=='small'else((XH+'_geo',XO),(XH+'_sf',sign*XO)):
						XG[g]+=x
						if index<=self.length:weights[g][index]+=x
		lm=log2(prime_min);XJ(XL,min(lm,XI),'small');XJ(max(XL,lm),min(XK,XI),'sieve');XJ(max(XL,XK),XI,'large');self.weights,self.normalizers=weights,dict(XG)
		XM=[sum(weights[k][j]for k in('small','sieve_geo','large_geo'))for j in range(self.length+1)]
		self.denominator=coefficients(XM,sum(XG[k]for k in('small','sieve_geo','large_geo')),self.length)
		self.selection_normalizer=exp(XG['sieve_sf']+XG['large_sf']-XG['sieve_geo']-XG['large_geo'])
		self._c={}
	def co(self,k):
		if k not in self._c:self._c[k]=coefficients(self.weights[k],self.normalizers[k],self.length)
		return self._c[k]
	def two(self,n,a,b=9**9,win=None,sf=False):
		"";a=min(a,n);key=n,b,win,sf
		if a<0 or b<0:return.0
		if key not in self._c:
			if key[1:]not in self._c:
				st=self.step;sm=self.co('small');D=min(self.length,2048);m=[.0]*(D+1)
				for l,v in enumerate(self.co('large_sf'if sf else'large_geo')[:min(b,D)+1]):
					if v and(win is None or l*st<=win[2]or any(self.X*k-win[2]<=l*st<=win[1]*k+win[2]for k in range(1,win[0]+1))):
						for t in range(D+1-l):m[l+t]+=v*sm[t]
				self._c[key[1:]]=m
			m=self._c[key[1:]];sv=self.co('sieve_sf'if sf else'sieve_geo');c=array('d');t=.0
			for x in range(min(n,len(m)-1)+1):t+=m[x]*sv[n-x];c.append(t)
			self._c[key]=c
		c=self._c[key];return c[min(a,len(c)-1)]
	def cap(self,nb):return floor(ceil(nb)/log2(self.minimum))
	def at(self,nb,sl,sb=None):
		if nb<=0 or not 0<=sl<500 or nb/self.step>self.length:raise ValueError('norm or slack outside model grid')
		n=round(nb/self.step);rb=(self.cap(nb)+1)/16;r=self.selection_normalizer*self.two(n,floor((sl-rb)/self.step),9**9 if sb is None else floor((sb-rb-.125)/self.step),None,True)/self.denominator[n]
		if not-1e-10<=r<=1+1e-9:raise ArithmeticError('retention outside probability range')
		return dict(retention=max(.0,min(1.,r)))
	def survivor_upper(self,nb,sl,ex=0,sb=None,lb=None,kmax=4):
		n=round(nb/self.step)
		if not 0<n<=self.length or not 0<=sl+ex<500 or ex<0:raise ValueError('invalid survivor range')
		a=floor((sl+ex)/self.step);return(self.two(n,a)if sb is None else self.two(n,a,floor((sb+.25)/self.step),(kmax,lb,(self.cap(nb)+2)/16+.25)))/self.denominator[n]
