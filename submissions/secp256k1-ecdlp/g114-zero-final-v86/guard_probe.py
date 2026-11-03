import json,sys,pathlib
P=pathlib.Path(__file__).resolve().parent
if (P/'files').is_dir():P=P/'files'
sys.path.insert(0,str(P));sys.path.insert(0,str(P/'kernels'));import a64,machine
json.load(sys.stdin);prog,lab=a64.a64_load((P/'kernels'/'route.s').read_text());K=256;counts=[];cases=0
for remaining in (1,2,17,127,255,256):
 nd=machine.Node();m=nd.mem;cap=512;base=nd.R2;end=base+32*cap;physical=end+32*K;cursor=end-32*remaining;src=nd.R1;lc=nd.T2;rc=lc|machine.RC
 m[lc]=machine.PB2+32;m[rc]=cursor;m[nd.E1]=src+32*K
 for i in range(K):
  for w in range(4):m[src+32*i+8*w]=0 if w==3 else i*4+w+1
 m[physical]=0xcafebabe;port=machine.Port(m,machine.MB2,machine.PB2);out={}
 n=a64.a64_run(prog,lab,'h2s',dict(x3=1,x13=nd.E1,x14=src,x15=32*K,x12=end,sp=nd.SP,x30=0xdeadbeef),m,out_regs=out,dev=port.dev)
 assert len(port.blocks)==1 and port.blocks[0][1]==cursor
 lv=machine.LinkView(m,port.blocks[0],[]);d={};delivery=a64.a64_run(prog,lab,'dlvk',dict(x0=machine.LINK,sp=nd.SP,x30=0xdeadbeef),lv,out_regs=d)
 assert len(lv.writes)==4*K and min(lv.writes)==cursor and max(lv.writes)==cursor+32*K-8 and max(lv.writes)<physical and m[physical]==0xcafebabe
 assert all(m[cursor+32*i+8*w]==m[src+32*i+8*w] for i in range(remaining) for w in range(4))
 fport=machine.Port(m,machine.MB2,machine.PB2)
 fn=a64.a64_run(prog,lab,'fin',dict(x0=machine.MB2,x7=lc,x3=rc,x8=machine.PB2+32,x1=1,x2=nd.E2,x6=end,sp=nd.SP,x30=0xdeadbeef),m,dev=fport.dev)
 assert m[rc]==end and fport.msgs[0][2]==end and len(fport.blocks)==0
 counts.append((remaining,n,delivery,fn));cases+=1
print(json.dumps({'all_checks':True,'cases':cases,'counts':counts,'scope':'finite native-simulation full-block H2S, dlvk, fin; no complete production proof'}))
