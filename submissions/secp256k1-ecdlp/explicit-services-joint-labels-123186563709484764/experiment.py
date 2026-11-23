"""Lossless stdlib runtime unpacking; all execution is inside the verifier.

runtime.json binds the readable files and the exact restored test/reference
sources. runtime-*.txt is ordinary base85(XZ(UTF-8 JSON)), without a delta codec.
The optional inspect mode returns one complete restored file and its digest.
"""
import base64
import contextlib
import hashlib
import io
import json
import lzma
import os
import re
from pathlib import Path, PurePosixPath
import runpy
import shutil
import sys
import tempfile

LIMIT=2*1024*1024
MANIFEST_SHA256="307a7ad4d3d3b9fbe79e213381340af3a744b28b389051fd7370fc600de4a646"

def unique(pairs):
    result={}
    for k,v in pairs:
        if k in result:raise ValueError('duplicate JSON member')
        result[k]=v
    return result

def valid_name(name):
    if not isinstance(name,str) or '\x00' in name or '\\' in name:
        raise ValueError('unsafe member name')
    p=PurePosixPath(name)
    if p.is_absolute() or any(x in ('','..','.') for x in name.split('/')):
        raise ValueError('unsafe member path')
    return name

def sha(data):return hashlib.sha256(data).hexdigest()

def decode(packed,manifest):
    if len(packed)>LIMIT:raise ValueError('packed limit')
    compressed=base64.b85decode(packed)
    dec=lzma.LZMADecompressor(format=lzma.FORMAT_XZ,memlimit=32*1024*1024)
    raw=dec.decompress(compressed,max_length=LIMIT+1)
    if len(raw)>LIMIT:raise ValueError('decompressed limit')
    if not dec.eof:raise ValueError('truncated or excess stream')
    if dec.unused_data:raise ValueError('trailing stream data')
    members=json.loads(raw.decode('utf-8'),object_pairs_hook=unique)
    if not isinstance(members,dict):raise ValueError('member object required')
    for name in members:valid_name(name)
    if set(members)!=set(manifest['members']):raise ValueError('member whitelist')
    if len(members)>32:raise ValueError('member count limit')
    result={}
    for name,text in members.items():
        if not isinstance(text,str):raise ValueError('member text required')
        data=text.encode('utf-8');meta=manifest['members'][name]
        if len(data)!=meta['bytes'] or len(data)>65536:raise ValueError('member length: '+name)
        if sha(data)!=meta['sha256']:raise ValueError('member hash: '+name)
        result[name]=data
    return result

def write_new(root,name,data):
    valid_name(name)
    target=root/name
    target.parent.mkdir(parents=True,exist_ok=True)
    # root is a private fresh directory. A repeated path fails closed rather
    # than following or replacing any pre-existing file or symlink.
    with target.open('xb') as f:f.write(data)

def controls(packed,manifest):
    got={}
    def reject(label,payload,meta=manifest):
        try:decode(payload,meta)
        except (ValueError,lzma.LZMAError,UnicodeError) as exc:got[label]=str(exc)
        else:raise AssertionError('damaged archive accepted: '+label)
    raw=base64.b85decode(packed)
    reject('truncated',base64.b85encode(raw[:-1]))
    reject('trailing',base64.b85encode(raw+b'x'))
    reject('concatenated',base64.b85encode(raw+lzma.compress(b'{}',preset=0)))
    damaged=bytearray(raw);damaged[len(damaged)//2]^=1
    reject('damaged',base64.b85encode(damaged))
    reject('expansion_limit',base64.b85encode(lzma.compress(b' '* (LIMIT+1),preset=0)))
    member=next(iter(manifest['members']))
    for field,value in (('bytes',manifest['members'][member]['bytes']+1),('sha256','0'*64)):
        meta=json.loads(json.dumps(manifest));meta['members'][member][field]=value
        reject('wrong_'+field,packed,meta)
    for name in ('../escape','/absolute','nul\x00name','unexpected.py'):
        payload=base64.b85encode(lzma.compress(json.dumps({name:''}).encode(),preset=0))
        reject('path_'+repr(name),payload)
    duplicate=b'{"x":"","x":""}'
    reject('duplicate',base64.b85encode(lzma.compress(duplicate,preset=0)))
    # UTF-8/newlines are not canonicalized by the codec.
    sample='\ufeffa\r\nb\n\u00e9\u03bb'
    data=sample.encode('utf-8')
    meta={'members':{'sample.txt':{'bytes':len(data),'sha256':sha(data)}}}
    samplepack=base64.b85encode(lzma.compress(json.dumps({'sample.txt':sample}).encode(),preset=0))
    assert decode(samplepack,meta)['sample.txt']==data
    return got

def main():
    source=Path('/snapshot/files') if Path('/snapshot/files/runtime.json').is_file() else Path(__file__).parent
    manifest_bytes=(source/'runtime.json').read_bytes()
    if sha(manifest_bytes)!=MANIFEST_SHA256:raise ValueError('manifest hash')
    manifest=json.loads(manifest_bytes,object_pairs_hook=unique)
    for name in manifest['chunks']:
        valid_name(name)
        if not re.fullmatch(r'runtime-[0-9]+[.]txt',name):raise ValueError('chunk name')
    if len(set(manifest['chunks']))!=len(manifest['chunks']):raise ValueError('duplicate chunks')
    packed=b''.join((source/name).read_bytes() for name in manifest['chunks'])
    if sha(packed)!=manifest['packed_sha256']:raise ValueError('packed hash')
    restored=decode(packed,manifest)
    cfg=json.load(sys.stdin)
    mode=cfg.pop('mode','run')
    if mode=='inspect':
        name=cfg['member']
        if name not in restored:raise ValueError('inspection whitelist')
        data=restored[name]
        text=data.decode('utf-8');offset=cfg.get('offset',0);limit=cfg.get('limit',60000)
        assert 0<=offset<=200000 and 1<=limit<=60000
        print(json.dumps({'member':name,'bytes':len(data),'sha256':sha(data),'manifest_sha256':MANIFEST_SHA256,'packed_sha256':sha(packed),'offset_characters':offset,'total_characters':len(text),'source':text[offset:offset+limit]}))
        return
    if mode not in ('run','roundtrip','costs'):raise ValueError('mode')
    negative=controls(packed,manifest)
    work=Path(tempfile.mkdtemp(prefix='curvesmash-runtime-'))
    try:
        for name,meta in manifest['plain'].items():
            valid_name(name);data=(source/name).read_bytes()
            if len(data)!=meta['bytes'] or sha(data)!=meta['sha256']:raise ValueError('plain hash: '+name)
            write_new(work,name,data)
        for name,data in restored.items():write_new(work,name,data)
        # Verify the restored files by read-back hashing.
        hashes={name:sha((work/name).read_bytes()) for name in restored}
        assert hashes=={name:meta['sha256'] for name,meta in manifest['members'].items()}
        try:write_new(work,next(iter(restored)),b'overwrite')
        except FileExistsError:negative['preexisting']='O_EXCL refused existing file'
        else:raise AssertionError('overwrote restored source')
        if mode=='roundtrip':result={'ok':True,'scope':'Exact source roundtrip and codec controls only; no attack suite executed.'}
        elif mode=='costs':
            sys.dont_write_bytecode=True
            sys.path[:]=[str(work),str(work/'kernels')]+[p for p in sys.path if p and '/snapshot' not in p and p!=str(source)]
            import costs
            result={'ok':True,'scope':'Exact costs only; current kernels/protocol are not checked by this mode.','costs':costs.summary()}
        else:
            # No package/cwd import fallback: the restored tree supplies every
            # nonstdlib module used by this suite. No source is rewritten here.
            sys.dont_write_bytecode=True
            sys.path[:]=[str(work),str(work/'kernels')]+[p for p in sys.path if p and '/snapshot' not in p and p!=str(source)]
            original_stdin=sys.stdin
            sys.stdin=io.StringIO(json.dumps(cfg));capture=io.StringIO()
            try:
                with contextlib.redirect_stdout(capture):runpy.run_path(str(work/'suite.py'),run_name='__main__')
            finally:sys.stdin=original_stdin
            result=json.loads(capture.getvalue())
            assert result['all_checks_pass']
        result['runtime_archive']={'packed_sha256':sha(packed),'members':hashes,'negative_controls':negative,
            'uncompressed_bytes':sum(len(b) for b in restored.values()),'codec':'stdlib XZ JSON UTF8, base85; no deltas'}
        print(json.dumps(result,separators=(',',':')))
    finally:shutil.rmtree(work)

if __name__=='__main__':main()
