"""Bounded lossless source restoration inside the verifier's isolated runner."""
import base64, contextlib, hashlib, io, json, pathlib, runpy, shutil, sys, tempfile, zlib
NAMES = ("kernels/a64.py", "kernels/generic.s", "kernels/leaf.s", "kernels/pointer_init.s", "kernels/route.s", "kernels/secp_fe.s", "kernels/ystage.s", "kernels/ytransport.s", "experiment.py")
# Two archives keep each file within 65,536 bytes: the eight kernel sources, then experiment.py.
ARCHIVES = (("runtime.json", NAMES[:8]), ("runtime2.json", NAMES[8:]))
HERE = next(p for p in (pathlib.Path(__file__).resolve().parent, pathlib.Path("/snapshot/files"), pathlib.Path.cwd()) if (p/"runtime.json").is_file())
def restore(archive, names):
    assert set(archive)=={"format","members"} and archive["format"]=="zlib-b85-kernels-v1"
    assert len(archive["members"])==len(names)
    result={}
    for m in archive["members"]:
        assert set(m)=={"name","size","sha256","data"}
        name=m["name"]
        assert name in names and name not in result
        assert type(m["size"]) is int and 0<m["size"]<=65536
        assert len(m["data"])<=65536
        d=zlib.decompressobj(); raw=d.decompress(base64.b85decode(m["data"]),65537)
        assert len(raw)==m["size"] and len(raw)<=65536 and d.eof and not d.unused_data and not d.unconsumed_tail
        assert hashlib.sha256(raw).hexdigest()==m["sha256"]
        raw.decode("utf-8"); result[name]=raw
    assert set(result)==set(names)
    return result
archives={f:json.loads((HERE/f).read_text()) for f,_ in ARCHIVES}; restored={}
for f,names in ARCHIVES: restored.update(restore(archives[f],names))
assert set(restored)==set(NAMES)
# Value/structure controls: a wrong declared digest, path member,
# wrong size and a valid zlib stream followed by trailing bytes must fail.
controls={}
for f,names in ARCHIVES:
    for kind in ("digest","name","size","trailing"):
        bad=json.loads(json.dumps(archives[f])); m=bad["members"][0]
        if kind=="digest": m["sha256"]="0"*64
        elif kind=="name": m["name"]="../escape.py"
        elif kind=="size": m["size"]+=1
        else: m["data"]=base64.b85encode(base64.b85decode(m["data"])+b"trailing").decode()
        try: restore(bad,names)
        except (AssertionError,ValueError,zlib.error): controls[f+":"+kind]=True
        else: raise AssertionError("codec control did not fail: "+f+":"+kind)
cfg=json.loads(sys.stdin.read()or"{}")
requested=cfg.pop("source","none")
assert requested=="none" or requested in NAMES
if requested!="none":
    print(json.dumps(dict(source=requested,sha256=hashlib.sha256(restored[requested]).hexdigest(),text=restored[requested].decode(),scope="Restored source bytes only; no experiment run.")))
else:
    with tempfile.TemporaryDirectory(prefix="flat64-runtime-") as directory:
        temp=pathlib.Path(directory)
        for f in HERE.iterdir():
            assert not f.is_symlink()
            if f.is_file() and f.name not in dict(ARCHIVES): shutil.copyfile(f,temp/f.name)
        for name,data in restored.items():
            f=temp/name; f.parent.mkdir(exist_ok=True); f.write_bytes(data)
        files=[f for f in temp.rglob("*") if f.is_file()]
        assert len(files)<=32 and all(f.stat().st_size<=65536 for f in files)
        sys.stdin=io.StringIO(json.dumps(cfg)); output=io.StringIO()
        with contextlib.redirect_stdout(output): runpy.run_path(str(temp/"experiment.py"),run_name="__main__")
        report=json.loads(output.getvalue())
        report["codec"]=dict(restored_members=len(restored),restored_sha256={k:hashlib.sha256(v).hexdigest()for k,v in restored.items()},negative_controls=controls)
        print(json.dumps(report,sort_keys=True))
