from pathlib import Path
import hashlib,json,os,subprocess,tempfile
root=Path.cwd()
pairs=[('evidence/issue-48-unforced-translator-match','attestations/native-vulkan/4bcd987/macos-unforced-translator-match'),('evidence/issue-100-kerr-png-reconstruction-pair','attestations/forensic-images/b739f92/kerr-png-reconstruction-pair'),('evidence/issue-100-kerr-exr-reconstruction-pair','attestations/forensic-images/b739f92/kerr-exr-reconstruction-pair')]
for tag,name in pairs:
 base=root/name;m=json.loads((base/'manifest.json').read_text());files=m['files'];rows=[]
 for item in subprocess.check_output(['git','ls-tree','-rz',tag]).split(b'\0'):
  if not item:continue
  metadata,path=item.split(b'\t',1);mode,typ,sha=metadata.decode().split();assert typ=='blob';p=root/os.fsdecode(path);rel=str(p.relative_to(base));expected=files.get(rel)
  if rel=='manifest.json':expected={'type':'file','sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
  assert expected is not None,rel
  if mode=='120000':assert p.is_symlink() and os.readlink(p)==expected['target'];data=os.readlink(p).encode()
  else:assert p.is_file() and not p.is_symlink();data=p.read_bytes();assert (mode=='100755')==bool(p.stat().st_mode&0o111)
  assert len(data)==expected['bytes'] and hashlib.sha256(data).hexdigest()==expected['sha256'],rel
  rows.append((sha,expected,rel))
 assert len(rows)==len(files)+1
 fd,temp=tempfile.mkstemp(prefix='git-evidence-readback-',dir=root/'out');os.close(fd)
 try:
  with open(temp,'wb') as f:subprocess.run(['git','cat-file','--batch'],input=''.join(sha+'\n' for sha,_,_ in rows).encode(),stdout=f,check=True)
  with open(temp,'rb') as f:
   for sha,expected,rel in rows:
    header=f.readline().decode().strip().split();assert header==[sha,'blob',str(expected['bytes'])],rel
    data=f.read(expected['bytes']);assert len(data)==expected['bytes'] and hashlib.sha256(data).hexdigest()==expected['sha256'] and f.read(1)==b'\n',rel
   assert f.read()==b''
 finally:Path(temp).unlink()
 print(json.dumps({'tag':tag,'raw_and_git_payloads_verified':len(rows),'manifest_sha256':hashlib.sha256((base/'manifest.json').read_bytes()).hexdigest()}),flush=True)
