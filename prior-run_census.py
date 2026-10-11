from pathlib import Path
import hashlib,json,os,subprocess,time,sys
r=Path.cwd();w=r/'out/native-ordinary-shader-info'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def win(p):return r'\\wsl.localhost\Ubuntu'+str(p.resolve()).replace('/','\\')
binding=json.loads((w/'preparation.json').read_text())
assert not subprocess.check_output(['git','status','--porcelain'])
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==binding['revision']
for n,h in binding['source_bindings'].items():assert sha(r/n)==h
inputs={str(p.relative_to(w)):sha(p) for p in w.rglob('*') if p.is_file()}
args=['/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe','-NoProfile','-NonInteractive','-ExecutionPolicy','Bypass','-File',win(w/'owner.ps1'),'-Stage',win(w),'-Output',win(w/'run'),'-Filter','native-extension-census','-Expected','native-extension-census','-Seconds','30','-TemporaryDirectory',win(w/'owner-temp')]
record={'revision':binding['revision'],'arguments':args,'input_bindings':inputs,'linux_owner_pid':os.getpid(),'linux_owner_start_ticks':Path('/proc/self/stat').read_text().split(') ',1)[1].split()[19],'scope':'Read-only native extension census; no full coupled or frame performance claim.'}
(w/'bridge.json').write_text(json.dumps(record,indent=2)+'\n')
clock=time.monotonic()
with (w/'bridge.stdout').open('wb') as out,(w/'bridge.stderr').open('wb') as err:
 p=subprocess.run(args,stdout=out,stderr=err)
record.update(returncode=p.returncode,elapsed_seconds=time.monotonic()-clock,completed=True)
for n,h in inputs.items():assert sha(w/n)==h
for n,h in binding['source_bindings'].items():assert sha(r/n)==h
assert not subprocess.check_output(['git','status','--porcelain'])
(w/'bridge.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'returncode':p.returncode,'elapsed_seconds':record['elapsed_seconds']}),flush=True)
sys.exit(p.returncode)
