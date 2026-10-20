"""One reused exact observer, three unique modules, one replaced archive member."""
from pathlib import Path
import hashlib,json,shutil,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/portable-checked-fp64-matched';C=R/'bin/linux-gcc'
def sha(p):return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
commands=json.loads((W/'commands.json').read_text())
for row in commands['compile']:subprocess.run(row['argv'],cwd=row['cwd'],check=True)
original=C/'src/sirius/backend/libsirius_backend.a';candidate=B/'candidate/libsirius_backend.a'
assert original.read_bytes()[:8]==b'!<arch>\n';before=sha(original)
shutil.copyfile(original,candidate)
subprocess.run(['/usr/bin/ar','rDs',str(candidate),str(B/'candidate/retained_compute.cpp.o')],check=True)
members=subprocess.check_output(['/usr/bin/ar','t',str(original)],text=True).splitlines()
assert members==subprocess.check_output(['/usr/bin/ar','t',str(candidate)],text=True).splitlines()
different=[]
for n in members:
 a=subprocess.check_output(['/usr/bin/ar','p',str(original),n]);b=subprocess.check_output(['/usr/bin/ar','p',str(candidate),n])
 if a!=b:different.append(n)
assert different==['retained_compute.cpp.o'] and sha(original)==before
for arm,row in commands['link'].items():subprocess.run(row['argv'],cwd=row['cwd'],check=True)
assert sha(original)==before
(W/'matched-result.json').write_text(json.dumps({'passed':True,'source_libraries_unmodified':True,'changed_archive_members':different,
 'shared_original_observer_reused_byte_exact':True,'baseline_ELF_sha256':sha(B/'baseline-tests'),'candidate_ELF_sha256':sha(B/'candidate-tests'),
 'source_scope':json.loads((W/'source-scope.json').read_text()),'scope':'Private matched original observer pair only; no scientific, integrated, benefit or release acceptance.'},indent=2)+'\n')
print(json.dumps({'passed':True,'changed_archive_members':different,'original_observer_compilations':0,'Vulkan_calls':0}))
