"""Relink two private providers from one build; replace exactly one archive member."""
from pathlib import Path
import hashlib, json, os, subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
S = R / 'bin/linux-gcc/lavapipe-vector-staged'
selection = json.loads((R / 'out/lavapipe-vector-candidate/baseline-link-selection.json').read_text())
B = Path(selection['cwd'])
original = B / selection['llvmpipe_archive']
candidate = S / 'candidate-object/lp_state_cs.c.o'
def identity(path):
    return {'path':str(path), 'bytes':path.stat().st_size,
            'sha256':hashlib.file_digest(path.open('rb'),'sha256').hexdigest()}
def run(argv, **kw):
    return subprocess.check_output(argv, cwd=B, **kw)

names = run(['/usr/bin/ar','t',str(original)], text=True).splitlines()
assert len(names) == len(set(Path(n).name for n in names))
members = []
for name in names:
    p = Path(name)
    if not p.is_absolute(): p = B / p
    p = p.resolve(strict=True)
    assert p.is_relative_to(B) and p.name.endswith('.c.o')
    members.append(p)
assert sum(p.name=='lp_state_cs.c.o' for p in members) == 1
original_driver = B / 'src/gallium/targets/lavapipe/libvulkan_lvp.so'
before = {str(p):identity(p) for p in [original, original_driver, *members, candidate]}
arms = {}
for arm in ['baseline','candidate']:
    out = S / arm
    out.mkdir()
    archive = out / 'libllvmpipe.a'
    inputs = [candidate if arm=='candidate' and p.name=='lp_state_cs.c.o' else p for p in members]
    argv = ['/usr/bin/ar','crsD',str(archive),*[str(p) for p in inputs]]
    subprocess.run(argv,cwd=B,check=True)
    assert archive.read_bytes()[:8] == b'!<arch>\n'
    listed = run(['/usr/bin/ar','t',str(archive)],text=True).splitlines()
    assert listed == [p.name for p in members]
    checked=[]
    for name,p in zip(listed,inputs):
        raw=run(['/usr/bin/ar','p',str(archive),name])
        assert raw == p.read_bytes(),name
        checked.append({'name':name,**identity(p)})
    target = out / 'libvulkan_lvp.so'
    link = list(selection['argv'])
    link[link.index('-o')+1]=str(target)
    assert link.count(selection['llvmpipe_archive'])==1
    link[link.index(selection['llvmpipe_archive'])]=str(archive)
    subprocess.run(link,cwd=B,check=True)
    icd=out/'lvp_icd.json'
    icd.write_text(json.dumps({'file_format_version':'1.0.0','ICD':{
        'library_path':str(target),'api_version':'1.4.305'}},indent=2)+'\n')
    arms[arm]={'archive':identity(archive),'members':checked,'link_argv':link,
               'driver':identity(target),'ICD':identity(icd)}
assert before == {str(p):identity(p) for p in [original, original_driver, *members, candidate]}
assert (S/'baseline/libvulkan_lvp.so').read_bytes()==original_driver.read_bytes(), 'matched baseline relink differs from original baseline output'
diff=[a['name'] for a,b in zip(arms['baseline']['members'],arms['candidate']['members']) if a['sha256']!=b['sha256']]
assert diff==['lp_state_cs.c.o'],diff
result={'passed':True,'original_archive':before[str(original)],'original_driver':before[str(original_driver)],
        'regular_archive_member_order_and_bytes_verified':True,'changed_members':diff,
        'baseline_relink_byte_equal_to_original':True,'arms':arms,
        'scope':'Private matched Mesa25.2.8/GCC14/LLVM21 providers; no installation, dispatch, arithmetic, benefit or release qualification.'}
(W/'matched-links.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['passed','changed_members','baseline_relink_byte_equal_to_original','scope']}))
