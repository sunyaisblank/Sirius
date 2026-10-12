"""Use the existing owned guard for this batch's export/stage or offline reader."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import hashlib,json,shutil,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
assert len(sys.argv)==2 and sys.argv[1] in ('download','stage','device-timestamps-analysis','mixed-transport-analysis','factored-endpoints-analysis')
action=sys.argv[1];source=json.loads((WORK/'source.json').read_text());head=source['source_revision']
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==head
assert subprocess.check_output(['git','rev-parse','@{u}'],text=True).strip()==head
assert not subprocess.check_output(['git','status','--porcelain'],text=True)
ci=json.loads((WORK/'ci.json').read_text());assert ci['accepted'] and ci['source_revision']==head
folder=WORK/(action+'-owner');assert not folder.exists();folder.mkdir()
for name,digest in {'run_owned.py':'b47e4bb22e63c16244802a5d121160a63b942f5c37614d26ded319fe4aaea57d','owned_guard.py':'747fe1e0bcaa0cad8fb867f22a30935b1c367dc0e2c0ad8204a87f6593899056'}.items():
    path=WORK/name;assert hashlib.sha256(path.read_bytes()).hexdigest()==digest;shutil.copyfile(path,folder/name)
paths=[p for p in WORK.iterdir() if p.is_file()]+list((WORK/'kernels').iterdir())
if action=='download':
    argv=['python3','-B',str(WORK.relative_to(ROOT)/'download_export.py')];limits={'seconds':420,'rss_kib':4194304,'reserve_kib':2097152}
    scope='Existing420s4GiB2GiBreserve exact official native-build export preservation only; no native/scientific dispatch or qualification verdict.'
elif action=='stage':
    assert json.loads((WORK/'artifact-extraction.json').read_text())['accepted']
    paths += [p for p in (WORK/'artifact').rglob('*') if p.is_file()]
    argv=['python3','-B',str(WORK.relative_to(ROOT)/'prepare_native_stage.py')];limits={'seconds':420,'rss_kib':4194304,'reserve_kib':2097152}
    scope='Existing420s4GiB2GiBreserve authenticated native export reconstitution and exact executable/module/provider/source/resource joins. No native invocation/scientific qualification.'
else:
    case=action.removesuffix('-analysis');paths += [p for p in (WORK/case).rglob('*') if p.is_file()]
    paths += [WORK/'baseline-products/bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h']
    argv=['python3','-B',str(WORK.relative_to(ROOT)/'analyze_words.py'),case];limits={'seconds':60,'rss_kib':524288,'reserve_kib':2097152}
    scope='Existing60s512MiB2GiBreserve offline finite original-control/capture/actual-word reader only. No dispatch or universal/ordinary-performance/fullqualification verdict.'
settings={action:{'revision':head,'kind':'diagnostic','expected_exit':0,'argv':argv,'limits':limits,'frozen_inputs':list(dict.fromkeys(str(p.relative_to(ROOT)) for p in paths if p.is_file())),'scope':scope}}
(folder/'actions.json').write_text(json.dumps(settings,indent=2)+'\n');print(folder)
