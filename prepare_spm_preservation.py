"""Adapt the already reviewed finished-work preservation and cleanup controllers."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import ast,json,subprocess
ROOT=Path.cwd();WORK=Path(__file__).resolve().parent
TAG='evidence/issue-44-native-executed-cost-6be3132'
def blob(n):return subprocess.check_output(['git','show',TAG+':'+n]).decode()
def replace(s,o,n):
 assert s.count(o)==1,(o[:100],s.count(o))
 return s.replace(o,n)
outputs={}
def save(n,s):assert n not in outputs and not (WORK/n).exists();outputs[n]=s
for n in ('collect_finished_births.py','run_finished_census.py'):save(n,blob(n))
ps=blob('fresh_finished_consumers.ps1')
ps=replace(ps,"$targets=@('native-executed-cost',","$targets=@('native-spm-cost',")
save('fresh_finished_consumers.ps1',ps)
protect=blob('protect_executed_cost.py')
protect=replace(protect,"tag = 'evidence/issue-44-native-executed-cost-'", "tag = 'evidence/issue-44-native-spm-'")
protect=protect.replace("/ 'native-executed-cost'","/ 'native-spm-cost'")
start=protect.index('    paths = sorted(p for p in WORK.iterdir()')
end=protect.index('    index = WORK /',start)
protect=protect[:start]+'''    paths = sorted(p for p in WORK.iterdir() if p.is_file() and
                   (p.suffix in ('.py','.ps1','.cs','.cpp','.h','.md','.txt','.json') or
                    p.name in ('libamdrdf-license',)))
'''+protect[end:]
protect=replace(protect,'Preserve reviewed original native profiler acquisition and analysis at ',
                'Preserve reviewed original native SPM availability acquisition at ')
protect=replace(protect,"prior_path = ROOT / 'attestations/native-vulkan/21e55ce/guarded-native-default-fma-rejected/manifest.json'",
                "prior_path = ROOT / 'attestations/native-vulkan/6be3132/native-executed-cost/manifest.json'")
protect=replace(protect,"for value in [*prior['executed_input_ledger'], *prior['diagnostics'],\n              *prior['volatile_product_and_native_stage_files']]:",
                "for value in [*prior['executed_input_ledger'], *prior['records']]:")
protect=replace(protect,"for value in source['files']:","recoveries[key(seal(prior_path))] = {'kind':'retained-prior-payload','path':str(prior_path.relative_to(ROOT))}\nfor value in source['files']:")
protect=replace(protect,"download = read(WORK / 'tools-download.json')\narchive = ROOT / download['archive']\narchive_value = expected(download)",
                "download = prior['tool_archive']\nassert download['recovery']['kind'] == 'local-archive-payload'\narchive = ROOT / download['recovery']['path']\narchive_value = expected(download)")
protect=replace(protect,'Whole finished acquisition/GUI diagnostic files and exact consumed input versions preserved. Prior failed epochs remain distinct. Source-only snapshot is pushed; raw archives and official tool ZIP remain local. No rendering/scientific/release completion or cleanup clearance.',
                'Whole finished SPM availability diagnostic files and exact consumed input versions preserved. Prior SPM-disabled traces remain immutable and separately classified. Useful source is pushed; raw trace/archive and existing official tool ZIP remain local. No event/source cost, rendering/scientific/release completion or cleanup clearance.')
save('protect_spm_cost.py',protect)
finish=blob('finish_executed_cost.py')
finish=finish.replace("'out/native-executed-cost'","'out/native-spm-cost'")
finish=finish.replace("/ 'native-executed-cost'","/ 'native-spm-cost'")
finish=replace(finish,"assert plan['remove'] == ['out/native-spm-cost', zip_origin]", "assert plan['remove'] == ['out/native-spm-cost']")
finish=replace(finish,"assert zip_origin == 'bin/windows-msvc/qualification-tools/rdts-2026-10-02-1872/RadeonDeveloperToolSuite-2026-10-02-1872.zip'",
               "assert zip_origin == 'attestations/native-vulkan/6be3132/native-executed-cost/payloads/128fdebb889d1fda1d04ba9e4792fbef1b18c086b4e901642262653a7a53917e'")
finish=replace(finish,"for p in [*[ROOT / name for name,v in original_tree.items() if isinstance(v,dict)], ROOT / zip_origin]:",
               "for p in [ROOT / name for name,v in original_tree.items() if isinstance(v,dict)]:")
finish=replace(finish,"('native-executed-cost', zip_origin)","('native-spm-cost', 'native-b81e061-executed-cost', 'rdts-2026-10-02-1872')")
finish=replace(finish,'Only finished scratch and redundant archived ZIP removed.',
               'Only finished SPM diagnostic scratch removed; existing protected ZIP preserved.')
finish=replace(finish,"(ROOT / zip_origin).unlink(); shutil.rmtree(WORK)","shutil.rmtree(WORK)")
finish=replace(finish,"assert not WORK.exists() and not (ROOT / zip_origin).exists()", "assert not WORK.exists() and seal(ROOT / zip_origin) == expected(manifest['tool_archive'])")
save('finish_spm_cost.py',finish)
for n,s in outputs.items():
 if n.endswith('.py'):ast.parse(s,filename=n)
for n,s in outputs.items():(WORK/n).write_text(s)
print(json.dumps({'prepared_preservation_sources':len(outputs)}))
