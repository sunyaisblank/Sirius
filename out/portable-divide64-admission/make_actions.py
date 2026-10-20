from pathlib import Path
import json,re,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent
inputs={str(p.relative_to(R)) for p in W.rglob('*') if p.is_file()}
inputs|={'bin/linux-gcc/src/sirius/backend/retained/retained_endpoint_portable_normal_sum.spv','/usr/bin/python3'}
for tool in json.loads((W/'endpoint-arguments.json').read_text())['tools'].values():
 p=Path(tool);assert p.is_file();inputs.add(str(p))
 result=subprocess.run(['/usr/bin/ldd',str(p)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
 if result.returncode==0:inputs|={str(Path(a).resolve(strict=True)) for a in re.findall(r'(/[^\s()]+)',result.stdout) if Path(a).is_file()}
a={'revision':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'kind':'diagnostic','expected_exit':0,
 'limits':{'seconds':120,'rss_kib':4194304,'reserve_kib':2097152},'argv':['/usr/bin/python3','-B',str(W/'compile.py')],
 'frozen_inputs':sorted(inputs),'software_provider':False,'scope':'One fresh complete current Endpoint candidate compile after repairing the old candidate-only capability-order guard. Exact current stock module/assembly reused. All original float/interface/mode/barrier/divide checks retained; no Vulkan call, arithmetic, JIT or benefit claim.'}
(W/'actions.json').write_text(json.dumps({'compile':a},indent=2)+'\n');print(json.dumps({'frozen_inputs':len(inputs),'fresh_compilations':1,'stock_compilation_reused':True}))
