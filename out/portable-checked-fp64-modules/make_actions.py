from pathlib import Path
import json,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent
inputs={str(p.relative_to(R)) for p in W.rglob('*') if p.is_file()}|{'out/portable-checked-fp64/endpoint-arguments.json','out/portable-checked-fp64/run/complete-output-agreement.json','out/portable-checked-fp64/host/owner.json','bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'}
tools=json.loads((R/'out/portable-checked-fp64/endpoint-arguments.json').read_text())['tools'];inputs|=set(tools.values())
for stem in ['endpoint','transport','dopri_phase']:inputs.add('bin/linux-gcc/src/sirius/backend/retained/retained_'+stem+'_portable_normal_sum.spv')
revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
(W/'actions.json').write_text(json.dumps({'compile':{'revision':revision,'kind':'diagnostic','expected_exit':0,'limits':{'seconds':120,'rss_kib':4194304,'reserve_kib':2097152},'argv':['/usr/bin/python3','-B',str(W/'compile.py')],'frozen_inputs':sorted(inputs),'software_provider':False,'scope':'Compile complete private Endpoint/Transport/DopriPhase normal modules with checked FP64 estimates. All program words, existing binary32 contracts, ABI/layout/modes/barriers and unchanged integer fallbacks required. No runtime/benefit qualification.'}},indent=2)+'\n')
