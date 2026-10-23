"""Seal the positive-frame reader's existing source/object/scene authorities."""
from pathlib import Path
import json
R=Path.cwd();W=Path(__file__).resolve().parent
assert not (W/'frame').exists()
p=W/'actions.json';actions=json.loads(p.read_text());a=actions['frame']
extra={'bin/linux-gcc/compile_commands.json',
 'bin/linux-gcc/tests/render/CMakeFiles/sirius_render_tests.dir/app_render_execution_test.cpp.o',
 'out/native-original-input-ownership/cold-frame/result.json',
 'attestations/qualification/08add7e5/native-original-input-ownership-finite-outcome/findings.json',
 'out/lavapipe-owner-bulk-zero/cold-frame/build/owner.json',
 'out/lavapipe-owner-bulk-zero/cold-frame/build/inputs-before.json',
 'out/lavapipe-owner-bulk-zero/cold-frame/build/inputs-after.json',
 'out/lavapipe-owner-bulk-zero/readers/control-8-acceptance.json',
 str(Path(__file__).resolve().relative_to(R))}
assert (W/'read_frame.py').is_file()
extra.add(str((W/'read_frame.py').relative_to(R)))
for key in extra:assert (R/key).is_file(),key
a['frozen_inputs']=sorted(set(a['frozen_inputs'])|extra)
a['scope']+=' Existing compile commands, actual required-runtime test object, original scene record and positive-frame reader explicitly sealed. No render argument, limit, governor or assertion changed.'
p.write_text(json.dumps(actions,indent=2)+'\n')
print(json.dumps({'original_frame_unchanged':True,'reader_authorities_explicitly_sealed':sorted(extra)}))
