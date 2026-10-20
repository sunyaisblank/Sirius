from pathlib import Path
import hashlib,json,re,shutil,zipfile
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/portable-checked-fp64-jit';B.mkdir()
bundle=R/'attestations/qualification/0109ab46/int64-admission-and-divide64-rejected'
manifest=json.loads((bundle/'MANIFEST.json').read_text())
with zipfile.ZipFile(bundle/'diagnostics-and-trial.zip') as z:
    for name in ['capture_jit.py','capture_features.py','read_features.py','feature-layout.json','probe-source.cpp']:
        (W/name).write_bytes(z.read('out/portable-divide64-jit/'+name))
    source=z.read('out/portable-divide64-jit/run_capture.py').decode().replace('bin/linux-gcc/portable-divide64-admission/candidate-endpoint-normal.spv','bin/linux-gcc/portable-checked-fp64-modules/checked-endpoint.spv').replace('bin/linux-gcc/portable-divide64-jit/short_input_jit_probe','bin/linux-gcc/portable-checked-fp64-jit/short_input_jit_probe')
    source+='\nassert json.loads((folder/"feature-verification.json").read_text())["features"]["shaderFloat64"] == 1\n'
    (W/'run_capture.py').write_text(source)
    old_inputs=json.loads(z.read('out/portable-divide64-jit/capture-baseline/inputs-before.json'))
    (W/'baseline-source-scope.json').write_bytes(z.read('out/portable-divide64-jit/capture-baseline/owner.json'))
    (W/'baseline-ordinary-divide-review.json').write_bytes(z.read('out/portable-divide64-jit/ordinary-divide-review.json'))
member='bin/linux-gcc/portable-divide64-jit/short_input_jit_probe'
with zipfile.ZipFile(bundle/'actual-products.zip') as z:
    raw=z.read(member);assert hashlib.sha256(raw).hexdigest()==manifest['files'][member]['sha256']
    (B/'short_input_jit_probe').write_bytes(raw);(B/'short_input_jit_probe').chmod(0o755)
for p,old in old_inputs.items():
    if p.startswith(('/usr/','/etc/')):
        path=Path(p);assert hashlib.sha256(path.read_bytes()).hexdigest()==old['sha256'],p
for name in ['owned_guard.py','run_owned.py']:shutil.copyfile(R/'out/portable-checked-fp64'/name,W/name)
(W/'reuse.json').write_text(json.dumps({'probe_byte_exact_original':manifest['files'][member],'original_probe_CXX_revision':'e0675ee973772aac04178f4a7cae3c2045e8cac2','current_public_logic_equal_except_LoadKernel_call_formatting':True,'passive_GDB_and_feature_scripts_unchanged':True,'observer_binary_recompiled':False,'baseline_capture_reused_original_actual_module_and_same_probe_factory_driver_LLVM_loader_observer':True,'all_saved_system_inputs_equal':True,'baseline_executed_again':False,'scope':'Ordinary complete module/code binding; short input guard does not execute retained numerical arithmetic.'},indent=2)+'\n')
