from pathlib import Path
import shutil,zipfile
R=Path.cwd();W=Path(__file__).resolve().parent
shutil.copytree(R/'src/sirius/kernels',W/'kernels')
shutil.copyfile(R/'out/portable-checked-fp64/portable_binary32.h',W/'kernels/portable_binary32.h')
s=(R/'scripts/build-retained-kernels.py').read_text()
s=s.replace('def portable_normal_sum_controls(', 'def stock_portable_normal_sum_controls(',1)
needle='        if capabilities != ["Shader"] or not integer_widths or set(integer_widths) != {"32"}:'
assert s.count(needle)==1
s=s.replace(needle,'        if sorted(capabilities) != ["Float64", "Shader"] or not integer_widths or set(integer_widths) != {"32"}:')
needle='def compile_shader('
assert s.count(needle)==1
s=s.replace(needle,(W/'checked_controls.py').read_text()+'\n\n'+needle)
(W/'candidate_builder.py').write_text(s)
with zipfile.ZipFile('attestations/qualification/0109ab46/int64-admission-and-divide64-rejected/diagnostics-and-trial.zip') as z:
    (W/'interface_reader.py').write_bytes(z.read('out/portable-divide64-admission/historical-screen.py'))
for n in ['owned_guard.py','run_owned.py']:
    shutil.copyfile(R/'out/portable-checked-fp64'/n,W/n)
