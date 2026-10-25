"""Keep exact object/SPIRV/FP joins; follow Rows to field control/layer exports."""
from pathlib import Path

basis = Path.cwd() / 'out/native-cpu-zero-term/read_linked_native.py'
source = basis.read_text()
before = "'sirius_simd_' + slug + '_packet'"
assert source.count(before) == 1
source = source.replace(before, "'sirius_native_' + slug + '_field_try_packet'")
before = "exec(compile(source, str(basis), 'exec'))"
assert source.count(before) == 1
source = source.replace(before, '''before = "assert not re.search(r'\\\\b(?:v?fm(?:add|sub)|v?fnm(?:add|sub))\\\\w*\\\\b', all_native_disassembly)"
assert source.count(before) == 1
after = """for slug in ('transport', 'endpoint', 'dopri_phase'):
    assert re.search(r'\\\\b(?:call|jmp)\\\\s+[0-9a-fx]+\\\\s+<sirius_simd_' + slug + '_field_layer>', all_native_disassembly)
""" + before
source = source.replace(before, after)
source = source.replace('callback_shared_mode_reaches_rows_pure_SIMD_and_scalar_all3', 'callback_reaches_rows_pure_field_control_layer_and_scalar_all3')
source = source.replace('SIMD8 and original scalar tails', 'new field control/layer exports and original D fallback')
exec(compile(source, str(basis), 'exec'))''')
exec(compile(source, str(basis), 'exec'))
