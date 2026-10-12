"""Supplement missing generated-header identities without rewriting receipts."""
from pathlib import Path
import ast
import hashlib
import importlib.util
import json
import re
import struct
import subprocess
import sys

sys.dont_write_bytecode = True
root = Path.cwd()
work = Path(__file__).resolve().parent
private = root / 'bin/linux-gcc/portable-denominator-preparation'
folder = Path(sys.argv[1])
identity = lambda p: {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

builder = load('retained_builder', root / 'scripts/build-retained-kernels.py')
program = load('retained_program', root / 'src/sirius/kernels/retained_program.py')
generated = folder / 'generated'
generated.mkdir()
builder.generate_owner_layout(program, generated)
layout = private / 'products/retained_owner_layout.h'
assert layout.read_bytes() == (generated / layout.name).read_bytes()
previous = root / 'attestations/software-vulkan/1d59f17/portable-exact-cancellation/files/bin/linux-gcc/portable-exact-cancellation/products'
assert layout.read_bytes() == (previous / layout.name).read_bytes()

lines = ['// Generated optional portable owner shaders; do not edit.', '#pragma once',
         '#include <array>', '#include <cstdint>', 'namespace sirius::backend::retained_owner {']
arrays = {}
for kind in ('Transport', 'Endpoint'):
    module = private / ('products/retained_' + kind.lower() + '_portable_stage7_owner.spv')
    data = module.read_bytes()
    words = struct.unpack('<%dI' % (len(data) // 4), data)
    name = 'k' + kind + 'Shader'
    arrays[name] = data
    lines.append(f'inline constexpr std::array<std::uint32_t, {len(words)}> {name}{{{{')
    for start in range(0, len(words), 16):
        lines.append(','.join(str(v) + 'u' for v in words[start:start + 16]) + ',')
    lines.append('}};')
lines.append('} // namespace sirius::backend::retained_owner')
owner_header = private / 'products/retained_owner_kernels.h'
assert owner_header.read_bytes() == ('\n'.join(lines) + '\n').encode()

# Read the actual optional array spans in the already-linked executable.
tree = ast.parse((work / 'verify_consumers.py').read_text())
function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'symbols')
source = ast.get_source_segment((work / 'verify_consumers.py').read_text(), function)
source = source.replace('16retained_program', '14retained_owner')
namespace = {'struct': struct, 're': re}
exec(source, namespace)
exe = (private / 'original_candidate_controls').read_bytes()
symbols = namespace['symbols'](exe)
for name, payload in arrays.items():
    _, offset, size = symbols[name]
    assert size == len(payload) and exe[offset:offset + size] == payload

fixture = root / 'bin/linux-gcc/tests/backend/retained_camera/program_fixture.h'
subprocess.run([sys.executable, '-B', str(root / 'tests/support/retained_camera/build_fixture.py'),
                '--output', str(generated / fixture.name)], check=True)
assert fixture.read_bytes() == (generated / fixture.name).read_bytes()
headers = {str(p.relative_to(root)): identity(p) for p in (owner_header, layout, fixture)}
report = {'passed': True, 'headers': headers, 'owner_arrays_full_elf_joins': 2,
          'layout_exact_from_sealed_program_source_and_archived_map': True,
          'fixture_exact_from_sealed_program_and_reference_cases': True,
          'historical_compiler_seals_omitted_these_three_headers': True,
          'historical_receipts_rewritten': False,
          'scope': 'Current producer/reference and executable joins plus future comparison freeze. '
                   'Does not assert retroactive pre/post compiler observation of omitted headers.'}
(work / 'generated-header-provenance.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report))
