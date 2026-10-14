"""Issue 101 preparation and one compiler-only normal-Transport screen.

No mode invokes Vulkan, a pipeline, an integrated control or a benefit trial.
The separately declared capture action is for root's later explicit operation.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import re
import resource
import subprocess
import sys

sys.dont_write_bytecode = True
ROOT = Path.cwd().resolve()
WORK = Path(__file__).resolve().parent
PRODUCTS = ROOT / 'bin/linux-gcc/computed-carry-exit-screen'
REV = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'
CANONICAL = ROOT / 'bin/linux-gcc/src/sirius/backend/retained/retained_transport_portable_normal_sum.spv'
PRIOR = ROOT / 'attestations/software-vulkan/4bcd987/guarded-product-high-jit'
PROBE = PRIOR / 'files/bin/linux-gcc/guarded-product-high-runtime/short_input_jit_probe'
OLD = ROOT / 'attestations/software-vulkan/4bcd987/power-two-divide-prerequisites/files/out/power-two-divide-screen'


def record(path):
    path = Path(path)
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': path.stat().st_size, 'sha256': digest}


def write(name, value):
    (WORK / name).write_text(json.dumps(value, indent=2) + '\n')


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def verify_frozen():
    frozen = json.loads((WORK / 'source-freeze.json').read_text())
    for name, expected in frozen['private_inputs'].items():
        assert record(WORK / name) == expected, name
    for name, expected in frozen['external_inputs'].items():
        assert record(name) == expected, name
    for path in (WORK / 'stock/src/sirius/kernels').iterdir():
        relative = 'src/sirius/kernels/' + path.name
        original = subprocess.check_output(['git', 'show', REV + ':' + relative])
        assert path.read_bytes() == original == (ROOT / relative).read_bytes(), relative
    builder = 'scripts/build-retained-kernels.py'
    assert (ROOT / builder).read_bytes() == subprocess.check_output(['git', 'show', REV + ':' + builder])


def pin_context():
    """Root pins the then-clean CI/diagnostic HEAD without changing product source."""
    verify_frozen()
    assert not PRODUCTS.exists(), PRODUCTS
    assert not (WORK / 'compile-normal-transport').exists()
    assert not (WORK / 'capture-candidate-transport-normal').exists()
    assert not (WORK / 'operator-context.json').exists(), 'pin once; do not overwrite a consumed context'
    assert not subprocess.check_output(['git', 'status', '--porcelain=v1'])
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    actions = json.loads((WORK / 'actions.json').read_text())
    assert all(a['revision'] == 'PENDING_OPERATOR_CLEAN_CONTEXT' for a in actions.values())
    for action in actions.values():
        action['revision'] = revision
    write('actions.json', actions)
    write('operator-context.json', {'revision': revision, 'status': '',
        'product_revision': REV, 'actions': record(WORK / 'actions.json'),
        'source_freeze': record(WORK / 'source-freeze.json'),
        'scope': 'One clean operator context; no compiler or device child launched by pin.'})
    print('Pinned clean operator context', revision, flush=True)


def accumulator(text):
    match = re.search(r'OpName (%\S+) "RTAccumulateTerm"', text)
    assert match
    body = re.search(r'^\s*' + re.escape(match[1]) + r' = OpFunction\b.*?OpFunctionEnd', text, re.M | re.S)
    assert body
    return body.group(0) + '\n'


def compile_normal_transport():
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    _, previous = resource.getrlimit(resource.RLIMIT_FSIZE)
    file_limit = 64 * 1024**2 if previous == resource.RLIM_INFINITY else min(previous, 64 * 1024**2)
    resource.setrlimit(resource.RLIMIT_FSIZE, (file_limit, file_limit))
    verify_frozen()
    context = json.loads((WORK / 'operator-context.json').read_text())
    assert context['actions'] == record(WORK / 'actions.json')
    assert context['source_freeze'] == record(WORK / 'source-freeze.json')
    assert not PRODUCTS.exists(), PRODUCTS
    (PRODUCTS / 'modules').mkdir(parents=True)
    builder = module('issue101_production_builder', ROOT / 'scripts/build-retained-kernels.py')
    programs = module('issue101_stock_programs', WORK / 'stock/src/sirius/kernels/retained_program.py')
    program = programs.build_transport_program(parallel=True)
    special = programs.build_schwarzschild_transport_program(parallel=True)
    encode = lambda p: [p['instructions'], p['registers'], len(p['outputs'])] + p['outputs'] + p['operations'] + p['layer_offsets']
    dimensions = (len(encode(program)), len(encode(special)), len(special['layer_offsets']) - 1)
    assert special['registers'] == program['registers'] and len(special['outputs']) == 40
    cache = (ROOT / 'bin/linux-gcc/CMakeCache.txt').read_text()
    tool = lambda n: re.search(r'^' + n + r':FILEPATH=(.*)$', cache, re.M)[1]
    compiler, assembler, disassembler, validator, optimizer = [tool(n) for n in
        ('SIRIUS_SLANGC', 'SIRIUS_SPIRV_AS', 'SIRIUS_SPIRV_DIS', 'SIRIUS_SPIRV_VAL', 'SIRIUS_SPIRV_OPT')]
    texts = {}
    results = {}
    for owner in ('stock', 'candidate'):
        source = WORK / owner / 'src/sirius/kernels/retained_transport.slang'
        output = PRODUCTS / 'modules' / (owner + '-transport-normal.spv')
        builder.compile_shader(source, output, compiler, assembler, disassembler, validator,
            program['registers'], 5, len(program['layer_offsets']) - 1,
            program.get('prefix_instructions', 0), portable=True, optimizer=optimizer,
            normal_sum32=True, transport_specialized=dimensions)
        assembly = output.with_suffix('.final.spvasm')
        subprocess.run([disassembler, str(output), '-o', str(assembly)], check=True)
        texts[owner] = assembly.read_text()
        body = accumulator(texts[owner])
        (PRODUCTS / 'modules' / (owner + '-RTAccumulateTerm.spvasm')).write_text(body)
        results[owner] = {'module': record(output), 'assembly': record(assembly),
            'accumulator_loop_count': body.count('OpLoopMerge'),
            'accumulator_sum_call_count': len(re.findall(r'OpFunctionCall .* %RetainedSum32\S* ', body)),
            'accumulator_radius_update_count': len(re.findall(r'OpFunctionCall .* %RPUpAdd ', body)),
            'accumulator_bound_call_count': len(re.findall(r'OpFunctionCall .* %RPBounded ', body))}
    stock = PRODUCTS / 'modules/stock-transport-normal.spv'
    candidate = PRODUCTS / 'modules/candidate-transport-normal.spv'
    assert stock.read_bytes() == CANONICAL.read_bytes()
    assert candidate.read_bytes() != stock.read_bytes()
    inventory = lambda t: sorted(re.findall(r'OpExecutionMode %\S+ (.*)|OpDecorate %\S+ ((?:Binding|DescriptorSet|ArrayStride) \d+)|OpMemberDecorate %\S+ (\d+ Offset \d+)', t))
    assert inventory(texts['stock']) == inventory(texts['candidate'])
    for instruction in ('OpControlBarrier', 'OpMemoryBarrier'):
        assert texts['stock'].count(instruction) == texts['candidate'].count(instruction)
    for owner in results:
        assert results[owner]['accumulator_loop_count'] == 1
        assert results[owner]['accumulator_sum_call_count'] == 1
        assert results[owner]['accumulator_radius_update_count'] == 1
        assert results[owner]['accumulator_bound_call_count'] == 5
    verify_frozen()
    write('compiler-result.json', {'production_policy_passed': True,
        'stock_equals_current_canonical': True, 'modules': results,
        'mode_binding_stride_member_offset_inventory_exact': True,
        'barrier_counts_exact': True, 'manual_accumulator_control_flow_review_required': True,
        'scope': 'Compiler-only complete normal Transport modules. Inventory checks support ABI review; they do not establish semantics, ordinary JIT behavior or benefit.'})
    print('Complete normal Transport pair passes original compiler policy; inspect both accumulator CFGs.', flush=True)


if __name__ == '__main__':
    {'pin-context': pin_context, 'compile': compile_normal_transport}[sys.argv[1]]()
