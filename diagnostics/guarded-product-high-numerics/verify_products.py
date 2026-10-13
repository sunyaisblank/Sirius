"""Compare actual final SPIR-V operands, typed constants, call graphs and wires."""
from pathlib import Path
import collections
import hashlib
import json
import re
import sys

sys.dont_write_bytecode = True
root = Path.cwd().resolve()
work = Path(__file__).resolve().parent
stage = root / 'bin/linux-gcc/guarded-product-high-numerics'
pilot = root / 'attestations/software-vulkan/4bcd987/guarded-product-high-screen/files/out/guarded-product-high-screen/modules/candidate-endpoint-core-final.spvasm'


class Module:
    def __init__(self, text):
        self.text = text
        self.defs = {ident: rhs.split() for ident, rhs in
                     re.findall(r'^\s*(%\S+) = (Op\S+[^\n]*)$', text, re.M)}
        self.functions = {re.match(r'(%\S+)', body)[1]: body for body in
                          re.findall(r'%\S+ = OpFunction .*?OpFunctionEnd', text, re.S)}
        self.names = {ident: name for ident, name in re.findall(r'OpName (%\S+) "([^"]+)"', text)}
        self.decorations = collections.defaultdict(list)
        self.members = collections.defaultdict(list)
        for ident, value in re.findall(r'^\s*OpDecorate (%\S+) (.*)$', text, re.M):
            self.decorations[ident].append(value)
        for ident, value in re.findall(r'^\s*OpMemberDecorate (%\S+) (.*)$', text, re.M):
            self.members[ident].append(value)

    def name(self, ident):
        return re.sub(r'_\d+$', '', self.names.get(ident, ident.removeprefix('%')))

    def function(self, name):
        found = [ident for ident in self.functions if self.name(ident) == name]
        assert len(found) == 1, (name, found)
        return found[0]

    def value(self, ident):
        row = self.defs[ident]
        assert row[0].startswith(('OpType', 'OpConstant', 'OpSpecConstant'))
        return [row[0], *[self.value(token) if token.startswith('%') else token for token in row[1:]],
                sorted(self.decorations.get(ident, [])), sorted(self.members.get(ident, []))]

    def body(self, ident):
        locals_ = {}
        def resolve(match):
            token = match[0]
            if token == ident:
                return 'SELF'
            if token in self.defs and self.defs[token][0].startswith(('OpType', 'OpConstant', 'OpSpecConstant')):
                return json.dumps(self.value(token), separators=(',', ':'))
            if token in self.functions:
                return json.dumps(['CALLEE', self.name(token), self.value(self.defs[token][3])], separators=(',', ':'))
            if token not in locals_:
                locals_[token] = len(locals_)
            return 'LOCAL' + str(locals_[token])
        return '\n'.join(' '.join(re.sub(r'%[^\s]+', resolve, line).split())
            for line in self.functions[ident].splitlines()
            if not re.match(r'^\s*Op(?:Line|NoLine)\b', line))

    def closure(self, ident):
        found = set()
        def visit(current):
            if current in found:
                return
            found.add(current)
            for target in re.findall(r'= OpFunctionCall %\S+ (%\S+)', self.functions[current]):
                visit(target)
        visit(ident)
        resolved = {}
        for current in found:
            key = json.dumps([self.name(current), self.value(self.defs[current][3])], separators=(',', ':'))
            body = self.body(current)
            # Slang can emit repeated typed helper specializations. Merge only
            # when their entire resolved bodies are identical, including calls.
            assert key not in resolved or resolved[key] == body, key
            resolved[key] = body
        return sorted(resolved.items())

    def wire(self):
        bindings = {}
        for ident, decorations in self.decorations.items():
            slots = [int(value.split()[1]) for value in decorations if value.startswith('Binding ')]
            if not slots:
                continue
            assert len(slots) == 1
            row = self.defs[ident]
            assert row[0] == 'OpVariable' and row[2] == 'StorageBuffer'
            bindings[slots[0]] = [sorted(decorations), self.value(row[1]), row[2:]]
        assert set(bindings) == {0, 1}
        return {'bindings': bindings,
            'local_size': re.findall(r'OpExecutionMode %\S+ LocalSize (\d+ \d+ \d+)', self.text),
            'memory_model': re.findall(r'OpMemoryModel (.*)', self.text),
            'entry_signature': self.value(self.defs[self.function('ComputeMain')][3])}


def normalized_module(text):
    ids = {}
    def replace(match):
        token = match[0]
        if token not in ids:
            ids[token] = len(ids)
        return '%' + str(ids[token])
    return '\n'.join(' '.join(re.sub(r'%[^\s]+', replace, line).split())
        for line in text.splitlines() if line.strip() and not line.lstrip().startswith(';'))


pilot_module = Module(pilot.read_text())
guarded = pilot_module.closure(pilot_module.function('PB32ProductResidual'))
reports = []
variants = [('residual', 'core'), ('residual', 'rte32')]
variants += [('limb', mode) for mode in ['integer', 'normal', 'normal_transport', 'normal_endpoint']]
for family, mode in variants:
    baseline_mode = 'core' if family == 'residual' else mode
    stock_path = stage / f'{family}-stock-{baseline_mode}.final.spvasm'
    candidate_path = stage / f'{family}-candidate-{mode}.final.spvasm'
    stock, candidate = [Module(path.read_text()) for path in [stock_path, candidate_path]]
    assert stock.wire() == candidate.wire()
    assert stock.body(stock.function('PB32ProductResidual')) == candidate.body(candidate.function('PB32ProductResidualFallback'))
    assert candidate.closure(candidate.function('PB32ProductResidual')) == guarded
    old_names = ['RetainedSum32', 'RPUpAdd']
    old_stock = sorted(stock.body(ident) for ident in stock.functions if stock.name(ident) in old_names)
    old_candidate = sorted(candidate.body(ident) for ident in candidate.functions if candidate.name(ident) in old_names)
    assert old_stock == old_candidate
    old_floating = [ident for ident in stock.functions if stock.name(ident) in old_names
                    and re.search(r'= OpF(?!unction)', stock.functions[ident])]
    float_types = re.findall(r'(%\S+) = OpTypeFloat 32', candidate.text)
    assert len(float_types) == 1
    assert not re.search(r'= OpExtInst ' + re.escape(float_types[0]) + r' ', candidate.text)
    requires_rte32 = mode not in ['core', 'integer']
    assert re.findall(r'OpCapability (\S+)', candidate.text) == (['Shader', 'RoundingModeRTE'] if requires_rte32 else ['Shader'])
    assert candidate.wire()['local_size'] == ['64 1 1']
    for owner in ['stock', 'candidate'] if not (family == 'residual' and mode == 'rte32') else ['candidate']:
        folder = work / f'compile-{family}-{owner}-{mode}'
        receipt = json.loads((folder / 'owner.json').read_text())
        assert receipt['passed'] and receipt['status'] == 'terminal' and not receipt['cleanup_errors']
        assert receipt['whole_inputs_unchanged'] and receipt['source_unchanged'] and not receipt['remaining_owned_processes']
    baseline_relation = None
    if family == 'limb':
        archived = root / 'attestations/software-vulkan/661e157/portable-limb-accumulation/files/bin/linux-gcc/portable-execution-attribution/limb-screen' / f'arithmetic-baseline-{mode}.spv'
        product = stage / f'limb-stock-{mode}.spv'
        if mode == 'integer':
            assert (stage / f'limb-stock-{mode}.optimized.spv').read_bytes() == archived.read_bytes()
            assert normalized_module((stage / f'limb-stock-{mode}.spvasm').read_text()) == normalized_module(stock.text)
            baseline_relation = 'Original optimized bytes exact; complete final assembler roundtrip instructions exact modulo one-to-one IDs.'
        else:
            assert product.read_bytes() == archived.read_bytes()
            baseline_relation = 'Actual final stock module byte-exact to original baseline producer.'
    reports.append({'family': family, 'mode': mode, 'wire': candidate.wire(),
        'guarded_full_call_closure_equal_to_preserved_pilot': True,
        'stock_residual_equal_to_candidate_fallback': True,
        'old_sum_and_upward_bodies_typed_constants_callees_equal': len(old_stock),
        'old_normal_floating_bodies_equal': len(old_floating),
        'baseline_relation': baseline_relation,
        'candidate_final_disassembly_sha256': hashlib.sha256(candidate_path.read_bytes()).hexdigest()})
target = work / 'actual-products-review.json'
assert not target.exists()
target.write_text(json.dumps({'passed': True, 'source_revision': '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc',
    'comparisons': reports, 'scope': 'Actual final modules with resolved typed constants/callee identities and complete guarded closure/fallback/wire joins. No numerical or runtime benefit verdict.'}, indent=2) + '\n')
print(json.dumps({'passed': True, 'candidate_final_modules': len(reports),
    'old_normal_floating_bodies': [report['old_normal_floating_bodies_equal'] for report in reports]}))
