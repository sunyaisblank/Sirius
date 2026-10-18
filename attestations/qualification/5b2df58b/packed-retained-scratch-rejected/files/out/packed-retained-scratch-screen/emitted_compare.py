"""Bounded saved-SPIR-V-text comparison only; no project imports or tools."""
from pathlib import Path
from collections import Counter
import functools
import hashlib
import json
import re

Q = Path(__file__).resolve().parent
ROOT = Q.parents[1]
BASE = ROOT / 'bin/linux-gcc/src/sirius/backend/retained'
ASM = Q / 'posttools/disassemble'
HELPERS = {'ScratchAddress', 'ScratchRead', 'ScratchWrite'}
TOKEN = re.compile(r'"(?:\\.|[^"\\])*"|[^\s]+')
IDENT = re.compile(r'%[^\s]+')


def pin(path):
    return {'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':')).encode()).hexdigest()


class Module:
    def __init__(self, path):
        self.path = path
        self.text = path.read_text()
        self.lines = [line.strip() for line in self.text.splitlines()
                      if line.strip() and not line.lstrip().startswith(';')]
        self.defs = dict(re.findall(r'^\s*(%\S+) = (.+)$', self.text, re.M))
        self.names = dict(re.findall(r'OpName (%\S+) "([^"]*)"', self.text))
        self.bodies = dict(re.findall(r'^\s*(%\S+) = (OpFunction\b.*?OpFunctionEnd)',
                                      self.text, re.M | re.S))
        counts = Counter()
        self.identities = {}
        self.by_identity = {}
        self.locals = {}
        for ident, body in self.bodies.items():
            name = self.names[ident]
            key = (name, counts[name])
            counts[name] += 1
            self.identities[ident] = key
            self.by_identity[key] = ident
            values = [ident] + re.findall(r'^\s*(%\S+) =', body, re.M)
            self.locals[ident] = {value: n for n, value in enumerate(values)}
        self.local_owner = {value: (self.identities[ident], n)
                            for ident, values in self.locals.items()
                            for value, n in values.items()}
        self.globals = {ident: expr for ident, expr in self.defs.items()
                        if expr.startswith('OpVariable ') and
                        self.tokens(expr)[2] != 'Function'}
        self.global_ids = {ident: self.names.get(ident, ident) for ident in self.globals}
        assert len(set(self.global_ids.values())) == len(self.global_ids)
        self.scratch = next(ident for ident in self.globals
                            if self.global_ids[ident] == 'retainedScratch')
        pointer = self.tokens(self.globals[self.scratch])[1]
        pt = self.tokens(self.defs[pointer])
        assert pt[:2] == ['OpTypePointer', 'Workgroup']
        self.array = pt[2]
        ar = self.tokens(self.defs[self.array])
        assert ar[0] == 'OpTypeArray' and self.tokens(self.defs[ar[1]]) == ['OpTypeInt', '32', '0']
        self.words = self.integer(ar[2])
        self.norm = functools.lru_cache(None)(self._norm)

    @staticmethod
    def tokens(line):
        return TOKEN.findall(line)

    def integer(self, ident):
        parts = self.tokens(self.defs[ident])
        assert parts[0] == 'OpConstant' and len(parts) == 3
        return int(parts[2], 0)

    def _norm(self, ident):
        if ident in self.identities:
            return ('function', self.identities[ident])
        if ident in self.globals:
            return ('global', self.global_ids[ident])
        if ident in self.local_owner:
            return ('function-local', self.local_owner[ident])
        expr = self.defs[ident]
        assert expr.startswith(('OpType', 'OpConstant', 'OpSpecConstant',
                                'OpExtInstImport ', 'OpUndef ')), (ident, expr)
        return tuple(self.norm(x) if x.startswith('%') else x for x in self.tokens(expr))

    def body(self, key):
        ident = self.by_identity[key]
        locals_ = self.locals[ident]
        rows = [ident + ' = ' + self.bodies[ident].splitlines()[0]]
        rows += self.bodies[ident].splitlines()[1:]
        return tuple(tuple(('local', locals_[x]) if x in locals_ else
                           self.norm(x) if x.startswith('%') else x
                           for x in self.tokens(line)) for line in rows if line.strip())

    def global_declarations(self):
        rows = {}
        for ident, expr in self.globals.items():
            if ident == self.scratch:
                rows[self.global_ids[ident]] = ('OpVariable', ('uint32-array', 'separate-length-check'), 'Workgroup')
            else:
                rows[self.global_ids[ident]] = tuple(self.norm(x) if x.startswith('%') else x
                                                   for x in self.tokens(expr))
        return rows

    def semantic_header(self):
        rows = []
        for line in self.lines:
            if not line.startswith(('OpCapability ', 'OpExtension ', 'OpMemoryModel ',
                                    'OpEntryPoint ', 'OpExecutionMode ', 'OpExecutionModeId ',
                                    'OpDecorate ', 'OpMemberDecorate ', 'OpDecorateId ')):
                continue
            rows.append(tuple(('scratch-array', 'uint32') if x == self.array else
                              self.norm(x) if x.startswith('%') else x
                              for x in self.tokens(line)))
        return rows

    def counters(self, key):
        body = self.bodies[self.by_identity[key]]
        return dict(sorted(Counter(re.findall(r'\b(Op[A-Za-z0-9]+)\b', body)).items()))

    def locus(self, key):
        ident = self.by_identity[key]
        found = re.search(r'^\s*' + re.escape(ident) + r' = OpFunction\b', self.text, re.M)
        return self.text.count('\n', 0, found.start()) + 1


def read_json(path):
    return json.loads(path.read_text())


def main():
    target = Q / 'independent-emitted-review.json'
    assert not target.exists()
    report = read_json(Q / 'root-compiler-products.json')
    producer = read_json(Q / 'candidate/owner.json')
    post = read_json(ASM / 'owner.json')
    for owner in (producer, post):
        assert owner['passed'] is True and owner['child_exit'] == 0
        assert owner['whole_inputs_unchanged'] is True and owner['source_unchanged'] is True
        assert owner['stop_reason'] is None and not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    assert pin(Q / 'candidate/owner.json') == report['owner'] and report['passed'] is True
    commands = read_json(ASM / 'results.json')
    assert len(commands) == 50 and all(type(x['returncode']) is int and x['returncode'] == 0 for x in commands)
    by_label = {x['label']: x for x in commands}
    assert len(by_label) == 50
    before = read_json(ASM / 'inputs-before.json')
    assert before == read_json(ASM / 'inputs-after.json')
    for row in commands:
        tool, product, option, output = row['argv']
        assert option == '-o'
        physical = Path(product)
        expected = ASM / (row['label'] + '.spvasm')
        assert Path(output) == expected and physical.is_file() and expected.is_file()
        for path in [Path(tool), physical]:
            key = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
            assert {k: before[key][k] for k in ['bytes', 'sha256']} == pin(path)
    source_review = ROOT / 'out/packed-retained-scratch-source/independent-source-review.json'
    assert read_json(source_review)['passed'] is True
    products = report['products']
    assert len(products) == 35 and sum(x['excluded'] for x in products) == 10
    cases = []
    patterns = {}
    for product in products:
        name = product['name']
        a_path, b_path = BASE / name, Q / 'candidate' / name
        assert pin(a_path) == product['baseline'] and pin(b_path) == product['candidate'], name
        if product['excluded']:
            assert a_path.read_bytes() == b_path.read_bytes(), name
            continue
        stem = Path(name).stem
        left, right = Module(ASM / ('baseline-' + stem + '.spvasm')), Module(ASM / ('candidate-' + stem + '.spvasm'))
        assert by_label['baseline-'+stem]['argv'][1] == str(a_path)
        assert by_label['candidate-'+stem]['argv'][1] == str(b_path)
        assert set(left.by_identity) == set(right.by_identity), name
        assert left.words % 5 == 0 and right.words == 4 * (left.words // 5), name
        assert left.global_declarations() == right.global_declarations(), ('globals', name)
        assert left.semantic_header() == right.semantic_header(), ('header/decorations', name)
        unchanged = []
        changed = []
        helper_keys = {}
        for key in left.by_identity:
            a, b = left.body(key), right.body(key)
            if key[0] in HELPERS:
                changed.append(key)
                assert a != b
                item = {'baseline_alpha_sha256': fingerprint(a), 'candidate_alpha_sha256': fingerprint(b),
                        'baseline_ops': left.counters(key), 'candidate_ops': right.counters(key)}
                pattern_key = fingerprint(item)
                variants = patterns.setdefault(key[0], {})
                assert pattern_key not in variants or variants[pattern_key] == item
                variants[pattern_key] = item
                helper_keys[key[0]] = pattern_key
            else:
                assert a == b, ('outside-helper-body', name, key)
                unchanged.append((key, fingerprint(a)))
        assert set(changed) == {(name, 0) for name in HELPERS}
        # The complete WriteValue body must be unchanged; additionally bind ordered component calls.
        write = right.bodies[right.by_identity[('WriteValue', 0)]]
        calls = [right.tokens(line.split(' = ', 1)[-1]) for line in write.splitlines()
                 if 'OpFunctionCall ' in line]
        writes = [x for x in calls if right.identities[x[2]] == ('ScratchWrite', 0)]
        bits = [x for x in calls if right.identities[x[2]] == ('RPBits', 0)]
        assert len(writes) == 5 and len(bits) in (0, 4) and len(calls) == len(writes) + len(bits)
        assert [right.integer(x[4]) for x in writes] == [0, 1, 2, 3, 4]
        assert len({(x[3], x[5], x[6]) for x in writes}) == 1
        cases.append({'name': name, 'baseline_asm': pin(left.path), 'candidate_asm': pin(right.path),
                      'scratch_uint_words': [left.words, right.words], 'unchanged_function_bodies': len(unchanged),
                      'unchanged_body_set_sha256': fingerprint(unchanged),
                      'candidate_helper_lines': {key[0]: right.locus(key) for key in sorted(changed)},
                      'helper_pattern_keys': helper_keys,
                      'write_value_line': right.locus(('WriteValue', 0)),
                      'interface_decorations_sha256': fingerprint(left.semantic_header())})
    assert len(cases) == 25 and set(patterns) == HELPERS
    # Word and declaration equality, not source-builder execution.
    a, b = (BASE / 'retained_kernels.h').read_text(), (Q / 'candidate/retained_kernels.h').read_text()
    programs = []
    for item in report['programs']:
        pattern = r'\bk' + re.escape(item['name'][1:]) + r'\{\{(.*?)\}\};'
        aa = re.search(pattern, a, re.S).group(1)
        bb = re.search(pattern, b, re.S).group(1)
        assert aa == bb
        words = re.findall(r'\b(\d+)u\b', aa)
        assert len(words) == item['words']
        programs.append({'name': item['name'], 'words': len(words), 'exact_body_sha256': hashlib.sha256(aa.encode()).hexdigest()})
    assert len(programs) == 7
    for name in report['owner_products_byte_exact']:
        assert (BASE / name).read_bytes() == (Q / 'candidate' / name).read_bytes()
    result = {'schema': 'sirius-packed-retained-scratch-independent-emitted-review-v1',
              'passed': True, 'source_revision': '5b2df58bb9f24e7db5f7fbc0c38bdd832a58900e',
              'scope': 'All25 changed SPV saved-text comparisons and35product/7table/exclusion joins only; source author also authored packing prototype; ROOT alone produced/disassembled. Not independent raw producer or a whole-program formal proof.',
              'input_pins': {'root_compiler_products': pin(Q / 'root-compiler-products.json'),
                             'compiler_owner': pin(Q / 'candidate/owner.json'),
                             'posttools_owner': pin(ASM / 'owner.json'), 'posttools_results': pin(ASM / 'results.json'),
                             'posttools_before': pin(ASM / 'inputs-before.json'), 'posttools_after': pin(ASM / 'inputs-after.json'),
                             'source_review': pin(source_review), 'reader': pin(Path(__file__))},
              'comparison': 'Full ordered function instructions normalized by local-definition ordinal, exact recursive type/constant/import definitions, function(name,occurrence), and separately exact global declarations. All function bodies except exactly3Scratch helpers equal; includes fullCFG/barriers/calls/WriteValue component0..4 order. Full semantic header/decorations equal except separately checked scratch-array length.',
              'products': {'total': 35, 'exact_exclusions': 10, 'changed_pairs': 25, 'programs': programs,
                           'owner_products_byte_exact': report['owner_products_byte_exact']},
              'helper_lowering_patterns': patterns, 'cases': cases,
              'emitted_semantics': {'address': 'uint32 4*index+(component==4?3:component), same signature; native generic branch and Function-local store/load, portable SSA Phi selecting3/component retained.',
                                    'read': 'Single packed load; component3 low31 mask, component4 logical>>31, components0..2 direct. Two conditional branches retained.',
                                    'write': 'Component3 masks value and clears high validity; component4 loads oldword, masks low31, ORs (value&1)<<31 and stores; others direct. Two conditional branches and read-modify-write retained.',
                                    'publication': 'All25 complete WriteValue bodies alpha-exact and five same-index component0..4 calls ordered. Outside-helper CFG/barrier equality composes with admitted source single-writer/layer ownership proof. No intermediate-value or opaque initialization equality claimed.'},
              'decision': 'EMITTED_CONTRACT_PASS_COST_UNPROVEN', 'runtime_admitted': False, 'production_adopted': False,
              'limits': ['Source4R versus5R Workgroup words verified, separate status/public types unchanged. No driver LDSallocation/occupancy/load-count/speed claim.',
                         'Generic helper branches, Functionlocal address temporary and extra shared RMW load remain at originalO0; no Slang store-forwarding observed. All25SPVgrow(+732B native12,+704B portable13); portableSSA replaces address Functionlocal withPhi, while packed RMW stays. Actual driver lowering/cost unknown.',
                         'Validity0/1 and error-sign eligibility, S7 exclusion and opaque initialization limitations reuse admitted independentSOURCEreview. This saved-text comparison is not numerical/device/failure-injection qualification.',
                         'Original scientific/reference/fullword/refusal/recovery controls and benefit gates still required before adoption; no frame retry or runtime released.']}
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'passed': True, 'changed_pairs': len(cases), 'receipt': pin(target)}))


if __name__ == '__main__':
    main()
