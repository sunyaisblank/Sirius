"""Original full binary32 arithmetic; private verified FP32 restoring proposal only.

Root operates every producer/numerical owner. The separate final reader never
runs inside an owner that it tries to accept. Stock results are immutable reuse.
"""
from pathlib import Path
import collections, hashlib, json, os, re, shlex, shutil, stat, struct, subprocess, sys
import xml.etree.ElementTree as ET
sys.dont_write_bytecode = True
W = Path(__file__).resolve().parent
R = W.parents[1]
T = R / 'bin/linux-gcc/portable-divide-fp32-estimate-restoring-arithmetic'
P = json.loads((W / 'packet.json').read_text())

def verify_static():
    frozen = json.loads((W / 'source-freeze.json').read_text())
    for name, expected in frozen['pins'].items():
        assert record(name) == expected, name
        path = Path(name)
        assert {'mode': stat.S_IMODE(path.stat().st_mode),
                'symlink_target': os.readlink(path) if path.is_symlink() else None} == frozen['filesystem'][name], name


def selected_archive(label):
    row = P['archives'][label]
    archive = Path(row['path'])
    assert content(archive / 'manifest.json') == row['manifest']
    manifest = json.loads((archive / 'manifest.json').read_text())['files']
    for member, expected in row['selected_members'].items():
        path = archive / member
        assert not path.is_symlink() and content(path) == manifest[member] == expected, member
    return archive


def verify_source():
    stock = (W / 'stock/src/sirius/kernels/portable_binary32.h').read_bytes()
    header = (W / 'candidate/src/sirius/kernels/portable_binary32.h').read_bytes()
    assert hashlib.sha256(stock).hexdigest() == '49d3008fb46cd95fe4cbc2f84d41336ac6b2e5ed0d82334c831357ce9a823a60'
    assert hashlib.sha256(header).hexdigest() == '4335001872a77a107b5da94f9a4a4b21de007080e911a82b8c585ed4bd86eeca'
    assert stock == (R / 'src/sirius/kernels/portable_binary32.h').read_bytes()
    insertion = (W / 'candidate-insertion.txt').read_bytes()
    before = (W / 'originals/canonical-prologue.txt').read_bytes()
    after = (W / 'originals/candidate-prologue.txt').read_bytes()
    assert header.count(insertion) == 1 and header.startswith(after) and stock.startswith(before)
    assert header.replace(insertion,b'',1).replace(after,before,1) == stock
    # The remaining prologue delta is exact explanatory comments, not arithmetic.
    assert [line for line in before.splitlines() if not line.startswith(b'//')] == [line for line in after.splitlines() if not line.startswith(b'//')]
    return {'canonical_header': record(W / 'stock/src/sirius/kernels/portable_binary32.h'),
            'candidate_header': record(W / 'candidate/src/sirius/kernels/portable_binary32.h'),
            'optional_insertion': record(W / 'candidate-insertion.txt'),
            'original_arithmetic_outside_optional_block_unchanged': True,
            'remaining_delta_exact_comment_prologue_only': True}


def verify_reference_acceptances():
    # Relevant immutable decisions and exact source/module receipts only. The
    # old 920-input canary capture is not recursively claimed as new consumption.
    archive = selected_archive('ordinary_jit')
    old = archive / P['archives']['ordinary_jit']['work_member']
    capture = json.loads((old / 'capture-acceptance.json').read_text())
    terminal = json.loads((old / 'independent-terminal-review.json').read_text())
    admission = json.loads((old / 'root-terminal-admission.json').read_text())
    cost_path = old / 'cost-review/cost-premise-review.json'
    cost = json.loads(cost_path.read_text())
    cost_admission = json.loads((old / 'root-cost-premise-admission.json').read_text())
    assert capture['passed'] and terminal['passed'] and admission['passed'] and cost['passed'] and cost_admission['passed']
    assert admission['source_revision'] == cost_admission['source_revision'] == P['runner_preparation_revision']
    for key, path in [('independent_terminal_review', old / 'independent-terminal-review.json'),
                      ('capture_acceptance', old / 'capture-acceptance.json')]:
        assert content(path) == {k:admission[key][k] for k in ('bytes','sha256')}
    assert cost['schema'] == 'sirius-fp32-estimate-restoring-ordinary-cost-premise-v1'
    assert cost['decision'] == cost_admission['decision'] == 'COST_PREMISE_SURVIVES_UNPROVEN'
    assert content(cost_path) == {k:cost_admission['cost_review'][k] for k in ('bytes','sha256')}
    assert cost['no_arithmetic_or_benefit_admission'] and not cost['retained_arithmetic_executed_by_canary']
    assert not cost['production_adopted'] and not cost['performance_or_coldness_or_hit_rate_claim']
    compiler = selected_archive('compiler') / P['archives']['compiler']['work_member']
    graph = selected_archive('graph') / P['archives']['graph']['work_member']
    graph_result = json.loads((graph / 'graph-result.json').read_text())
    graph_acceptance = json.loads((graph / 'graph-acceptance.json').read_text())
    assert graph_result['passed'] and graph_acceptance['passed']
    assert graph_acceptance['original_failed_passed'] is False and graph_acceptance['original_failed_child_exit'] == 1
    actual = graph_result['actual_Divide']
    assert [step['shift'] for step in actual['restoring_steps']] == [4,3,2,1,0]
    assert actual['corrective_UDiv_absent'] and actual['only_original_fallback_UDiv_remains']
    assert (compiler / 'candidate/src/sirius/kernels/portable_binary32.h').read_bytes() == (W / 'candidate/src/sirius/kernels/portable_binary32.h').read_bytes()
    verify_source()
    return {'ordinary_manifest': record(archive / 'manifest.json'),
            'capture_acceptance': record(old / 'capture-acceptance.json'),
            'independent_terminal_review': record(old / 'independent-terminal-review.json'),
            'root_terminal_admission': record(old / 'root-terminal-admission.json'),
            'ordinary_cost_review': record(cost_path), 'ordinary_cost_decision': cost['decision'],
            'root_cost_admission': record(old / 'root-cost-premise-admission.json'),
            'graph_result': record(graph / 'graph-result.json'), 'graph_acceptance': record(graph / 'graph-acceptance.json'),
            'source': verify_source(),
            'scope': 'Relevant immutable compiler/graph and ordinary native-cost prerequisite only. Original failed compiler remains false/1. No retained arithmetic, benefit, hit-rate or native accuracy premise is imported.'}


def verify_host_reuse():
    # Reuse the original numeric GPU TU object, not a refusal canary or the old
    # Float64/Int64 private Vulkan TU. Fresh linking uses canonical device/base.
    archive = selected_archive('stock')
    old = Path(P['stock_work'])
    owner = accepted_folder(old / 'compile-gpu')
    seal = json.loads((old / 'compile-gpu/output-seal.json').read_text())
    for key, path in [('source_freeze',old / 'source-freeze.json'),
                      ('operator_context',old / 'operator-context.json'),
                      ('inputs_before',old / 'compile-gpu/inputs-before.json')]:
        assert content(path) == {k:seal[key][k] for k in ('bytes','sha256')}
    original = R / 'bin/linux-gcc/portable-divide64-arithmetic/gpu.o'
    assert content(P['numeric_gpu_object']) == {k:seal['outputs'][str(original)][k] for k in ('bytes','sha256')}
    assert seal['recipe'] == json.loads((old / 'packet.json').read_text())['recipes']['compile-gpu']
    before = json.loads((old / 'compile-gpu/inputs-before.json').read_text())
    for row in P['numeric_gpu_source_joins']:
        expected = before[row['original_input_key']]
        assert content(row['actual']) == {k:expected[k] for k in ('bytes','sha256')}
        if 'canonical' in row: assert Path(row['canonical']).read_bytes() == Path(row['actual']).read_bytes()
    assert (old / 'gpu.cpp').read_bytes() == (W / 'gpu.cpp').read_bytes()
    dep = Path(P['numeric_gpu_dependency']).read_text().replace('\\\n','')
    names = shlex.split(dep.split(':',1)[1])
    assert set(names) == {row['original_path'] for row in P['numeric_gpu_source_joins']}
    canonical = selected_archive('canonical_producer')
    link_folder = canonical / 'stock-comparison/files/out/guarded-product-high-runtime/link-probe'
    accepted_folder(link_folder)
    inputs = json.loads((link_folder / 'inputs-before.json').read_text())
    for row in P['canonical_producer_joins']:
        assert content(row['actual']) == {k:inputs[row['original_input_key']][k] for k in ('bytes','sha256')}
        assert content(row['actual']) == content(row['accepted_current_producer'])
    for row in P['canonical_producer_source_joins']:
        assert content(row['path']) == {k:inputs[row['key']][k] for k in ('bytes','sha256')}
    return {'numeric_gpu_object':record(P['numeric_gpu_object']),
            'original_compile_owner':record(old / 'compile-gpu/owner.json'),
            'original_output_seal':record(old / 'compile-gpu/output-seal.json'),
            'canonical_device_object':record(P['device_object']), 'canonical_base':record(P['base_archive']),
            'canonical_prior_link_owner':record(link_folder / 'owner.json'),
            'scope':'Exact original numeric GPU source/object and .d; fresh canonical Vulkan object/base relink, no private Int64/Float64 TU, no new numeric-process feature-argument trace is claimed.'}
def canonical_json(value):
    # Exact durable JSON comparison from the accepted separate pair reader.
    # Normalize integer object keys only through their actual JSON persistence.
    normalized = json.loads(json.dumps(value, allow_nan=False))
    return json.dumps(normalized, sort_keys=True, separators=(',', ':'), allow_nan=False)


def record(path):
    path = Path(path)
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': path.stat().st_size, 'sha256': digest,
            'resolved_path': str(path.resolve(strict=True))}


def content(path):
    return {k: v for k, v in record(path).items() if k != 'resolved_path'}


def write(path, value):
    path = Path(path)
    assert not path.exists(), path
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')


def corpus_checks():
    fixture = (W / 'portable_binary32_reference.bin').read_bytes()
    assert hashlib.sha256(fixture).hexdigest() == P['original_fixture_sha256']
    assert struct.unpack_from('<III', fixture) == (0x32334250, 1, 174268)
    rows = list(struct.iter_unpack('<6I', fixture[12:]))
    assert len(rows) == 174268
    counts = {str(key): value for key, value in sorted(collections.Counter(row[0] for row in rows).items())}
    assert counts == P['corpus_operation_cardinalities'] and set(counts) == {str(n) for n in range(14)}
    assert (W / 'input.bin').read_bytes() == struct.pack('<I', 174268) + b''.join(
        struct.pack('<4I', op, a, b, 0) for op, a, b, hi, lo, valid in rows)
    assert (W / 'expected.bin').read_bytes() == b''.join(
        struct.pack('<4I', hi, lo, valid, 0) for op, a, b, hi, lo, valid in rows)
    divide = [(a, b) for op, a, b, hi, lo, valid in rows if op == 3]
    finite = [(a, b) for a, b in divide if
              0 < (a & 0x7fffffff) < 0x7f800000 and 0 < (b & 0x7fffffff) < 0x7f800000]
    def significand(word):
        word &= 0x7fffffff
        fraction = word & 0x7fffff
        return (fraction | (1 << 23)) if word >> 23 else fraction << (24 - fraction.bit_length())
    coverage = {'divide_records': len(divide), 'optional_finite_nonzero_branch_records': len(finite),
                'normalization_doubling_records': sum(significand(a) < significand(b) for a, b in finite),
                'normalization_kept_records': sum(significand(a) >= significand(b) for a, b in finite),
                'non_power_of_two_denominator_significand_records': sum(significand(b) != (1 << 23) for a, b in finite)}
    assert coverage == P['division_branch_coverage']
    return {'original_records': 174268, 'operation_cardinalities': counts,
            'division_branch_coverage': coverage,
            'scope': 'Exact original ordered fixture/provider byte joins and source-derived branch partitions; not exhaustive machine proof.'}


def accepted_folder(folder):
    folder = Path(folder)
    owner = json.loads((folder / 'owner.json').read_text())
    assert owner['status'] == 'terminal' and owner['passed'] is True and owner['child_exit'] == 0
    assert owner['stop_reason'] is None and owner['cleanup_errors'] == []
    assert owner['remaining_owned_processes'] == []
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert owner['observed_births'] and all(row['absent'] for row in owner['observed_births'])
    assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
    return owner


def accepted(label):
    return accepted_folder(W / label)


def verify_output_seal(label):
    accepted(label)
    seal = json.loads((W / label / 'output-seal.json').read_text())
    assert seal['label'] == label and seal['recipe'] == P['recipes'][label]
    assert seal['source_freeze'] == record(W / 'source-freeze.json')
    assert seal['operator_context'] == record(W / 'operator-context.json')
    assert seal['inputs_before'] == record(W / label / 'inputs-before.json')
    assert seal['outputs'] == {name: record(name) for name in P['outputs'][label]}
    return record(W / label / 'output-seal.json')


def verify_archive(path, expected):
    path = Path(path)
    assert content(path / 'manifest.json') == expected
    manifest = json.loads((path / 'manifest.json').read_text())
    for name, pin in manifest['files'].items():
        item = path / name
        assert not item.is_symlink() and content(item) == pin, item
    return record(path / 'manifest.json')


def verify_cache_recovery():
    proof = json.loads(Path(P['cache_recovery']).read_text())
    base = Path(proof['base']['path']); data = base.read_bytes()
    assert proof['base']['git_commit'] == '1dd9bf1f0c7ec8e91adb89cc3563eb705d353709'
    assert hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest() == proof['base']['git_blob']
    assert content(base) == {key: proof['base'][key] for key in ('bytes', 'sha256')}
    manifest = Path(proof['base']['manifest_path'])
    assert content(manifest)['sha256'] == proof['base']['manifest_sha256']
    assert json.loads(manifest.read_text())['files'][proof['base']['manifest_member']] == content(base)
    assert proof['raw_old_cache_found'] is False and proof['canonical_cache_modified'] is False
    assert proof['original_input'] == str(R / 'bin/linux-gcc/CMakeCache.txt')
    before, after = (proof['literal_replacement'][key].encode() for key in ('before','after'))
    assert data.count(before) == proof['literal_replacement']['count'] == 1
    recovered = Path(P['historical_exact_cache'])
    assert data.replace(before, after, 1) == recovered.read_bytes()
    assert content(recovered) == proof['recovered_pin'] == {'bytes':60385,'sha256':'c7d8e854b5219e3df5d0da448d90a5d8e6f07a6cf11f216614085fc3d460de02'}
    assert proof['original_605_frozen_pin'] == {**content(recovered), 'resolved_path':str(R / 'bin/linux-gcc/CMakeCache.txt')}
    return record(recovered)


def verify_stock_reuse():
    archive = Path(P['stock_archive']); old = Path(P['stock_work']); products = Path(P['stock_build'])
    manifest = verify_archive(archive, P['stock_archive_manifest'])
    frozen = json.loads((old / 'source-freeze.json').read_text())['pins']
    assert len(frozen) == len(P['stock_static_join_rows']) == 737
    for row in P['stock_static_join_rows']:
        assert row['pin'] == frozen[row['original']]
        assert content(row['actual']) == {k: row['pin'][k] for k in ('bytes', 'sha256')}
        if row['original'] == str(R / 'bin/linux-gcc/CMakeCache.txt'):
            assert record(row['actual']) == verify_cache_recovery()
        elif row['actual'] == row['original']:
            assert record(row['actual']) == row['pin']
        else:
            assert row['pin']['resolved_path'] == row['original'] and not Path(row['actual']).is_symlink()
    result = json.loads((old / 'arithmetic-result.json').read_text())
    independent = json.loads((old / 'independent-terminal-review.json').read_text())
    assert result['passed'] and independent['accepted']
    assert result['original_records'] == 174268 and result['complete_output_words_per_provider'] == 697072
    assert result['host_explicit64_and_both_full_gpu_outputs_equal_original_expected']
    assert len(result['owners']) == 17
    for name, pin in result['owners'].items():
        accepted_folder(old / name)
        assert content(old / name / 'owner.json') == {k: pin[k] for k in ('bytes', 'sha256')}
    owner = accepted_folder(old / 'execute-gpu-stock')
    provider_mappings(owner)
    before = json.loads((old / 'execute-gpu-stock/inputs-before.json').read_text())
    for original, actual in ((R / 'bin/linux-gcc/portable-divide64-arithmetic/gpu', products / 'gpu'),
                             (R / 'bin/linux-gcc/portable-divide64-arithmetic/stock.spv', products / 'stock.spv'),
                             (R / 'out/portable-divide64-arithmetic/input.bin', old / 'input.bin'),
                             (R / 'out/portable-divide64-arithmetic/expected.bin', old / 'expected.bin')):
        key = str(original.relative_to(R))
        assert {k: before[key][k] for k in ('bytes', 'sha256')} == content(actual)
    for name in ('portable_binary32_reference.bin', 'input.bin', 'expected.bin', 'probe.slang', 'gpu.cpp'):
        assert (W / name).read_bytes() == (old / name).read_bytes()
    assert (old / 'execute-gpu-stock/actual.bin').read_bytes() == (W / 'expected.bin').read_bytes()
    rows = [json.loads(line) for line in (old / 'execute-gpu-stock/stdout.log').read_text().splitlines() if line.startswith('{')]
    assert len(rows) == 1 and rows[0]['records'] == 174268 and rows[0]['complete_words'] == 697072 and rows[0]['mismatches'] == 0
    assert rows[0]['driver'] == 'llvmpipe'
    return {'archive_manifest': manifest, 'accepted_arithmetic_result': record(old / 'arithmetic-result.json'),
            'independent_terminal_review': record(old / 'independent-terminal-review.json'),
            'runtime_context': record(old / 'runtime-context.json'), 'build_review': record(old / 'build-review.json'),
            'owner': record(old / 'execute-gpu-stock/owner.json'), 'consumed_inputs': record(old / 'execute-gpu-stock/inputs-before.json'),
            'actual': record(old / 'execute-gpu-stock/actual.bin'), 'gpu': record(products / 'gpu'),
            'stock_module': record(products / 'stock.spv'), 'original_input': record(old / 'input.bin'),
            'original_expected': record(old / 'expected.bin'), 'all737_static_inputs_and17_terminal_owners_joined': True,
            'scope': 'Exact accepted complete stock arithmetic/source/build/producer/corpus/provider receipt reused. No fresh stock compile or dispatch is claimed.'}


def dependencies(name):
    text = (T / name).read_text().replace('\\\n', '')
    return [Path(path).resolve(strict=True) for path in shlex.split(text.split(':', 1)[1])]


def provider_mappings(owner):
    expected = set(P['provider_identities'])
    observed = set()
    births = {(row['pid'], row['saved_start_ticks']) for row in owner['observed_births']}
    for row in owner['provider_mappings']:
        assert (row['pid'], row['start_ticks']) in births
        assert row['owned_session'] == owner['child_pid']
        for line in row['maps']:
            fields = line.split(None, 5)
            assert len(fields) == 6 and fields[5] in expected
            identity = P['provider_identities'][fields[5]]
            device = [int(part, 16) for part in fields[3].split(':')]
            assert device == identity['device'] and int(fields[4]) == identity['inode']
            observed.add(fields[5])
    assert observed == expected


def join_runtime_inputs(label, owner):
    action = json.loads((W / 'actions.json').read_text())[label]
    folder = W / label
    assert owner['argv'] == [arg.replace('{folder}', str(folder)) for arg in action['argv']]
    before = json.loads((folder / 'inputs-before.json').read_text())
    for name in action['frozen_inputs']:
        path = Path(name)
        key = str(path.relative_to(R)) if path.is_relative_to(R) else str(path)
        assert before[key] == record(path), ('numerical consumed input differs from current sealed input', label, name)


def primitive_policy(text, annotated):
    from interface_signature import function_body
    caps = re.findall(r'OpCapability (\S+)',text)
    assert caps == ['Shader']
    assert set(re.findall(r'OpTypeInt (\d+) [01]',text)) == {'32'}
    floating = re.findall(r'(%\S+) = OpTypeFloat (\d+)',text)
    assert len(floating) == 1 and floating[0][1] == '32'
    divide = function_body(text,'PB32Divide')
    proposals = re.findall(r'^\s*(%\S+) = OpFDiv '+re.escape(floating[0][0])+r' .*$',text,re.M)
    assert len(proposals) == 1 and re.findall(r'= (OpF(?!unction)\S+) ',text) == ['OpFDiv']
    proposal = proposals[0]
    assert re.search(r'^\s*'+re.escape(proposal)+r' = OpFDiv ',divide,re.M)
    assert not re.search(r'\b(?:RelaxedPrecision|FPFastMathMode|FPFastMathDefault|OpFmaKHR|OpConvertFToU|OpConvertFToS|OpConvertUToF|OpConvertSToF|NoUnsignedWrap|NoSignedWrap)\b',text)
    assert re.findall(r'OpExecutionMode %\S+ (.*)',text) == ['LocalSize 64 1 1']
    assert text.count('OpControlBarrier') == text.count('OpMemoryBarrier') == 0
    decorations = re.findall(r'OpDecorate (%\S+) NoContraction',text)
    assert decorations == ([proposal] if annotated else [])
    return proposal


def annotate_primitive():
    raw = (T / 'candidate.raw.spvasm').read_text()
    proposal = primitive_policy(raw,False)
    annotation = '\n               OpDecorate '+proposal+' NoContraction'
    declarations = list(re.finditer(r'^.*OpDecorate .*?$',raw,re.M))
    assert declarations
    point = declarations[-1].end()
    text = raw[:point]+annotation+raw[point:]
    assert text.count(annotation) == 1 and text.replace(annotation,'',1) == raw
    assert primitive_policy(text,True) == proposal
    (T / 'candidate.spvasm.input').write_text(text)


def semantic_instructions(text):
    ids = re.findall(r'^\s*(%\S+) =',text,re.M)
    assert len(ids) == len(set(ids))
    aliases = {name:'id'+str(i) for i,name in enumerate(ids)}
    result = []
    for line in text.splitlines():
        lexer = shlex.shlex(line,posix=False); lexer.whitespace_split = True; lexer.commenters = ';'
        tokens = list(lexer)
        if tokens: result.append([aliases[token] if token.startswith('%') else token for token in tokens])
    return result


def module_checks():
    from interface_signature import interface_signature, function_body
    import graph, helper_normalize
    texts = {name:(T / (name+'.spvasm')).read_text() for name in ('stock','candidate')}
    assert (T / 'stock.spv').read_bytes() == Path(P['original_stock_module']).read_bytes()
    raw = (T / 'candidate.raw.spvasm').read_text()
    assembly = (T / 'candidate.spvasm.input').read_text()
    proposal = primitive_policy(raw,False)
    annotation = '\n               OpDecorate '+proposal+' NoContraction'
    point = list(re.finditer(r'^.*OpDecorate .*?$',raw,re.M))[-1].end()
    assert assembly == raw[:point]+annotation+raw[point:]
    assert assembly.replace(annotation,'',1) == raw
    assert semantic_instructions(assembly) == semantic_instructions(texts['candidate'])
    primitive_policy(texts['candidate'],True)
    assert re.findall(r'OpCapability (\S+)',texts['stock']) == ['Shader']
    assert not re.findall(r'OpTypeFloat (\d+)',texts['stock'])
    assert interface_signature(texts['stock']) == interface_signature(texts['candidate'])
    helpers = helper_normalize.compare_helpers(texts['stock'],texts['candidate'])
    actual = graph.inspect(texts['stock'],texts['candidate'],function_body)
    assert [step['shift'] for step in actual['restoring_steps']] == [4,3,2,1,0]
    assert actual['corrective_UDiv_absent'] and actual['only_original_fallback_UDiv_remains']
    results = {}
    for name,text in texts.items():
        widths = dict(re.findall(r'(%\S+) = OpTypeInt (\d+) [01]',text))
        assert set(widths.values()) == {'32'}
        assert re.findall(r'OpExecutionMode %\S+ (.*)',text) == ['LocalSize 64 1 1']
        assert text.count('OpControlBarrier') == text.count('OpMemoryBarrier') == 0
        entry = function_body(text,'ComputeMain')
        defs = dict(re.findall(r'^\s*(%\S+) = (.+)$',text,re.M))
        names = dict(re.findall(r'OpName (%\S+) "([^"]*)"',text))
        rows,location,current = {},{},None
        for line in entry.splitlines():
            line=line.strip();found=re.match(r'(%\S+) = OpLabel$',line)
            if found:current=found[1];rows[current]=[]
            if current is not None:
                rows[current].append(line)
                if re.match(r'%\S+ = ',line):location[line.split()[0]]=current
        calls=[key for key in location if defs[key].split()[0]=='OpFunctionCall' and names.get(defs[key].split()[2])=='PB32Divide']
        assert len(calls)==1
        call=calls[0];block=location[call];values=defs[call].split();assert len(values)==5
        branches=[(label,line.split()) for label,lines in rows.items() for line in lines if line.startswith('OpBranchConditional ') and line.split()[2]==block]
        assert len(branches)==1
        header,branch=branches[0];condition=defs[branch[1]].split()
        assert condition[0]=='OpIEqual' and len(condition)==4
        constants=[value for value in condition[2:] if defs[value].split()[0]=='OpConstant' and defs[value].split()[2]=='3']
        assert len(constants)==1
        operation=next(value for value in condition[2:] if value!=constants[0])
        for value in (operation,*values[3:]):
            loaded=defs[value].split();assert loaded[:2]==['OpLoad',values[1]]
            pointer=defs[loaded[2]].split();assert pointer[0]=='OpAccessChain' and names.get(pointer[2])=='inputs'
        merge=[line.split()[1] for line in rows[header] if line.startswith('OpSelectionMerge ')]
        assert len(merge)==1 and len(rows[block])==4
        assert rows[block][2].split()[0]=='OpStore' and rows[block][2].split()[2]==call
        assert rows[block][-1]=='OpBranch '+merge[0]
        predecessors=set()
        for label,lines in rows.items():
            for line in lines:
                words=line.split()
                targets=(words[1:] if words[:1]==['OpBranch'] else words[2:4] if words[:1]==['OpBranchConditional'] else [words[2],*words[4::2]] if words[:1]==['OpSwitch'] else [])
                if block in targets:predecessors.add(label)
        assert predecessors=={header}
        divide=function_body(text,'PB32Divide')
        assert [widths[k] for k in re.findall(r'= OpUDiv (%\S+) ',divide)]==['32'] and divide.count('OpLoopMerge')==1
        assert all(widths[k]=='32' for k in re.findall(r'= Op[US](?:Div|Mod|Rem) (%\S+) ',text))
        results[name]={'module':record(T/(name+'.spv')),'assembly':record(T/(name+'.spvasm')),
            'capabilities':re.findall(r'OpCapability (\S+)',text),'integer_widths':sorted(set(widths.values())),
            'float_widths':re.findall(r'OpTypeFloat (\d+)',text),
            'op3_actual_Divide_call':{'id':call,'condition':branch[1],'case3_true_block_only':True,'result_stored':True},
            'divide_estimate_graph':actual if name=='candidate' else None,
            'original_complete_prefix_and_fallback_suffix_normalized_equal':True,
            'all_non_Divide_helper_instances_normalized_equal':helpers,
            'whole_assembler_input_readback_instructions_equal':len(semantic_instructions(assembly)) if name=='candidate' else None,
            'original_integer_public_interface_and_LocalSize64_retained':True,
            'no_Endpoint_RTE_lines_added_to_primitive':True}
    return results


def build_review():
    seals={label:verify_output_seal(label) for label in P['build_labels']}
    pins=json.loads((W/'source-freeze.json').read_text())['pins']
    for path in dependencies('host.d'):
        assert str(path) in pins and record(path)==pins[str(path)],path
    assert W/'candidate/src/sirius/kernels/portable_binary32.h' in dependencies('host.d')
    assert W/'stock/src/sirius/kernels/portable_binary32.h' not in dependencies('host.d')
    assert '-DSIRIUS_PB32_DIVIDE_FP32_ESTIMATE=1' in P['recipes']['compile-host'][0]
    assert '-DSIRIUS_PB32_DIVIDE_FP32_ESTIMATE=1' in P['recipes']['compile-shader-candidate'][0]
    needed=json.loads(Path(P['original_symbols']).read_text())['required_strong_sirius_definitions']
    strong=lambda text:{row.split()[2] for row in text.splitlines() if len(row.split())==3 and row.split()[1] in ('T','D','B','R')}
    symbols={}
    for name,binary in [('device-object',T/'vulkan_device.cpp.o'),('gpu',T/'gpu')]:
        text=subprocess.check_output(['/usr/bin/nm','-g','--defined-only',str(binary)],text=True)
        assert set(needed)<=strong(text);symbols[name]=text
    link_map=(T/'gpu.map').read_text()
    assert 'libsirius_backend.a' not in link_map and 'libsirius_core.a' not in link_map
    argv=P['recipes']['link-gpu'][0]
    assert argv==P['canonical_gpu_link_recipe']
    assert str(T/'vulkan_device.cpp.o') in argv and str(T/'numeric-gpu.o') in argv
    assert (T/'resources/tests/backend/portable_binary32_reference.bin').read_bytes()==(W/'portable_binary32_reference.bin').read_bytes()
    return {'passed':True,'modules':module_checks(),'output_seals':seals,'original_corpus':corpus_checks(),
            'canonical_numeric_object_relink':verify_host_reuse(),'stock_arithmetic_receipt':verify_stock_reuse(),
            'relevant_ordinary_prerequisites':verify_reference_acceptances(),'source':verify_source(),
            'required_strong_definitions':needed,'nm_outputs':symbols,
            'private_Int64_or_Float64_device_object_used':False,
            'original_backend_or_core_archive_extracted':False,
            'host_uses_candidate_header_and_explicit_macro':True,'probe_and_all_original_assertions_unchanged':True}


def pin_context():
    verify_static();verify_source();verify_host_reuse();verify_stock_reuse();verify_reference_acceptances();corpus_checks()
    assert not T.exists() and not (W/'operator-context.json').exists()
    assert not subprocess.check_output(['git','status','--porcelain=v1'],cwd=R)
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()
    assert revision==P['runner_preparation_revision']
    actions=json.loads((W/'actions.json').read_text())
    assert all(action['revision']==revision for action in actions.values())
    write(W/'operator-context.json',{'revision':revision,'actions':record(W/'actions.json'),
          'source_freeze':record(W/'source-freeze.json'),'preparation':record(W/'preparation.json')})
    print('Pinned original FP32-restoring arithmetic preparation; no compiler/provider launched.',flush=True)


def phase(label):
    verify_static()
    context=json.loads((W/'operator-context.json').read_text())
    for key,name in [('actions','actions.json'),('source_freeze','source-freeze.json'),('preparation','preparation.json')]:assert context[key]==record(W/name)
    assert label in P['outputs']
    for prior in P['dependencies'][label]:verify_output_seal(prior)
    assert all(not Path(path).exists() for path in P['outputs'][label])
    if label=='prepare-build':
        assert not T.exists();(T/'resources/tests/backend').mkdir(parents=True)
        (T/'resources/tests/backend/portable_binary32_reference.bin').write_bytes((W/'portable_binary32_reference.bin').read_bytes())
        for destination,original in P['reused_files'].items():shutil.copyfile(original,destination)
    elif label=='verify-build':write(W/'build-review.json',build_review())
    else:
        if label=='assemble-shader-candidate':annotate_primitive()
        for argv in P['recipes'][label]:subprocess.run(argv,cwd=R,check=True)
    write(W/label/'output-seal.json',{'label':label,'recipe':P['recipes'][label],
          'source_freeze':record(W/'source-freeze.json'),'operator_context':record(W/'operator-context.json'),
          'inputs_before':record(W/label/'inputs-before.json'),'outputs':{name:record(name) for name in P['outputs'][label]}})


def pin_runtime():
    verify_static();assert not (W/'runtime-context.json').exists()
    verify_output_seal('verify-build');build=json.loads((W/'build-review.json').read_text())
    assert build['passed'] and canonical_json(build['modules'])==canonical_json(module_checks())
    write(W/'runtime-context.json',{'build_review':record(W/'build-review.json'),
          'all_build_output_seals':{label:verify_output_seal(label) for label in P['outputs']},
          'gpu':record(T/'gpu'),'host':record(T/'host'),'modules':build['modules'],
          'original_input':record(W/'input.bin'),'original_expected':record(W/'expected.bin'),
          'canonical_numeric_object_relink':verify_host_reuse(),
          'scope':'Original complete binary32 primitive corpus only. Untrusted native FP32 proposal; exact integer guard/correction/fallback, no hit-rate, retained-program, benefit or qualification claim.'})
    print('Original arithmetic runtime inputs pinned; no device launched.',flush=True)


def review():
    verify_static();context=json.loads((W/'operator-context.json').read_text())
    assert context['actions']==record(W/'actions.json') and context['preparation']==record(W/'preparation.json')
    runtime=json.loads((W/'runtime-context.json').read_text())
    assert runtime['build_review']==record(W/'build-review.json')
    assert runtime['all_build_output_seals']=={label:verify_output_seal(label) for label in P['outputs']}
    assert runtime['gpu']==record(T/'gpu') and runtime['host']==record(T/'host')
    assert canonical_json(runtime['modules'])==canonical_json(module_checks())
    assert runtime['original_input']==record(W/'input.bin') and runtime['original_expected']==record(W/'expected.bin')
    assert runtime['canonical_numeric_object_relink']==verify_host_reuse()
    for label in ('execute-host','execute-gpu-candidate'):
        owner=accepted(label);join_runtime_inputs(label,owner)
        if label=='execute-gpu-candidate':
            provider_mappings(owner)
            assert (W/label/'actual.bin').read_bytes()==(W/'expected.bin').read_bytes()
            rows=[json.loads(line) for line in (W/label/'stdout.log').read_text().splitlines() if line.startswith('{')]
            assert len(rows)==1 and rows[0]['records']==174268 and rows[0]['complete_words']==697072 and rows[0]['mismatches']==0 and rows[0]['driver']=='llvmpipe'
    xml=ET.parse(W/'execute-host/gtest.xml').getroot()
    assert int(xml.get('tests','-1'))==1 and all(int(xml.get(k,'0'))==0 for k in ('failures','errors','disabled'))
    assert not any(xml.findall('.//'+tag) for tag in ('failure','error','skipped'))
    cases=xml.findall('.//testcase');assert len(cases)==1
    case=cases[0];assert case.get('classname')+'.'+case.get('name')=='PortableBinary32.MatchesIndependentExactCorpus'
    assert case.get('status')=='run' and case.get('result')=='completed'
    properties={node.get('name'):node.get('value') for node in case.findall('./properties/property')}
    assert properties['independent_exact_records_checked']=='174268' and properties['exact_raw_word_mismatches']=='0'
    write(W/'arithmetic-result.json',{'schema':'sirius-fp32-estimate-restoring-original-arithmetic-result-v1',
          'passed':True,'original_records':174268,'original_operation_count':14,'complete_output_words_per_provider':697072,
          'candidate_macro':P['candidate_macro'],'candidate_header':record(W/'candidate/src/sirius/kernels/portable_binary32.h'),
          'candidate_shader':record(T/'candidate.spv'),'modules':module_checks(),
          'stock_result_reused':True,'stock_arithmetic_receipt':verify_stock_reuse(),
          'relevant_ordinary_prerequisites':verify_reference_acceptances(),'build_review':record(W/'build-review.json'),
          'fresh_candidate_host_and_full_gpu_plus_reused_stock_equal_original_expected':True,
          'runtime_context':record(W/'runtime-context.json'),'original_corpus':corpus_checks(),
          'owners':{label:record(W/label/'owner.json') for label in [*P['outputs'],'execute-host','execute-gpu-candidate']},
          'scope':'Unchanged original174268/14-operation primitive corpus; fresh candidate host/full GPU plus exact stock receipt. Actual proposal acceptance/hit-rate unobserved; no retained-program, eight-control, stage, performance, adoption or qualification verdict.'})
    print('Original174268 candidate host/provider outputs accepted; exact complete stock receipt reused.',flush=True)


if __name__=='__main__':
    if sys.argv[1]=='phase':phase(sys.argv[2])
    else:{'pin-context':pin_context,'pin-runtime':pin_runtime,'review':review}[sys.argv[1]]()

