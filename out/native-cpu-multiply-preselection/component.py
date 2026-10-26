"""One exact native-multiply lifetime change; reuse the original ABI fixture."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
BUILD = ROOT / 'bin/linux-gcc/native-cpu-multiply-preselection'
BASIS = ROOT / 'out/native-cpu-value-fact-stores'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,
                    sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def write(path, value):
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2) + '\n')


def verified_inputs():
    prepared = doc(BASIS / 'prepared-cpp.json')
    for name, expected in prepared['inputs'].items():
        assert pin(ROOT / name) == expected, name
    built = doc(BASIS / 'cpp-build-result.json')
    assert built['passed']
    for row in built['compilations']:
        assert pin(ROOT / row['object']['path']) == {
            key: row['object'][key] for key in ('bytes', 'sha256')}
    return prepared, built


phase, = sys.argv[1:]
if phase == 'derive':
    prepared, built = verified_inputs()
    original = (BASIS / 'retained8.cpp').read_text()
    before = '''    U work=active&good&~is_zero&~a_one&~b_one&~a_minus&~b_minus;Triple general=invalid();
    if(any(work)) {
        U domain=nz(ta&u(1))&nz(tb&u(1))&bounded(a.hi,1518338048U)&bounded(b.hi,1518338048U);
        if(!all_on(domain,work)) return multiply_triple(a,b,active);
        Triple aa=choose(work,a,constant_zero()),bb=choose(work,b,constant_zero());
        general=facts_multiply_general(aa,bb,work);
    }
    Triple selected=choose(a_minus,negate(b),choose(b_minus,negate(a),general));
    selected=choose(a_one,b,choose(b_one,a,selected));
    selected=choose(is_zero,constant_zero(),selected);
    return choose(good,selected,invalid());'''
    after = '''    // Keep the all-lane partition: inactive general lanes still receive the
    // original general result when another lane executes the arithmetic.
    U general_lanes=good&~is_zero&~a_one&~b_one&~a_minus&~b_minus;
    U work=active&general_lanes;
    Triple selected=choose(a_minus,negate(b),choose(b_minus,negate(a),invalid()));
    selected=choose(a_one,b,choose(b_one,a,selected));
    selected=choose(is_zero,constant_zero(),selected);
    selected=choose(good,selected,invalid());
    Triple general=invalid();
    if(any(work)) {
        U domain=nz(ta&u(1))&nz(tb&u(1))&bounded(a.hi,1518338048U)&bounded(b.hi,1518338048U);
        if(!all_on(domain,work)) return multiply_triple(a,b,active);
        Triple aa=choose(work,a,constant_zero()),bb=choose(work,b,constant_zero());
        general=facts_multiply_general(aa,bb,work);
    }
    return choose(general_lanes,general,selected);'''
    assert original.count(before) == 1
    candidate = original.replace(before, after, 1)
    assert candidate.replace(after, before, 1) == original
    destination = WORK / 'retained8.cpp'
    with destination.open('x') as stream:
        stream.write(candidate)
    write(WORK / 'derivation.json', dict(
        basis=dict(path=str((BASIS / 'retained8.cpp').relative_to(ROOT)), **pin(BASIS / 'retained8.cpp')),
        candidate=dict(path=str(destination.relative_to(ROOT)), **pin(destination)),
        inverse_restores_exact_basis=True,
        scope='Pure bitwise identity/refusal preselection in cached_multiply only. '
              'Original general-lane partition, identity precedence, unknown/legacy metadata, '
              'valid words, signed zeros, inactive lanes, fallback and FP expressions/order remain. '
              'No new ABI, ISA, call, persistent allocation, raw stage, host or shader change. '
              'Source reasoning only; no speed or scientific/release claim.'))
elif phase in ('build', 'link', 'run'):
    prepared, built = verified_inputs()
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    proof = doc(WORK / 'derivation.json')
    assert pin(WORK / 'retained8.cpp') == {
        key: proof['candidate'][key] for key in ('bytes', 'sha256')}
    inputs = set(prepared['inputs']) | {
        str(path.relative_to(ROOT)) for path in WORK.iterdir() if path.is_file()}
    if phase == 'build':
        assert not BUILD.exists()
        BUILD.mkdir()
        argv = list(prepared['compilations'][0]['argv'])
        argv[argv.index('-c') + 1] = str(WORK / 'retained8.cpp')
        argv[argv.index('-o') + 1] = str(BUILD / 'retained8.o')
        scope = 'One strict original-flags C++ helper compilation for exact bitwise preselection. No workload.'
    else:
        previous = doc(WORK / ('build' if phase == 'link' else 'link') / 'owner.json')
        assert previous['passed'] and previous['status'] == 'terminal'
        assert previous['whole_inputs_unchanged'] and previous['source_unchanged']
        assert not previous['cleanup_errors'] and not previous['remaining_owned_processes']
        driver = next(row['object'] for row in built['compilations'] if row['name'] == 'driver')
        inputs |= {driver['path'], str((BUILD / 'retained8.o').relative_to(ROOT))}
        if phase == 'link':
            argv = [prepared['link'][0], str(ROOT / driver['path']),
                    str(BUILD / 'retained8.o'), '-o', str(BUILD / 'abi-check')]
            scope = 'Link unchanged original scalar/alias/inactive/tag driver object with the one changed helper.'
        else:
            inputs.add(str((BUILD / 'abi-check').relative_to(ROOT)))
            argv = [str(BUILD / 'abi-check')]
            scope = 'One original unchanged focused scalar compatibility fixture:555660 calls, '
            scope += '66679200 bank words/13335840 metadata words, masks/aliases/signed zeros/refusals/FP controls. '
            scope += 'Supporting finite compatibility only; no independent math, provider, benefit, frame or release qualification.'
    actions = doc(WORK / 'actions.json') if (WORK / 'actions.json').exists() else {}
    assert phase not in actions
    actions[phase] = dict(revision=revision, kind='diagnostic', expected_exit=0,
                          limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
                          argv=argv, frozen_inputs=sorted(inputs), software_provider=False,
                          scope=scope)
    (WORK / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=phase, executed=False, revision=revision)))
elif phase == 'inspect':
    owner = doc(WORK / 'build/owner.json')
    assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    prepared, built = verified_inputs()
    baseline = ROOT / next(row['object']['path'] for row in built['compilations'] if row['name'] == 'retained8')
    def inspect(path):
        raw = subprocess.check_output(['objdump', '-d', '--insn-width=16', str(path)], text=True)
        functions = {}
        name = None
        mnemonics = set()
        for line in raw.splitlines():
            match = re.match(r'^[0-9a-f]+ <(.+)>:', line)
            if match:
                name = match[1]
                functions[name] = []
            fields = line.split('\t')
            if len(fields) < 3 or not re.match(r'^\s*[0-9a-f]+:$', fields[0]):
                continue
            encoded = fields[1].split()
            assert encoded[0] != '62', 'Unexpected EVEX'
            text = fields[2].strip()
            if not text:
                continue
            mnemonic = text.split()[0]
            mnemonics.add(mnemonic)
            if name:
                functions[name].append(text)
        assert not re.search(r'%(?:zmm\d+|k[0-7]|mm[0-7]|st)\b', raw)
        assert not any(m.startswith(('vfm', 'vfnm', 'f')) for m in mnemonics)
        return dict(identity=pin(path), functions={
            name: dict(instructions=len(rows),
                       stack_memory_instructions=sum(bool(re.search(r'\([^)]*%(?:rsp|rbp)[^)]*\)', row)) for row in rows),
                       prologue=rows[:14]) for name, rows in functions.items()})
    result = dict(baseline=inspect(baseline), candidate=inspect(BUILD / 'retained8.o'),
                  scope='Static compiled code only. Stack-operand instructions are not measured spills, latency, dynamic costs or speedup.')
    write(WORK / 'code-inspection.json', result)
    print(json.dumps({kind: {name: {key: value for key, value in row.items() if key != 'prologue'}
                           for name, row in result[kind]['functions'].items()
                           if name in ('sirius_native_tagged_multiply8', 'sirius_native_tagged_divide8')}
                      for kind in ('baseline', 'candidate')}, indent=2))
elif phase == 'read':
    owner = doc(WORK / 'run/owner.json')
    assert owner['passed'] and owner['status'] == 'terminal'
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
    result = doc(WORK / 'run/stdout.log')
    expected = dict(passed=True, calls=555660, compared_words=66679200,
                    tag_words=13335840, active_lane_results=1666980,
                    case_tuples=49, input_alias_and_inactive_words_checked=True,
                    FP_controls_unchanged=True)
    assert result == expected, result
    write(WORK / 'compatibility.json', dict(result=result, executable=pin(BUILD / 'abi-check'),
          object=pin(BUILD / 'retained8.o'), driver_recompiled=False,
          scope=owner['scope']))
    print(json.dumps(result))
else:
    raise ValueError(phase)
