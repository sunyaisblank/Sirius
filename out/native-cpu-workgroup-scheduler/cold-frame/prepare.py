"""Adapt the sealed original single-link/single-frame authority to provider E."""
from pathlib import Path
import hashlib

basis = Path.cwd() / 'out/native-cpu-row-ownership/cold-frame/prepare.py'
assert hashlib.sha256(basis.read_bytes()).hexdigest() == '5c702abf6582dba1a79b47628576db51c96bd6544f8bac46136f562de7f34f6f'
source = basis.read_text()


def change(before, after):
    global source
    assert source.count(before) == 1, before
    source = source.replace(before, after)


change('P = W.parent', "P = R / 'out/native-cpu-workgroup-scheduler'")
change("B = C / 'native-cpu-row-ownership'", "B = C / 'native-cpu-workgroup-scheduler'")
change("FRAME = B / 'frame'; SOURCE = 'bfa095baf0ff9592cf648cfd49a73d5f6814c05c'",
       "FRAME = B / 'frame'\nARCHIVE = C / 'native-cpu-row-ownership/matched/libsirius_backend.a'\n"
       "SOURCE = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()")
source = source.replace("B / 'matched/libsirius_backend.a'", 'ARCHIVE')
for before, after in [
    ("P / 'readers/controls-1-8-acceptance.json'", "P / 'controls/readers/controls-1-8-acceptance.json'"),
    ("P / 'readers/benefit.json'", "P / 'benefit/readers/benefit.json'"),
    ("P / 'readers/decimal-agreement.json'", "P / 'benefit/readers/decimal-agreement.json'"),
    ("P / 'providers.json'", "P / 'benefit/providers.json'"),
    ("P / 'runtime/control-8/inputs-before.json'", "P / 'control-8/inputs-before.json'"),
    ("P / 'compiler-result.json'", "R / 'out/native-cpu-row-ownership/compiler-result.json'"),
    ("P / 'matched-result.json'", "R / 'out/native-cpu-row-ownership/matched-result.json'"),
]:
    assert before in source, before
    source = source.replace(before, after)

gates = '''
benefit = doc(P / 'benefit/readers/benefit.json')
decimal = doc(P / 'benefit/readers/decimal-agreement.json')
assert benefit['audit_passed'] and benefit['all22_numeric_pass'] and benefit['checks_passed'] == 22
assert decimal['agrees_on_all22'] and decimal['all22_pass']
c8 = doc(P / 'control-8-acceptance.json')
assert c8['passed'] and c8['all32_readbacks_match_preserved_originals']
assert c8['native_first_packets_handled'] == 3
observed = doc(P / 'control-8/native-observation.json')['calls']
assert len(observed) == 3 and all(x['groups'] == 1 and x['handled'] for x in observed)
assert doc(P / 'provider/linked-native-proof.json')['passed']
compiled = doc(R / 'out/native-cpu-static-program/batch-compile-result.json')
assert len(compiled['objects']) == 3 and all(x['exit'] == 0 for x in compiled['objects'])
provider = doc(P / 'benefit/providers.json')['providers']['candidate']
for key in ('library', 'ICD'):
    assert pin(Path(provider[key])) == provider[key + '_identity']
assert pin(ARCHIVE) == {'bytes': 6871532, 'sha256': '07693aba054a519f1127220bcb3afd94e42d10f416c7c6358fda336af046f215'}
host = C / 'native-cpu-row-ownership/matched/audited_backend_tests'
assert pin(host)['sha256'] == '7f2d1fa88953b954ab116233458a4cb962470f863fda9bbb3a1cad7d4a7d6a26'
saved_recipe = R / 'out/native-original-input-ownership/cold-frame/link-recipe.json'
assert pin(saved_recipe)['sha256'] == '101e87ba20534ce8911f2e861721ea34ec69ca690315f6faa13c60c12bdcf137'
evidence_inputs = {R / name for name in doc(P / 'controls/runtime/actions.json')['control-7']['frozen_inputs']}
evidence_inputs |= {host, saved_recipe, Path('/usr/bin/g++-14'),
    P / 'controls/readers/controls-1-8-acceptance.json', P / 'benefit/readers/benefit.json',
    P / 'benefit/readers/decimal-agreement.json', P / 'benefit/providers.json',
    P / 'control-8-acceptance.json', P / 'control-8/native-observation.json',
    P / 'provider/linked-native-proof.json', R / 'out/native-cpu-static-program/batch-compile-result.json',
    R / 'out/native-cpu-row-ownership/cold-frame/prepare.py',
    R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
    R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py',
    R / 'out/native-cpu-arithmetic-sharing/readers/native_evidence.py'}
evidence_inputs |= {p for p in W.iterdir() if p.is_file()}
for i in range(1, 8):
    for name in ('owner.json', 'gtest.xml', 'inputs-before.json', 'inputs-after.json'):
        evidence_inputs.add(P / 'controls/runtime' / ('control-' + str(i)) / name)
assert all(p.is_file() for p in evidence_inputs)

def module_joins(executable):
    data = executable.read_bytes()
    modules = doc(R / 'out/native-cpu-row-ownership/compiler-result.json')['modules']
    result = {}
    for kind in ('Transport', 'Endpoint', 'DopriPhase'):
        item = modules[kind]
        path = R / item['module']
        assert pin(path) == {k: item[k] for k in ('bytes', 'sha256')}
        assert data.count(path.read_bytes()) == 1, kind
        result[kind] = {'module': item['module'], **pin(path), 'complete_occurrences': 1}
    return result

'''
change('mode, = sys.argv[1:]\n', gates + 'mode, = sys.argv[1:]\n')
start = source.index("    lines = subprocess.check_output(['ninja'")
end = source.index('    original_exe = C / target', start)
source = source[:start] + '''    saved = doc(saved_recipe)
    args = list(saved['argv'])
    previous_archive = str(C / 'native-original-input-ownership/matched/libsirius_backend.a')
    assert args.count(previous_archive) == 1
    args = [str(ARCHIVE) if x == previous_archive else x for x in args]
    output_index = args.index('-o') + 1
    original_inputs = [Path(x) if Path(x).is_absolute() else C / x
        for i, x in enumerate(args) if i not in (0, output_index) and not x.startswith('-')]
    assert len(original_inputs) == 40 and all(p.is_file() for p in original_inputs)
''' + source[end:]
change("    assert args.count('src/sirius/backend/libsirius_backend.a') == 1\n"
       "    args = [str(ARCHIVE) if x == 'src/sirius/backend/libsirius_backend.a' else x for x in args]\n",
       '    assert args.count(str(ARCHIVE)) == 1\n')
change("    volume = C / 'src/sirius/app/resources'\n", "    volume = Path(saved['original_resources'])\n")
change('    assert resources and all(p.resolve().is_relative_to(volume.resolve()) for p in volume.rglob(\'*\') if p.is_file())\n',
       "    assert resources == saved['resources']\n"
       "    assert resources and all(p.resolve().is_relative_to(volume.resolve()) for p in volume.rglob('*') if p.is_file())\n")
change("    inputs |= {p for p in volume.rglob('*') if p.is_file()}\n",
       "    inputs |= {p for p in volume.rglob('*') if p.is_file()}\n    inputs |= evidence_inputs\n")
change("    inputs |= {p for p in W.iterdir() if p.is_file()}\n",
       "    inputs |= {p for p in W.iterdir() if p.is_file()}\n    inputs |= evidence_inputs\n")
change("elif mode == 'build':\n",
       "elif mode == 'build':\n    assert not (W / 'build-result.json').exists()\n"
       "    assert not (FRAME / 'sirius_render_tests').exists()\n")
change("    write(W / 'build-result.json', {'passed': True, 'source': SOURCE,\n",
       "    joins = module_joins(FRAME / 'sirius_render_tests')\n"
       "    write(W / 'build-result.json', {'passed': True, 'source': SOURCE, 'module_joins': joins,\n")
change("    temp = R / 'renders/native-cpu-row-bfa095ba/cold-frame-tmp'\n",
       "    assert module_joins(exe) == doc(W / 'build-result.json')['module_joins']\n"
       "    temp = R / ('renders/native-cpu-workgroup-' + SOURCE[:8] + '/cold-frame-tmp')\n")
exec(compile(source, str(basis), 'exec'))
