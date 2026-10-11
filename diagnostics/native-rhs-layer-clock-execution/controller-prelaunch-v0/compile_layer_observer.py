"""Root-owned bounded existing-toolchain admission; no Vulkan or native execution."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import ast, hashlib, json, os, re, signal, subprocess, time
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
PRIOR = ROOT / 'attestations/native-vulkan/2053ef7/native-shader-clock-admission/files'
COMPILER = Path('/opt/slang/bin/slangc')
LIBRARY = Path('/opt/slang/lib/libslang-compiler.so.0.2026.12.0.1')
TOOLS = ROOT / 'bin/linux-gcc/qualification-tools/spirv-tools-v2026.3/build/tools'
BASELINE = '2053ef79c90752d37b0075b8137cfb9b571b6820'
def read(p):
    return json.loads(p.read_text())
def seal(p):
    return {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
def birth(pid):
    try:
        return Path('/proc', str(pid), 'stat').read_text().rsplit(')', 1)[1].split()[19]
    except (FileNotFoundError, ProcessLookupError):
        return None
def git(*args):
    return subprocess.check_output(['git', *args], text=True, timeout=30).strip()
assert git('rev-parse', 'HEAD') == BASELINE == git('rev-parse', '@{u}')
assert not git('status', '--porcelain')
assert read(WORK / 'source-review.json')['accepted']
assert read(WORK / 'compiler-controller-review.json')['accepted']
prior = read(PRIOR / 'phase-compiler.json')
assert prior['accepted']
preparation = read(WORK / 'preparation.json')
assert preparation['recovered_prior_observer_byte_exact'] and preparation['row_words'] == 8226
for name, expected in preparation['input_seals'].items():
    assert seal(Path(name)) == expected
for name, expected in preparation['output'].items():
    assert seal(Path(name)) == expected
paths = [Path(__file__), WORK / 'prepare_layer_observer.py', WORK / 'phase/retained_transport.slang',
         WORK / 'preparation.json', WORK / 'source-review.json', WORK / 'compiler-controller-review.json',
         PRIOR / 'phase-compiler.json', PRIOR / 'clock-features.json', COMPILER, LIBRARY,
         TOOLS / 'spirv-val', TOOLS / 'spirv-dis', TOOLS / 'spirv-as', ROOT / 'scripts/build-retained-kernels.py']
paths.extend((ROOT / 'src/sirius/kernels').glob('retained_*.slang'))
for p in (COMPILER, LIBRARY, TOOLS / 'spirv-val', TOOLS / 'spirv-dis', TOOLS / 'spirv-as'):
    assert seal(p) == prior['input_seals'][str(p)]
features = read(PRIOR / 'clock-features.json')
assert features['clock_features']['shaderDeviceClock'] and features['logical_devices_created'] == 0
record = {'accepted': False, 'baseline_revision': BASELINE, 'input_seals': {str(p): seal(p) for p in paths},
          'seconds_per_command': 60, 'rss_bytes': 536870912, 'runs': [], 'errors': [],
          'scope': 'Native-default Transport layer observer compiler/validator admission only. No logical device/pipeline/dispatch. Inclusive perturbed lane0 source cuts; no exclusive arithmetic/barrier/access cost or ordinary rendering benefit/qualification.'}
signals = []
for s in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
    signal.signal(s, lambda n, f: signals.append(n))
def save():
    record['signals'] = list(signals)
    (WORK / 'compiler.json').write_text(json.dumps(record, indent=2) + '\n')
def run(name, args):
    r = {'name': name, 'arguments': list(map(str, args)), 'errors': [], 'outer_stop': False, 'peak_rss_bytes': 0}
    child = None
    ticks = None
    started = time.monotonic()
    try:
        assert not signals
        with (WORK / (name + '.stdout')).open('xb') as out, (WORK / (name + '.stderr')).open('xb') as err:
            child = subprocess.Popen(r['arguments'], stdout=out, stderr=err, start_new_session=True)
            r.update(pid=child.pid, start_ticks=None)
            ticks = birth(child.pid)
            r['start_ticks'] = ticks
            assert ticks, 'Child birth unavailable'
            while child.poll() is None:
                try:
                    rss = int(Path('/proc', str(child.pid), 'statm').read_text().split()[1]) * os.sysconf('SC_PAGE_SIZE')
                    r['peak_rss_bytes'] = max(r['peak_rss_bytes'], rss)
                except (FileNotFoundError, ProcessLookupError):
                    pass
                if signals or time.monotonic() - started > 60 or r['peak_rss_bytes'] > 536870912:
                    raise RuntimeError('Fixed signal/time/RSS guard triggered')
                try:
                    child.wait(timeout=.05)
                except subprocess.TimeoutExpired:
                    pass
            assert child.returncode == 0
    except BaseException as e:
        r['errors'].append(repr(e))
        r['outer_stop'] = True
    finally:
        if child is not None:
            try:
                if child.poll() is None:
                    assert ticks and birth(child.pid) == ticks and os.getpgid(child.pid) == child.pid
                    os.killpg(child.pid, signal.SIGKILL)
            except BaseException as e:
                r['errors'].append('Owned stop: ' + repr(e))
            try:
                r['returncode'] = child.wait(timeout=10)
            except BaseException as e:
                r['errors'].append('Owned reap: ' + repr(e))
            r['owned_birth_absent'] = bool(ticks) and birth(child.pid) != ticks
            r['owned_group_absent'] = False
            try:
                os.killpg(child.pid, 0)
            except ProcessLookupError:
                r['owned_group_absent'] = True
            except BaseException as e:
                r['errors'].append('Group census: ' + repr(e))
        r['elapsed_seconds'] = time.monotonic() - started
        r['accepted'] = (not r['errors'] and not signals and r.get('returncode') == 0
                         and r.get('owned_birth_absent') and r.get('owned_group_absent'))
        record['runs'].append(r)
        save()
    assert r['accepted'], r
try:
    args = list(prior['reused_transport_compiler_run']['arguments'])
    args[1] = str(WORK / 'phase/retained_transport.slang')
    args[args.index('-I') + 1] = str(ROOT / 'src/sirius/kernels')
    raw = WORK / 'phase/retained_transport.compiler.spv'
    assembly = WORK / 'phase/retained_transport.spvasm'
    final = WORK / 'phase/retained_transport.spv'
    assert not any(p.exists() for p in (raw, assembly, final))
    args[-1] = str(raw)
    run('transport-slang', args)
    run('transport-disassemble', [TOOLS / 'spirv-dis', raw, '-o', assembly])
    builder = ROOT / 'scripts/build-retained-kernels.py'
    tree = ast.parse(builder.read_text())
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'native_controls')
    env = {'re': re}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(builder), 'exec'), env)
    text = env['native_controls'](assembly.read_text())
    assembly.write_text(text)
    run('transport-assemble', [TOOLS / 'spirv-as', '--target-env', 'spv1.5', assembly, '-o', final])
    run('transport-validate', [TOOLS / 'spirv-val', '--target-env', 'vulkan1.2', final])
    assert text.count('OpReadClockKHR') == 2
    assert 'OpCapability ShaderClockKHR' in text and 'OpExtension "SPV_KHR_shader_clock"' in text
    assert 'OpTypeInt 64' not in text and 'OpCapability Int64' not in text
    assert re.search(r'OpExecutionMode %\S+ LocalSize 64 1 1', text)
    assert re.search(r'OpDecorate %phaseObservations Binding 2', text)
    calls = {}
    for name, expected in (('ObservePhase', 4), ('ObserveLayer', 3)):
        symbol = re.search(r'OpName (%\S+) "' + name + '"', text)[1]
        calls[name] = len(re.findall(r'OpFunctionCall %\S+ ' + re.escape(symbol) + r'\b', text))
        assert calls[name] == expected
    unsigned32 = re.search(r'(%\S+) = OpTypeInt 32 0\s*$', text, re.MULTILINE)[1]
    devices = set(re.findall(r'(%\S+) = OpConstant ' + re.escape(unsigned32) + r' 1\s*$', text, re.MULTILINE))
    assert devices and all(line.split()[-1] in devices for line in text.splitlines() if 'OpReadClockKHR' in line)
    record['product'] = {'spv': seal(final), 'assembly': seal(assembly), 'clock_instruction_count': 2,
                         'helper_calls': calls, 'clock_scope': 'Device(1)', 'logical_binding': 2,
                         'row_words': 8226, 'maximum_layers': 117, 'layer_entry_words': 10}
    assert {str(p): seal(p) for p in paths} == record['input_seals']
    assert git('rev-parse', 'HEAD') == BASELINE and not git('status', '--porcelain')
    assert not signals
    record['accepted'] = True
except BaseException as e:
    record['errors'].append(repr(e))
finally:
    save()
print(json.dumps({k: record[k] for k in ('accepted', 'runs', 'errors')}))
assert record['accepted']
