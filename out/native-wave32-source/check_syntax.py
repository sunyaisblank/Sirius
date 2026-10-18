import concurrent.futures
import hashlib
import json
import os
import resource
import shlex
import signal
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
rows = json.loads((ROOT / 'bin/linux-gcc/compile_commands.json').read_text())
recipe = next(row for row in rows if row['file'].endswith('/vulkan/vulkan_device.cpp'))
arguments = shlex.split(recipe['command'])
index = arguments.index('-o')
del arguments[index:index + 2]
arguments.remove('-c')
arguments[arguments.index(recipe['file'])] = str(WORK / 'vulkan_device.cpp')
arguments.append('-fsyntax-only')


def pin(path):
    path = Path(path)
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


inputs = [WORK / 'vulkan_device.cpp', ROOT / recipe['file'],
          ROOT / 'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h',
          Path(arguments[0]), ROOT / '.clang-format',
          ROOT / 'src/sirius/backend/vulkan/vulkan_device.h',
          ROOT / 'src/sirius/backend/vulkan/vulkan_portability.h',
          ROOT / 'bin/linux-gcc/compile_commands.json']
before = {str(p): pin(p) for p in inputs}


def bounds():
    resource.setrlimit(resource.RLIMIT_AS, (4 * 1024 ** 3, 4 * 1024 ** 3))


def run(mode):
    argv = arguments.copy()
    if mode == 'without-retained':
        argv.remove('-DSIRIUS_HAS_RETAINED_COMPUTE=1')
        argv.remove('-I' + str(ROOT / 'bin/linux-gcc/src/sirius/backend/retained'))
    started = time.monotonic()
    with (WORK / (mode + '.stdout')).open('wb') as out, (WORK / (mode + '.stderr')).open('wb') as err:
        process = subprocess.Popen(argv, cwd=recipe['directory'], stdout=out, stderr=err,
                                   start_new_session=True, preexec_fn=bounds)
        timed_out = False
        try:
            code = process.wait(timeout=60)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process.pid, signal.SIGKILL)
            code = process.wait()
    return {'mode': mode, 'argv': argv, 'returncode': code, 'timeout': timed_out,
            'elapsed_seconds': time.monotonic() - started,
            'address_space_limit_per_process_bytes': 4 * 1024 ** 3,
            'deadline_seconds': 60}


with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    results = list(pool.map(run, ['with-retained', 'without-retained']))
after = {str(p): pin(p) for p in inputs}
assert before == after
record = {'source_revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
          'inputs_before': before, 'inputs_after': after, 'results': results,
          'scope': 'Strict existing GCC preset recipe, syntax only; two supported source modes, no build product mutation/runtime/scientific/release qualification'}
record['passed'] = all(row['returncode'] == 0 and not row['timeout'] for row in results)
(WORK / 'syntax.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'passed': record['passed'], 'results': results}))
raise SystemExit(0 if record['passed'] else 1)
