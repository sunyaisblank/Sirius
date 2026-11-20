"""Run the supplemental P4 kernels with ordinary, task-local Mesa disk caching."""

from pathlib import Path
import shutil
import subprocess
import sys

from run_science import (GUARD, ICD, MODEL, PROVIDER, WORK, collect_plan, pin, read, require,
                         verify_junit, write)


def main():
    original = WORK / 'continued'
    prior = read(original / 'p4-kernel/owner.json')
    require(prior['status'] == 'terminal' and prior['stop_reason'] == 'diagnostic_timeout' and
            not prior['cleanup_errors'] and not prior['remaining_owned_processes'],
            'Expected the terminal supplemental cold batch timeout')
    require(shutil.disk_usage('/').free >= 100_000_000_000, 'Less than 100 GB available on /')
    groups, _ = collect_plan(read(MODEL))

    folder = WORK / 'p4-kernel-recovery'
    folder.mkdir()
    cache = folder / 'shader-cache'
    cache.mkdir()
    label = 'p4-kernel-cached'
    shutil.copyfile(GUARD / 'run_owned.py', folder / 'run_owned.py')
    (folder / 'owned_guard.py').symlink_to(GUARD / 'owned_guard.py')

    action = read(original / 'actions.json')['p4-kernel']
    action['disable_software_shader_cache'] = False
    env = action['build_environment']
    env.pop('MESA_SHADER_CACHE_DISABLE')
    env['MESA_SHADER_CACHE_DIR'] = str(cache)
    action['limits'] = {**action['limits'], 'seconds': 1200}
    action['frozen_inputs'].append(str(Path(__file__).resolve()))
    action['scope'] = 'Original registered supplemental P4 shader workloads with task-local Mesa disk caching'
    write(folder / 'actions.json', {label: action})
    code = subprocess.run([sys.executable, '-B', str(folder / 'run_owned.py'), label]).returncode

    run = folder / label
    owner_path = run / 'owner.json'
    owner = read(owner_path)
    junit = run / 'junit.xml'
    result = {'owner': pin(owner_path), 'stdout': pin(run / 'stdout.log'),
              'stderr': pin(run / 'stderr.log'), 'selection_from': 'p4-kernel',
              'status': {key: owner[key] for key in ('passed', 'status', 'wall_seconds',
                         'peak_sampled_rss_kib', 'limits', 'stop_reason', 'child_exit', 'source_unchanged',
                         'whole_inputs_unchanged', 'cleanup_errors', 'remaining_owned_processes')}}
    if junit.is_file():
        result['junit'] = pin(junit)
    if not code and owner['source_unchanged'] and owner['whole_inputs_unchanged'] and owner['provider_mappings']:
        result['cases'] = verify_junit(junit, groups['p4-kernel'][0])
        result['provider'] = pin(PROVIDER)
        result['icd'] = pin(ICD)
        result['accepted'] = True
    write(folder / 'result.json', result)
    return result.get('accepted', False)


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
