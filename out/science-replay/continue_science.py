"""Continue independent cohorts after a whole-cohort P4 time guard stopped scheduling."""

from pathlib import Path
import shutil
import subprocess
import sys

from run_science import (GUARD, MODEL, ROOT, WORK, collect_plan, pin, read, require,
                         verify_junit, write)


def main():
    result_path = WORK / 'result.json'
    result = read(result_path)
    initial = read(WORK / 'runs/actions.json')
    p4 = result['cohorts']['p4-host']
    require(p4['status']['status'] == 'terminal' and p4['status']['stop_reason'] == 'diagnostic_timeout',
            'Continuation requires the recorded terminal P4 timeout')
    require(not p4['status']['remaining_owned_processes'] and not p4['status']['cleanup_errors'],
            'Continuation requires a drained P4 session')
    require(not subprocess.check_output(['git', 'status', '--porcelain']), 'Source tree changed')
    require(subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip() == result['revision'],
            'Selected source revision changed')
    require(shutil.disk_usage('/').free >= 100_000_000_000, 'Less than 100 GB available on /')

    groups, _ = collect_plan(read(MODEL))
    original = ['p6-host', 'p6-kernel', 'p4-kernel', 'original-connected-render', 'original-64x64-writers']
    attempts = {label: (label, initial[label]) for label in original}
    recovery = dict(initial['p4-host'])
    recovery['limits'] = {**recovery['limits'], 'seconds': 900}
    render = Path(recovery['build_environment']['TMPDIR']).with_name('p4-host-recovered')
    render.mkdir()
    recovery['build_environment'] = {**recovery['build_environment'], 'TMPDIR': str(render)}
    attempts['p4-host-recovered'] = ('p4-host', recovery)

    work = WORK / 'continued'
    work.mkdir()
    shutil.copyfile(GUARD / 'run_owned.py', work / 'run_owned.py')
    (work / 'owned_guard.py').symlink_to(GUARD / 'owned_guard.py')
    for _, action in attempts.values():
        action['frozen_inputs'] = [*action['frozen_inputs'], str(Path(__file__).resolve()),
                                   str(WORK / 'runs/actions.json')]
    write(work / 'actions.json', {label: action for label, (_, action) in attempts.items()})

    for label, (origin, action) in attempts.items():
        require(label not in result['cohorts'], f'Existing outcome for {label}')
        folder = work / label
        code = subprocess.run([sys.executable, '-B', str(work / 'run_owned.py'), label]).returncode
        owner_path = folder / 'owner.json'
        owner = read(owner_path)
        cohort = {'owner': pin(owner_path), 'stdout': pin(folder / 'stdout.log'),
                  'stderr': pin(folder / 'stderr.log'), 'selection_from': origin,
                  'status': {key: owner[key] for key in ('passed', 'status', 'wall_seconds',
                             'peak_sampled_rss_kib', 'limits', 'stop_reason', 'child_exit', 'source_unchanged',
                             'whole_inputs_unchanged', 'cleanup_errors', 'remaining_owned_processes')}}
        junit = folder / 'junit.xml'
        if junit.is_file():
            cohort['junit'] = pin(junit)
        result['cohorts'][label] = cohort
        write(result_path, result)
        require(owner['status'] == 'terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
                and not owner['remaining_owned_processes'] and not owner['cleanup_errors'],
                f'Invalid terminal state in {label}')
        if code or not owner['passed']:
            print(f'{label}: preserved unsuccessful result; other independent cohorts can proceed', flush=True)
            continue
        try:
            cohort['cases'] = verify_junit(junit, groups[origin][0])
            if action['software_provider']:
                require(owner['provider_mappings'], f'Selected software Vulkan mapping missing in {label}')
                cohort['provider'] = pin(Path(action['runtime_provider']))
                cohort['icd'] = pin(Path(action['runtime_icd']))
        except (OSError, RuntimeError) as exc:
            cohort['capture_error'] = str(exc)
        else:
            cohort['accepted'] = True
        write(result_path, result)

    return all(result['cohorts'][label].get('accepted') for label in attempts)


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
