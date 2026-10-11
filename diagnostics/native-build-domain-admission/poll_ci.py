"""Observe the one exact native-build evidence run without launching work."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import datetime
import hashlib
import json
import subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
source = json.loads((WORK / 'source.json').read_text())
run = json.loads((WORK / 'ci-run.json').read_text())
assert run['source_revision'] == source['source_revision']
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True, timeout=30).strip() == source['source_revision']
assert not subprocess.check_output(['git', 'status', '--porcelain'], timeout=30)
base = f"repos/sunyaisblank/Sirius/actions/runs/{run['run_id']}"

def fetch(endpoint):
    return json.loads(subprocess.check_output(['gh', 'api', endpoint], timeout=40))

with ThreadPoolExecutor(max_workers=2) as pool:
    meta_future = pool.submit(fetch, base)
    jobs_future = pool.submit(fetch, base + '/jobs?per_page=100')
    meta, jobs = meta_future.result(), jobs_future.result()
assert meta['id'] == run['run_id'] and meta['head_sha'] == source['source_revision']
assert meta['event'] == 'workflow_dispatch' and meta['head_branch'] == source['branch']
assert jobs['total_count'] == len(jobs['jobs']) == 7
expected = {'integration-no-render', 'integration-windows-no-render', 'integration-macos-no-render'}
skipped = {'linux-gate', 'linux-sanitizers', 'windows-build', 'macos-build'}
assert {j['name'] for j in jobs['jobs']} == expected | skipped
record = {
    'source_revision': source['source_revision'], 'run_id': run['run_id'],
    'observed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'status': meta['status'], 'conclusion': meta['conclusion'],
    'workflow_dispatch_request_sha256': hashlib.sha256((WORK / 'workflow-dispatch-request.json').read_bytes()).hexdigest(),
    'metadata': meta, 'jobs': jobs,
    'scope': 'Authoritative observation of the one exact hosted run. Source/build/non-render evidence only; no full Mandatory, native runtime, frame or release qualification.',
}
if meta['status'] == 'completed':
    record['accepted'] = meta['conclusion'] == 'success' and all(
        j['status'] == 'completed' and j['conclusion'] == ('success' if j['name'] in expected else 'skipped')
        for j in jobs['jobs']
    )
stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%f')
destination = WORK / f'ci-observation-{stamp}.json'
assert not destination.exists()
destination.write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({
    'run_id': run['run_id'], 'status': meta['status'], 'conclusion': meta['conclusion'],
    'accepted': record.get('accepted'), 'receipt': str(destination.relative_to(ROOT)),
    'jobs': [{'name': j['name'], 'status': j['status'], 'conclusion': j['conclusion'],
              'current_step': next((s['name'] for s in j['steps'] if s['status'] == 'in_progress'), None)}
             for j in jobs['jobs']],
}))
