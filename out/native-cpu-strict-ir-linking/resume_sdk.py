"""Continue verified-length network work without discarding received bytes."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-strict-ir-linking'
ASSET = 'LLVM-23.1.2-Linux-X64.tar.xz'
PREFIX = B / (ASSET + '.part')
TAIL = B / (ASSET + '.tail.part')
ARCHIVE = B / ASSET
TOTAL = 2013958032


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


phase, = sys.argv[1:]
source = json.loads((W / 'toolchain-source.json').read_text())
if phase in ('prepare-resume', 'prepare-verify'):
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    actions = json.loads((W / 'actions.json').read_text())
    case = phase.removeprefix('prepare-')
    assert case not in actions and not (W / case).exists()
    inputs = [*W.glob('*.py'), W / 'toolchain-source.json']
    if case == 'resume':
        failed = json.loads((W / 'acquire/owner.json').read_text())
        assert failed['status'] == 'terminal' and not failed['passed']
        assert failed['source_unchanged'] and failed['whole_inputs_unchanged']
        assert not failed['remaining_owned_processes'] and not TAIL.exists() and not ARCHIVE.exists()
        assert PREFIX.stat().st_size == 988196919
        write(W / 'resume-plan.json', dict(prefix=dict(path=str(PREFIX.relative_to(R)), **pin(PREFIX)),
            expected_total_bytes=TOTAL, requested_range=f'{PREFIX.stat().st_size}-{TOTAL-1}',
            source=source, prior_failed_acquisition_owner=dict(path=str((W/'acquire/owner.json').relative_to(R)), **pin(W/'acquire/owner.json')),
            reason='Original fixed network time limit interrupted acquisition; retain immutable received prefix and fetch only the missing byte range, with no whole-download retry.'))
        inputs += [PREFIX, W / 'resume-plan.json', Path('/usr/bin/curl').resolve()]
        argv = ['/usr/bin/python3', '-B', str(W / 'resume_sdk.py'), 'resume']
    else:
        owner = json.loads((W / 'resume/owner.json').read_text())
        assert owner['status'] == 'terminal' and owner['passed'] and owner['whole_inputs_unchanged']
        assert pin(ARCHIVE)['sha256'] == source['sha256'] and not (B / 'toolchain').exists()
        inputs += [ARCHIVE, W / 'download-result.json', Path('/usr/bin/ldd').resolve()]
        argv = ['/usr/bin/python3', '-B', str(W / 'verify_sdk.py')]
    actions[case] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='diagnostic', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152), argv=argv,
        frozen_inputs=[str(p.relative_to(R)) if p.is_relative_to(R) else str(p) for p in inputs],
        scope='Explicit continuation of SDK setup after a preserved network-timeout receipt. Immutable received prefix, exact requested suffix and full published checksum; extraction/version checks remain separate. No hidden retry, setup-cost subtraction, physics/native execution, global installation or governor change.')
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=case, executed=False)))
elif phase == 'resume':
    plan = json.loads((W / 'resume-plan.json').read_text())
    assert pin(PREFIX) == {k: plan['prefix'][k] for k in ('bytes', 'sha256')}
    with TAIL.open('xb') as stream:
        subprocess.run(['/usr/bin/curl', '--fail', '--location', '--silent', '--show-error',
                        '--retry', '0', '--max-time', '280', '--proto', '=https',
                        '--range', plan['requested_range'], source['url']], stdout=stream, check=True)
    assert TAIL.stat().st_size == TOTAL - PREFIX.stat().st_size, 'Server did not provide exactly the missing range'
    with ARCHIVE.open('xb') as destination:
        for part in (PREFIX, TAIL):
            with part.open('rb') as stream:
                shutil.copyfileobj(stream, destination, 8 * 1024 * 1024)
    assert pin(ARCHIVE) == dict(bytes=TOTAL, sha256=source['sha256']), 'Complete published checksum mismatch'
    write(W / 'download-result.json', dict(passed=True, archive=dict(path=str(ARCHIVE.relative_to(R)), **pin(ARCHIVE)),
        immutable_prefix=plan['prefix'], suffix=dict(path=str(TAIL.relative_to(R)), **pin(TAIL)),
        whole_download_retries=0, initial_acquisition_failed=True, published_full_checksum_verified=True))
    print(json.dumps(dict(download_complete=True, published_checksum_verified=True, bytes=TOTAL)))
else:
    raise ValueError(phase)
