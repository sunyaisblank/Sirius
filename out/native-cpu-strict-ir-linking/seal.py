"""Seal this compiler gate and retire only recoverable SDK download fragments."""
from pathlib import Path
from collections import Counter
import hashlib
import json
import os
import re
import stat
import subprocess
import sys
import zipfile

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-strict-ir-linking'
A = R / 'attestations/qualification/97b3a761/native-cpu-strict-ir-linking-compiler-gate'
PHASES = ('emit', 'acquire', 'resume', 'verify', 'plugin', 'convert',
          'plugin-sha256', 'convert-sha256', 'math', 'link', 'optimize')


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,
                    sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def owners():
    result = {}
    for phase in PHASES:
        path = W / phase / 'owner.json'
        owner = json.loads(path.read_text())
        assert owner['status'] == 'terminal' and owner['source_unchanged']
        assert owner['whole_inputs_unchanged'] and not owner['remaining_owned_processes']
        assert not owner['cleanup_errors']
        for birth in owner['observed_births']:
            proc = Path('/proc') / str(birth['pid']) / 'stat'
            try:
                current = proc.read_text().rsplit(')', 1)[1].split()[19]
            except FileNotFoundError:
                continue
            assert current != birth['saved_start_ticks'], (phase, birth)
        result[phase] = dict(path=str(path.relative_to(R)), **pin(path),
            passed=owner['passed'], source=owner['source'], seconds=owner['wall_seconds'],
            peak_rss_kib=owner['peak_sampled_rss_kib'], stop_reason=owner['stop_reason'])
    assert [n for n in PHASES if not result[n]['passed']] == ['acquire', 'convert']
    return result


phase, = sys.argv[1:]
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
history = owners()
if phase == 'seal':
    assert not A.exists()
    before = json.loads((W / 'link-result.json').read_text())
    after = json.loads((W / 'optimize-result.json').read_text())
    counts = {}
    for label, result in (('linked', before), ('optimized', after)):
        text = (R / result['ir']['path']).read_text()
        assert pin(R / result['ir']['path']) == {k: result['ir'][k] for k in ('bytes', 'sha256')}
        counts[label] = dict(definitions=result['definitions'],
            private_math_calls=sum(result['private_math_calls'].values()),
            packed_payload_allocas=len(re.findall(r'alloca \[120 x i32\]', text)),
            packed_tag_allocas=len(re.findall(r'alloca \[24 x i32\]', text)),
            context_allocas=text.count('alloca %v8_varying_KernelContext_0'),
            ordinary_fp_operations=result['ordinary_fp_operations'])
    assert counts['linked']['private_math_calls'] == 1581
    assert counts['optimized']['private_math_calls'] == 1580
    assert counts['linked']['packed_payload_allocas'] == counts['optimized']['packed_payload_allocas'] == 1581
    assert counts['linked']['packed_tag_allocas'] == counts['optimized']['packed_tag_allocas'] == 1575
    download = json.loads((W / 'download-result.json').read_text())
    archive = R / download['archive']['path']
    assert pin(archive) == {k: download['archive'][k] for k in ('bytes', 'sha256')}
    fragments = []
    offset = 0
    for key in ('immutable_prefix', 'suffix'):
        entry = download[key]
        path = R / entry['path']
        assert pin(path) == {k: entry[k] for k in ('bytes', 'sha256')}
        metadata = path.lstat()
        assert stat.S_ISREG(metadata.st_mode) and metadata.st_uid == os.getuid() and metadata.st_nlink == 1
        fragments.append(dict(**entry, archive_offset=offset, uid=metadata.st_uid,
            device=metadata.st_dev, inode=metadata.st_ino, nlink=1,
            allocated_bytes=metadata.st_blocks * 512))
        offset += entry['bytes']
    assert offset == download['archive']['bytes']
    write(W / 'retirement-plan.json', dict(archive=download['archive'], fragments=fragments,
        recovery='Read exact archive_offset/bytes from the complete verified local SDK archive; compare fragment SHA256 before reuse. Historical frozen-input paths are retired only after this mapping is sealed and pushed. SDK and complete archive stay cached.'))
    paths = [p for p in W.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    paths += [p for p in B.iterdir() if p.is_file() and p.suffix in ('.bc', '.ll', '.so')]
    paths += [R / 'out/native-cpu-exact-product-zero-residual/retained8.cpp',
              R / 'out/native-cpu-exact-product-zero-residual/optional.cpp',
              R / 'out/native-cpu-inline-private-bank/packinghelper.isph',
              R / 'out/native-cpu-value-fact-stores/value_facts.isph',
              R / 'out/native-cpu-local-identities/identity_store.isph',
              R / 'bin/linux-gcc/native-cpu-inline-private-bank/ispc-source/dopri_phase.ispc',
              R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8/src/util/u_math.h',
              R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8/src/gallium/drivers/llvmpipe/lp_state_cs.c',
              R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
              R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py']
    entries = {str(p.relative_to(R)): pin(p) for p in sorted(set(paths))}
    A.mkdir(parents=True)
    bundle = A / 'compiler-gate.zip'
    with zipfile.ZipFile(bundle, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zipped:
        for name in entries:
            zipped.write(R / name, name)
    with zipfile.ZipFile(bundle) as zipped:
        assert set(zipped.namelist()) == set(entries)
        for name, identity in entries.items():
            with zipped.open(name) as stream:
                assert hashlib.file_digest(stream, 'sha256').hexdigest() == identity['sha256'], name
            assert zipped.getinfo(name).file_size == identity['bytes']
            assert pin(R / name) == identity, name
    manifest = dict(status='COMPILER_PIPELINE_ESTABLISHED_PACKED_ESCAPE_REMAINS',
        execution_source_revision='97b3a7611d3bb2dadb297fcfcfc246e6b2869d9e',
        custody_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        owners=history, structure=counts,
        decision='Ordinary O2 internalization propagates fixed indices but retains1580/1581 math calls and all1581 payload banks. Stop before native codegen, numerical fixtures, provider rebuild or renderer. A material interface change is required; no speed or equivalence inferred.',
        bundle=dict(path=str(bundle.relative_to(R)), **pin(bundle), members=entries),
        setup=dict(official_sdk=json.loads((W / 'toolchain.json').read_text()),
                   initial_network_failure_retained=True, failed_plugin_load_retained=True,
                   setup_and_failure_costs_preserved=True),
        cleanup_plan=dict(path=str((W / 'retirement-plan.json').relative_to(R)), **pin(W / 'retirement-plan.json')),
        limitations='Compiler/verifier/static IR gates only. No new scientific, numerical, frame, hardware-domain, public-provider or release acceptance. Original scene/counts/layout/allocations/arithmetic/work limits and old public objects unchanged.')
    write(A / 'MANIFEST.json', manifest)
    print(json.dumps(dict(sealed=str(A.relative_to(R)), bundle_bytes=bundle.stat().st_size,
                          phase_passes=Counter(row['passed'] for row in history.values()), structure=counts)))
elif phase == 'retire':
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD']) == subprocess.check_output(['git', 'rev-parse', '@{upstream}'])
    subprocess.run(['git', 'ls-files', '--error-unmatch', str((A / 'MANIFEST.json').relative_to(R)),
                    str((A / 'compiler-gate.zip').relative_to(R))], check=True, stdout=subprocess.DEVNULL)
    manifest = json.loads((A / 'MANIFEST.json').read_text())
    bundle = manifest['bundle']
    assert pin(R / bundle['path']) == {k: bundle[k] for k in ('bytes', 'sha256')}
    plan = json.loads((W / 'retirement-plan.json').read_text())
    assert pin(W / 'retirement-plan.json') == {k: manifest['cleanup_plan'][k] for k in ('bytes', 'sha256')}
    archive = R / plan['archive']['path']
    assert pin(archive) == {k: plan['archive'][k] for k in ('bytes', 'sha256')}
    for fragment in plan['fragments']:
        with archive.open('rb') as stream:
            stream.seek(fragment['archive_offset'])
            digest = hashlib.sha256()
            remaining = fragment['bytes']
            while remaining:
                data = stream.read(min(1024 * 1024, remaining))
                assert data
                digest.update(data)
                remaining -= len(data)
            assert digest.hexdigest() == fragment['sha256']
    targets = {str((R / f['path']).resolve()) for f in plan['fragments']}
    holders = []
    denied = []
    for proc in Path('/proc').glob('[0-9]*'):
        try:
            for fd in (proc / 'fd').iterdir():
                try:
                    destination = os.readlink(fd)
                except FileNotFoundError:
                    continue
                if destination in targets:
                    holders.append(dict(pid=proc.name, fd=fd.name, path=destination))
            for line in (proc / 'maps').read_text().splitlines():
                if any(line.endswith(' ' + target) for target in targets):
                    holders.append(dict(pid=proc.name, mapping=line))
        except FileNotFoundError:
            continue
        except PermissionError:
            denied.append(proc.name)
    assert not holders, holders
    free_before = os.statvfs(R).f_bavail * os.statvfs(R).f_frsize
    allocated = 0
    for fragment in plan['fragments']:
        path = R / fragment['path']
        current = path.lstat()
        assert stat.S_ISREG(current.st_mode) and current.st_nlink == fragment['nlink'] == 1
        assert (current.st_uid, current.st_dev, current.st_ino) == (fragment['uid'], fragment['device'], fragment['inode'])
        assert pin(path) == {k: fragment[k] for k in ('bytes', 'sha256')}
        allocated += current.st_blocks * 512
        path.unlink()
    free_after = os.statvfs(R).f_bavail * os.statvfs(R).f_frsize
    write(A / 'retirement.json', dict(retired=plan['fragments'],
        allocated_bytes_unlinked=allocated, volume_free_before=free_before,
        volume_free_after=free_after, observed_volume_free_delta=free_after - free_before,
        open_consumer_matches=holders, inaccessible_processes=denied,
        recoverability='Every fragment rehashed from exact full-archive range before unlink; sealed manifest and full SDK/archive remain.',
        retained_archive=dict(path=str(archive.relative_to(R)), **pin(archive)),
        cache_retained='Verified12.4GB expanded SDK and full2.014GB archive; no global install or sibling workspace.',
        all_owned_phases_terminal=True, useful_sources_and_bundle_pushed_before_retirement=True))
    print(json.dumps(dict(retired_files=len(plan['fragments']), allocated_bytes_unlinked=allocated,
                          observed_volume_free_delta=free_after - free_before)))
else:
    raise ValueError(phase)
