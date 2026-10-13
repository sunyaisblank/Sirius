from pathlib import Path
import hashlib, json, shutil, subprocess, sys

R = Path.cwd()
W = Path(__file__).resolve().parent
T = R / 'bin/linux-gcc/vulkan-descriptor-reuse'
A = R / 'attestations/software-vulkan/8a27321/vulkan-descriptor-reuse'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
identity = lambda p: {'bytes': p.stat().st_size, 'sha256': sha(p)}
load = lambda p: json.loads(p.read_text())
write = lambda p, x: p.write_text(json.dumps(x, indent=2) + '\n')
relative = lambda p: str(p.relative_to(R))
candidate = load(W / 'audited-benefit/retention-gate-frozen.json')['candidate_revision']
assert len(sys.argv) == 2 and sys.argv[1] in ['stage', 'finalise']
assert not (A / 'manifest.json').exists()

if sys.argv[1] == 'stage':
    assert not A.exists()
    owners = sorted(W.rglob('owner.json'))
    assert len(owners) == 31
    for p in owners:
        o = load(p)
        assert o['status'] == 'terminal' and o['passed'] and o['child_exit'] == 0
        assert not o['cleanup_errors'] and not o['stop_reason']
        assert o['remaining_owned_processes'] == []
        assert all(b['absent'] for b in o['observed_births'])
        assert o['whole_inputs_unchanged'] and o['source_unchanged']
        assert load(p.parent / 'inputs-before.json') == load(p.parent / 'inputs-after.json')
    assert load(W / 'control-acceptance.json')['passed']
    assert load(W / 'audited-witness-acceptance.json')['passed']
    assert not (W / 'benefit/software-comparison.json').exists()
    benefit = load(W / 'audited-benefit/software-comparison.json')
    assert len(benefit['checks']) == 20 and len(benefit['accounting_checks']) == 2
    source = load(W / 'preserved-source.json')
    assert source['parent'] == candidate
    assert subprocess.check_output(['git', 'ls-remote', 'origin', 'refs/tags/' + source['tag']], text=True).strip() == source['remote_readback']
    A.mkdir(parents=True)
    shutil.copytree(W, A / 'files' / relative(W))
    raw = {tuple(identity(p).values()): p for p in A.rglob('*') if p.is_file()}
    prior = {}
    manifests = subprocess.check_output(['rg', '--files', 'attestations', '--glob', 'manifest.json'], text=True).splitlines()
    for name in manifests:
        path = R / name
        items = load(path).get('files', [])
        if not isinstance(items, list):
            continue
        for record in items:
            if isinstance(record, dict) and all(k in record for k in ['path', 'bytes', 'sha256']):
                target = (path.parent / record['path']).resolve()
                if target.is_relative_to(path.parent.resolve()):
                    prior.setdefault((record['bytes'], record['sha256']), target)
    def archive_recovery(p, expected):
        key = (expected['bytes'], expected['sha256'])
        if key in raw:
            target = raw[key]
            assert identity(target) == expected
        elif key in prior and prior[key].is_file() and identity(prior[key]) == expected:
            target = prior[key]
        else:
            assert p.is_file() and identity(p) == expected, str(p)
            target = A / 'input-versions' / expected['sha256']
            target.parent.mkdir(exist_ok=True)
            if not target.exists():
                shutil.copyfile(p, target)
            assert identity(target) == expected
        raw[key] = target
        return {'kind': 'archive', 'path': relative(target)}
    products = []
    for p in sorted(T.rglob('*')):
        if p.is_file():
            value = identity(p)
            products.append({'origin': relative(p), **value, 'recovery': archive_recovery(p, value)})
    inputs = {}
    git_cache = {}
    for owner_path in owners:
        owner = load(owner_path)
        revision = owner['source']['revision']
        if revision not in git_cache:
            output = subprocess.check_output(['git', 'ls-tree', '-r', '-z', revision])
            git_cache[revision] = {row.split(b'\t', 1)[1].decode(): row.split(b'\t', 1)[0].split()[2].decode()
                                   for row in output.split(b'\0') if row}
        for origin, value in load(owner_path.parent / 'inputs-before.json').items():
            expected = {k: value[k] for k in ['bytes', 'sha256']}
            key = (origin, expected['bytes'], expected['sha256'])
            if key in inputs:
                inputs[key]['consumers'].append(relative(owner_path))
                continue
            if origin in git_cache[revision]:
                blob = git_cache[revision][origin]
                data = subprocess.check_output(['git', 'cat-file', 'blob', blob])
                assert len(data) == expected['bytes'] and hashlib.sha256(data).hexdigest() == expected['sha256']
                recovery = {'kind': 'git', 'revision': revision, 'blob': blob, 'path': origin}
            else:
                p = R / origin
                # External tool/provider files are shared installed resources,
                # retained rather than removed with either disposable root.
                if not p.is_relative_to(R):
                    assert p.is_file() and identity(p) == expected
                    recovery = {'kind': 'retained_external', 'path': str(p)}
                else:
                    recovery = archive_recovery(p, expected)
            inputs[key] = {'origin': origin, **expected, 'recovery': recovery,
                           'consumers': [relative(owner_path)]}
    write(A / 'input-version-recovery.json', {'records': list(inputs.values()),
        'scope': 'Every observed input/version joined to immutable Git or exact archived bytes; external installed tool/provider versions checked and retained. Historical incomplete and corrected compiler receipts remain distinct.'})
    write(A / 'private-product-recovery.json', {'records': products,
        'scope': 'Every complete private-root file recoverable byte-exactly; identical previously archived products reused without duplicate payloads.'})
    write(A / 'disposition.json', {'candidate_revision': candidate, 'controls_pass': True,
        'performance_pass': benefit['pass'], 'original20_pass': benefit['original20_pass'],
        'inclusive_api_cost_pass': benefit['provisioning_accounting_pass'],
        'owners': len(owners), 'source_tag': source['tag'], 'source_commit': source['commit'],
        'initial_collection': 'Four normal raw executions retained but missing command/API timing schema; no20+2 verdict. Only corrected audited-benefit determines adoption.',
        'scope': 'Finite descriptor binding candidate; correctness prerequisites and actual frozen cost verdict only. No complete frame, full scientific, native/runtime or release qualification claim.'})
    print(json.dumps({'staged': relative(A), 'owners': len(owners), 'input_versions': len(inputs), 'private_products': len(products)}))
else:
    assert (A / 'disposition.json').is_file()
    for p in sorted(W.rglob('*')):
        if not p.is_file():
            continue
        target = A / 'files' / relative(p)
        if target.exists():
            assert target.read_bytes() == p.read_bytes(), str(p)
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(p, target)
    for ledger in ['input-version-recovery.json', 'private-product-recovery.json']:
        for record in load(A / ledger)['records']:
            recovery = record['recovery']
            if recovery['kind'] == 'git':
                data = subprocess.check_output(['git', 'cat-file', 'blob', recovery['blob']])
                assert len(data) == record['bytes'] and hashlib.sha256(data).hexdigest() == record['sha256']
            else:
                assert identity(R / recovery['path']) == {k: record[k] for k in ['bytes', 'sha256']}
    files = [{'path': str(p.relative_to(A)), **identity(p)} for p in sorted(A.rglob('*')) if p.is_file()]
    write(A / 'manifest.json', {'source_revision': candidate, 'files': files,
        'passed': load(A / 'disposition.json')['performance_pass'],
        'scope': 'Sealed complete actual #97 candidate evidence, including adverse cost results, historical compiler versions and exact recovery joins. Adoptability follows unchanged original20+2 gates.'})
    print(json.dumps({'sealed': relative(A), 'files': len(files), 'bytes': sum(p['bytes'] for p in files), 'manifest_sha256': sha(A / 'manifest.json')}))
