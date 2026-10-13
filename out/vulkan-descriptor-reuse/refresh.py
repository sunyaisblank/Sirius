from pathlib import Path
import hashlib, json, shlex, shutil, subprocess

R = Path.cwd()
W = R / 'out/vulkan-descriptor-reuse'
T = R / 'bin/linux-gcc/vulkan-descriptor-reuse'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert subprocess.check_output(['git', 'status', '--porcelain'], text=True) == ''
if not (W / 'actions-e75f156.json').exists():
    shutil.copyfile(W / 'actions.json', W / 'actions-e75f156.json')
else:
    assert (W / 'actions-e75f156.json').read_bytes() == (W / 'actions.json').read_bytes()
old_products = W / 'prior-products'
old_products.mkdir(exist_ok=True)
for suffix in ['o', 'd']:
    prior = old_products / ('vulkan_smoke_test-e75f156.' + suffix)
    if prior.exists():
        assert prior.read_bytes() == (T / ('vulkan_smoke_test.' + suffix)).read_bytes()
    else:
        shutil.copyfile(T / ('vulkan_smoke_test.' + suffix), prior)
reuse = []
recompile = {}
for name in ['vulkan_device', 'retained_compute_test', 'retained_dopri_test']:
    owner_path = W / ('compile-' + name) / 'owner.json'
    owner = json.loads(owner_path.read_text())
    assert owner['passed'] and owner['status'] == 'terminal'
    before = json.loads((owner_path.parent / 'inputs-before.json').read_text())
    dependency_path = T / (name + '.d')
    words = shlex.split(dependency_path.read_text().replace('\\\n', '').split(':', 1)[1])
    dependencies = {}
    missing_freeze = []
    for word in words:
        p = Path(word).resolve()
        key = str(p.relative_to(R)) if p.is_relative_to(R) else str(p)
        identity = {'bytes': p.stat().st_size, 'sha256': sha(p)}
        if key not in before:
            missing_freeze.append(str(p))
        else:
            assert all(before[key][k] == v for k, v in identity.items()), (key, identity)
        dependencies[key] = identity
    obj = T / (name + '.o')
    if missing_freeze:
        for suffix in ['o', 'd']:
            shutil.copyfile(T / (name + '.' + suffix), old_products / (name + '-e75f156.' + suffix))
        recompile[name] = {'missing_precompile_freeze': missing_freeze, 'actual_dependencies': list(dependencies)}
        continue
    reuse.append({'object': {'path': str(obj.relative_to(R)), 'bytes': obj.stat().st_size, 'sha256': sha(obj)},
                  'actual_owner': str(owner_path.relative_to(R)), 'owner_sha256': sha(owner_path),
                  'actual_source_revision': owner['source']['revision'],
                  'dependency_file_sha256': sha(dependency_path),
                  'byte_identical_current_non_system_dependencies': dependencies})
actions = json.loads((W / 'actions.json').read_text())
actions['compile-layout-witness'] = dict(actions['compile-vulkan_smoke_test'], revision=revision)
for name, record in recompile.items():
    action = dict(actions['compile-' + name], revision=revision)
    action['frozen_inputs'] = list(dict.fromkeys(action['frozen_inputs'] + record['actual_dependencies']))
    action['scope'] = 'Recompile with complete actual non-system dependency freeze after discovering generated inputs missing from the earlier private compiler action.'
    actions['compile-complete-' + name] = action
for key, action in actions.items():
    if not key.startswith('compile-'): action['revision'] = revision
(W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
(W / 'current-object-reuse.json').write_text(json.dumps({'candidate_revision': revision, 'objects': reuse, 'recompile_for_complete_freeze': recompile,
    'scope': 'Reuse exactly unchanged non-system compiler dependencies; rebuild only the changed smoke TU. All immutable prior source/action versions remain separately recoverable.'}, indent=2) + '\n')
print(json.dumps({'candidate_revision': revision, 'unchanged_objects_reused': len(reuse), 'changed_TUs': ['vulkan_smoke_test'], 'incomplete_freeze_TUs': list(recompile)}))
