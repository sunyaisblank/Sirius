"""Retire three fixed preserved FP64-estimate scratch roots only.

Default creates a reviewable audit. --remove requires its exact reviewed hash
and repeats preservation, terminal-owner, declared-input and local-consumer
checks. No compiler/test/provider invocation, Git write, tag or push.
"""
from pathlib import Path
import argparse
import fnmatch
import json
import os
import re
import shlex
import shutil
import stat

import mature_cleanup as mature

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CONFIG = json.loads((HERE / 'configuration.json').read_text())
TARGETS = tuple(ROOT / name for name in CONFIG['targets'])
RECEIPT = ROOT / 'attestations/cleanup/2026-10-15/fp64-estimate-terminal-scratch'
require, identity = mature.require, mature.identity


def load(path):
    return json.loads(Path(path).read_text())


def inside(path):
    return any(path == target or path.is_relative_to(target) for target in TARGETS)


def source_join(expected):
    mature.clean_source()
    require(expected == CONFIG['current_product_source'], 'exact clean4498 source required')
    require(mature.git('rev-list', '--left-right', '--count', 'HEAD...@{upstream}').split() == [b'0', b'0'],
            'canonical branch is not synchronized with its upstream')
    branch = mature.git('symbolic-ref', '--short', 'HEAD').decode().strip()
    require(branch == CONFIG['branch'], 'exact identified development branch required')
    remote_ref = mature.git('rev-parse', '--symbolic-full-name', '@{upstream}').decode().strip()
    require(remote_ref.startswith('refs/remotes/origin/'), 'origin tracking branch required')
    remote_branch = 'refs/heads/' + remote_ref.removeprefix('refs/remotes/origin/')
    require(mature.git('ls-remote', 'origin', remote_branch).decode().split() == [expected, remote_branch],
            'actual source branch is not pushed at the exact expected commit')
    worktrees = mature.git('worktree', 'list', '--porcelain').decode().splitlines()
    require([line for line in worktrees if line.startswith('worktree ')] == ['worktree ' + str(ROOT)],
            'sole canonical worktree required')
    require(not mature.git('ls-files', '-z', '--', *CONFIG['targets']), 'a cleanup root contains tracked files')
    return {'revision': expected, 'clean': True, 'branch': branch, 'upstream_0_0': True,
            'exact_origin_branch_verified': remote_branch, 'sole_worktree': str(ROOT)}


def snapshots():
    result = {}
    for target in TARGETS:
        files = mature.files_under(target)
        directories = {'.': stat.S_IMODE(target.stat().st_mode)}
        for base, children, _ in os.walk(target, followlinks=False):
            for child in children:
                path = Path(base) / child
                directories[path.relative_to(target).as_posix()] = stat.S_IMODE(path.stat().st_mode)
        rows = {name: dict(identity(path), mode=stat.S_IMODE(path.stat().st_mode)) for name, path in files.items()}
        result[str(target.relative_to(ROOT))] = {'files': rows, 'directories': directories,
            'regular_files': len(rows), 'regular_bytes': sum(row['bytes'] for row in rows.values())}
    return result


def preservation():
    # Each complete archive retains the original raw/Git/remote/manifest check.
    # Scratch is matched below against their union, not prematurely per archive.
    driver = CONFIG['driver_archive']
    rows = [*CONFIG['prior_evidence'], [driver['tag'], driver['commit'], driver['bundle'], driver['manifest_sha256'], []]]
    evidence = [mature.preserved(row) for row in rows]
    versions = {}
    for row in rows:
        bundle = ROOT / row[2]
        for member, expected in load(bundle / 'manifest.json')['files'].items():
            for target in TARGETS:
                prefix = 'files/' + str(target.relative_to(ROOT)) + '/'
                if member.startswith(prefix):
                    path = str(target.relative_to(ROOT)) + '/' + member[len(prefix):]
                    versions.setdefault(path, []).append({'commit': row[1], 'member': member,
                        'bytes': expected['bytes'], 'sha256': expected['sha256'],
                        'mode': stat.S_IMODE((bundle / member).stat().st_mode)})
    owned = {}
    for target in TARGETS:
        files = mature.files_under(target)
        prefix = str(target.relative_to(ROOT)) + '/'
        require(set(files) == {name[len(prefix):] for name in versions if name.startswith(prefix)},
                'complete fixed-root union file inventory differs: ' + str(target))
        for name, path in files.items():
            entries = versions[prefix + name]
            final = [entry for entry in entries if entry['commit'] == driver['commit']]
            chosen = final[0] if final else entries[0]
            if not final:
                require(all({key: entry[key] for key in ('bytes', 'sha256', 'mode')} ==
                            {key: chosen[key] for key in ('bytes', 'sha256', 'mode')} for entry in entries),
                        'conflicting prior versions need an authoritative final archive member')
            actual = dict(identity(path), mode=stat.S_IMODE(path.stat().st_mode))
            require(actual == {key: chosen[key] for key in ('bytes', 'sha256', 'mode')},
                    'owned bytes/mode differ from preserved union: ' + str(path))
            owned[prefix + name] = {'identity_and_mode': actual, 'authoritative_commit': chosen['commit'],
                'preserved_versions': entries}
    bundle = ROOT / driver['bundle']
    manifest = load(bundle / 'manifest.json')['files']
    records = {}
    for role, member in driver['required_records'].items():
        require(member in manifest, 'required final review is not archived: ' + role)
        records[role] = load(bundle / member)
    require(records['capture_acceptance']['passed'] is True
            and records['capture_acceptance']['original_failed_compiler_wrapper_preserved'] is True,
            'original completed capture or failed compiler history is missing')
    require(records['independent_terminal_review']['passed'] is True
            and records['independent_terminal_review']['source_revision'] == CONFIG['current_product_source'],
            'actual independent terminal review is missing')
    require(records['summary']['source_revision'] == CONFIG['current_product_source']
            and records['summary']['qualification_claimed'] is False
            and records['cost_review'], 'source/qualification scope or completed cost review missing')
    # A rejected cost premise remains a valid terminal cleanup disposition.
    # Preserve the actual cost review; do not reinterpret it as numerical proof.
    stop = ROOT / CONFIG['prelaunch_stop']['root']
    finding = load(stop / 'finding.json')
    require(identity(stop / 'finding.json') == CONFIG['prelaunch_stop']['finding_pin']
            and finding['operator_or_helper_execution'] is False
            and finding['source_or_owner_failure_log_created'] is False
            and not list(stop.rglob('owner.json')), 'unconsumed prelaunch stop was relabeled/executed')
    stop_files = mature.files_under(stop)
    require(set(stop_files) == set(CONFIG['prelaunch_stop']['files']) and len(stop_files) == 10,
            'exact ten prelaunch stop files changed')
    for name, expected in finding['preimages'].items():
        require(identity(stop / name) == expected, 'seven-file prelaunch preimage changed')
    mandatory_final = [name for name in owned if name.startswith(CONFIG['prelaunch_stop']['root'] + '/')
        or name.startswith(CONFIG['driver_root'] + '/') or name.startswith(CONFIG['driver_bin_root'] + '/')]
    mandatory_final += CONFIG['final_top_level_files']
    require(all('files/' + name in manifest for name in mandatory_final),
            'final driver archive lacks whole driver/bin/prelaunch/current release coverage')
    return {'evidence': evidence, 'owned_union': owned,
        'final_review_members': driver['required_records'], 'unconsumed_prelaunch_files': 10,
        'mode_scope': 'Current owned file modes match preserved raw copy2 file modes and are captured in the audit/receipt. Generic collector Git100644 is not claimed as original executable-mode provenance.'}


def terminal_records():
    records = mature.retired_owners()
    require({row['path'] for row in records} == set(CONFIG['expected_owners']) and len(records) == 5,
            'the exact five reviewed owners changed')
    births = []
    for row in records:
        path = ROOT / row['path']
        owner = load(path)
        expected = CONFIG['expected_owners'][row['path']]
        require(owner['passed'] is expected['passed'] and owner['child_exit'] == expected['child_exit']
                and owner['source'] == {'revision': CONFIG['current_product_source'], 'status': ''}
                and owner['stop_reason'] is None and len(owner['observed_births']) == expected['births'],
                'original passing/failed owner disposition or birth count changed')
        require(owner['wall_seconds'] < owner['limits']['seconds']
                and owner['peak_sampled_rss_kib'] <= owner['limits']['rss_kib']
                and owner['minimum_sampled_available_kib'] >= owner['limits']['reserve_kib'],
                'original owner limits were exceeded')
        for birth in owner['observed_births']:
            births.append((owner.get('boot_id'), birth['pid'], str(birth['saved_start_ticks'])))
        before, after = path.parent / 'inputs-before.json', path.parent / 'inputs-after.json'
        require(load(before) == load(after), 'owner whole-input records differ')
        row.update(inputs_before=identity(before), inputs_after=identity(after), recorded_births=len(owner['observed_births']))
    require(len(births) == len(set(births)) == 9, 'exact nine distinct recorded births required')
    return records


def declared_inputs():
    """Fresh bounded scan of actual out/ declaration paths, not prose history.

    This is the existing frozen-input/argv/pin-path check expanded to discover a
    newly prepared screen. It is not a general Python dependency proof.
    """
    def option_path(value):
        value = value.split('=', 1)[-1]
        for prefix in ('-isystem', '-include', '-I', '-L'):
            if value.startswith(prefix):
                return value[len(prefix):]
        return value

    def path_leaf(value):
        if not isinstance(value, str):
            return None
        value = option_path(value)
        return value if value.startswith('/') or re.match(r'^(?:\./)*(?:out|bin)/', value) else None

    files, visited = [], 0
    for base, directories, names in os.walk(ROOT / 'out', followlinks=False):
        base = Path(base)
        for name in directories:
            candidate = base / name
            require(not candidate.is_symlink() or not inside(candidate.resolve()),
                    'external directory alias points into a fixed retirement root')
        directories[:] = [name for name in directories if not inside(base / name)
                          and base / name != HERE and not (base / name).is_symlink()]
        visited += len(directories) + len(names)
        require(visited <= 50000, 'out/ declaration census exceeds its explicit entry bound')
        for name in names:
            if any(fnmatch.fnmatchcase(name, pattern) for pattern in CONFIG['metadata_name_patterns']):
                path = base / name
                require(not path.is_symlink() and path.is_file() and path.stat().st_size <= 20 * 1024 * 1024,
                        'unexpected declaration file type/size')
                files.append(path)
    require(len(files) <= 200, 'declaration file count exceeds its explicit bound')
    records = {}
    for path in sorted(files):
        value = load(path)
        tokens = []

        def walk(node):
            if isinstance(node, dict):
                for name, child in node.items():
                    if path_leaf(name):
                        tokens.append(path_leaf(name))
                    if name in {'pins', 'files', 'private_inputs', 'external_inputs'} and isinstance(child, dict):
                        tokens.extend(child)
                    if name in {'frozen_inputs', 'argv', 'commands'}:
                        def strings(item):
                            if isinstance(item, str):
                                tokens.extend(shlex.split(item))
                            elif isinstance(item, list):
                                for element in item:
                                    strings(element)
                        strings(child)
                    if name in {'path', 'resolved_path', 'source_freeze', 'preparation', 'actions_file'} and isinstance(child, str):
                        tokens.append(child)
                    walk(child)
            elif isinstance(node, list):
                for child in node:
                    walk(child)
            elif path_leaf(node):
                tokens.append(path_leaf(node))
        walk(value)
        for token in tokens:
            if not isinstance(token, str) or not token:
                continue
            token = option_path(token)
            if not token:
                continue
            candidate = Path(token)
            candidates = [candidate] if candidate.is_absolute() else [ROOT / candidate, path.parent / candidate]
            require(not any(inside(item) or inside(item.resolve()) for item in candidates),
                    'remaining packet consumes retired scratch: ' + str(path) + ': ' + token)
        records[str(path.relative_to(ROOT))] = identity(path)
    return {'declaration_files': records, 'metadata_files_examined': len(records),
        'live_cleanup_root_path_found': False,
        'scope': 'Fresh selected out/ metadata frozen-input, pin/path, scalar path leaves, path-like map keys and actual recipe/argv paths including joined -I/-L/-include/-isystem options; excludes the three fixed retired roots and this cleanup apparatus. No history exemption or general undeclared program-dependency proof.',
        'next_screen_rule': 'Use retained canonical builds and immutable d547/56076/3dcc/final-driver archives only; any external future input into the fixed three roots blocks retirement.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-source', required=True)
    parser.add_argument('--remove', action='store_true')
    parser.add_argument('--accepted-audit-sha256')
    args = parser.parse_args()
    fixed = ['out/portable-divide-fp64-estimate-screen',
        'bin/linux-gcc/portable-divide-fp64-estimate-screen',
        'out/portable-divide-fp64-estimate-property-query']
    require(CONFIG['released'] is True and CONFIG['targets'] == fixed
            and args.expected_source == CONFIG['current_product_source'] == '4498cc41168ad29fc4c111819a4cb7ccaa54a8bf',
            'held candidate is not a final exact three-root release')
    driver = CONFIG['driver_archive']
    require(isinstance(driver['tag'], str) and driver['tag'].startswith('evidence/issue-44-')
            and isinstance(driver['bundle'], str) and driver['bundle'].startswith('attestations/software-vulkan/4bcd987/')
            and re.fullmatch('[0-9a-f]{40}', driver['commit'] or '')
            and re.fullmatch('[0-9a-f]{64}', driver['manifest_sha256'] or '')
            and set(driver['required_records']) == {'capture_acceptance', 'independent_terminal_review', 'cost_review', 'summary'}
            and all(isinstance(member, str) and member for member in driver['required_records'].values()),
            'actual preserved/pushed driver archive and review member coordinates required')
    freeze = load(HERE / 'source-freeze.json')
    require(freeze['schema'] == 'sirius-fp64-estimate-retirement-final-source-freeze-v1'
            and freeze['released'] is True, 'final independently reviewable source freeze required')
    require({str((HERE / name).relative_to(ROOT)) for name in
        ['retire.py', 'mature_cleanup.py', 'configuration.json', 'commands.json']} <= set(freeze['pins']),
        'required final cleanup source is not pinned')
    for name, expected in freeze['pins'].items():
        require(identity(ROOT / name) == expected, 'cleanup preparation source changed')
    mature.ROOT, mature.SOURCE, mature.TARGETS = ROOT, args.expected_source, TARGETS
    mature.prefix_target = lambda path, targets=TARGETS: any(path == item or path.is_relative_to(item) for item in targets)
    source = source_join(args.expected_source)
    previous = None
    if args.remove:
        previous = load(HERE / 'audit.json')
        require(identity(HERE / 'audit.json')['sha256'] == args.accepted_audit_sha256
                and previous['status'] == 'verified-not-removed' and previous['source'] == source
                and previous['driver_archive'] == driver, 'exact independently reviewed audit required')
        require(not RECEIPT.exists(), 'fresh fixed cleanup receipt root required')
    else:
        require(args.accepted_audit_sha256 is None and not (HERE / 'audit.json').exists(), 'fresh audit required')
    proof = preservation()
    terminal = terminal_records()
    dependencies = declared_inputs()
    roots = snapshots()
    census = mature.consumer_census()
    document = {'schema': 'sirius-fp64-estimate-fixed-terminal-retirement-v1', 'status': 'verified-not-removed',
        'source': source, 'driver_archive': driver, 'preservation': proof, 'roots': roots,
        'terminal_owners': terminal, 'recorded_owner_count': 5, 'recorded_birth_count': 9,
        'historical_pre_canary_inventory_counts_not_reused': True, 'declared_inputs': dependencies,
        'initial_census': census, 'removed': [], 'qualification_claimed': False,
        'global_consumer_absence_claimed': False, 'canonical_and_other_builds_retained': True}
    if not args.remove:
        mature.write(HERE / 'audit.json', document)
        print(json.dumps({'audit': identity(HERE / 'audit.json'), 'fixed_roots': CONFIG['targets'],
            'terminal_owners': 5, 'recorded_births': 9, 'permission_gaps': len(census['permission_gaps'])}))
        return
    require(previous['preservation'] == proof and previous['roots'] == roots
            and previous['declared_inputs'] == dependencies, 'preserved roots/declarations changed after audit review')
    source_join(args.expected_source)
    document['final_census'] = mature.consumer_census()
    require(snapshots() == roots, 'owned bytes/modes/layout changed immediately before fixed removal')
    RECEIPT.mkdir(parents=True)
    for name in ['retire.py', 'mature_cleanup.py', 'configuration.json', 'source-freeze.json', 'commands.json', 'audit.json',
                 'ownership-dependencies.json', 'retire.preimage.py', 'adaptation.patch', 'source-lineage.json',
                 'inventory-source-freeze.json', 'initial-unconsumed-release.json']:
        shutil.copy2(HERE / name, RECEIPT / name)
    mature.write(RECEIPT / 'receipt.json', document)
    for target in TARGETS:
        require(target in TARGETS and target.is_dir() and not target.is_symlink() and target.resolve() == target,
                'fixed owned root was replaced')
        shutil.rmtree(target)
        require(not target.exists(), 'fixed root still exists')
        document['removed'].append(str(target.relative_to(ROOT)))
        mature.write(RECEIPT / 'receipt.json', document)
    source_join(args.expected_source)
    document['status'] = 'preserved-terminal-roots-removed'
    document['source_unchanged'] = True
    mature.write(RECEIPT / 'receipt.json', document)
    mature.write(RECEIPT / 'manifest.json', {'files': {path.name: identity(path) for path in sorted(RECEIPT.iterdir()) if path.is_file()}})
    print(json.dumps({'receipt': identity(RECEIPT / 'receipt.json'), 'removed': document['removed'],
        'Git_commit_tag_push_performed': False, 'canonical_and_other_builds_retained': True}))


if __name__ == '__main__':
    main()
