"""Retire six fixed preserved scratch roots, never any bin/ build.

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
RECEIPT = ROOT / 'attestations/cleanup/2026-10-15/minimum-step-floor-terminal-scratch'
require, identity = mature.require, mature.identity


def load(path):
    return json.loads(Path(path).read_text())


def inside(path):
    return any(path == target or path.is_relative_to(target) for target in TARGETS)


def source_join(expected):
    mature.clean_source()
    require(mature.git('merge-base', '--is-ancestor', CONFIG['current_product_source'], expected) == b'',
            'final source does not retain the integrated focused repair')
    require(mature.git('rev-list', '--left-right', '--count', 'HEAD...@{upstream}').split() == [b'0', b'0'],
            'canonical branch is not synchronized with its upstream')
    branch = mature.git('symbolic-ref', '--short', 'HEAD').decode().strip()
    require(branch.startswith('development/'), 'identified development branch required')
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


def preservation(commit, manifest):
    focused = CONFIG['focused']
    rows = [*CONFIG['prior_evidence'], [focused['tag'], commit, focused['bundle'], manifest, focused['targets']]]
    require({name for row in rows for name in row[4]} == set(CONFIG['targets']), 'fixed preservation set differs')
    result = [mature.preserved(row) for row in rows]
    bundle = ROOT / focused['bundle']
    summary = load(bundle / 'result.json')
    require(summary['schema'] == 'sirius-minimum-step-floor-focused-controls-v1'
            and summary['source_revision'] == CONFIG['current_product_source']
            and summary['result']['passed'] is True
            and summary['result']['completed_case_count'] == summary['result']['required_case_count'] == 6
            and summary['qualification_claimed'] is False
            and summary['admitted_external_domains'] == 0 and summary['pending_external_domains'] == 8,
            'focused archive is not the original complete finite result')
    independent = load(bundle / 'files/out/minimum-step-floor-current/independent-terminal-review.json')
    require(independent['passed'] is True and independent['source_revision'] == CONFIG['current_product_source'],
            'actual independent terminal review missing or different')
    for row in CONFIG['prior_evidence']:
        old = load(ROOT / row[2] / 'result.json')
        require(old['configure_accepted'] is False and old['compilation_started'] is False
                and old['numerical_execution_started'] is False and old['qualification_claimed'] is False,
                'prior failed configure was relabeled as accepted')
    return result


def terminal_records():
    # Keep the mature PID/birth recheck and explicit same-boot visibility scope.
    records = mature.retired_owners()
    require({row['path'] for row in records} == set(CONFIG['expected_owners']),
            'the exact ten reviewed owners changed')
    for row in records:
        path = ROOT / row['path']
        owner = load(path)
        expected = CONFIG['expected_owners'][row['path']]
        require(owner['passed'] is expected['passed'] and owner['child_exit'] == expected['child_exit']
                and owner['source'] == {'revision': expected['source_revision'], 'status': ''}
                and owner['stop_reason'] is None, 'original passing/refused owner disposition changed')
        before, after = path.parent / 'inputs-before.json', path.parent / 'inputs-after.json'
        require(load(before) == load(after), 'owner whole-input records differ')
        row['inputs_before'] = identity(before)
        row['inputs_after'] = identity(after)
    return records


def declared_inputs():
    """Fresh bounded scan of actual out/ declaration paths, not prose history.

    This is the existing frozen-input/argv/pin-path check expanded to discover a
    newly prepared screen. It is not a general Python dependency proof.
    """
    files, visited = [], 0
    for base, directories, names in os.walk(ROOT / 'out', followlinks=False):
        base = Path(base)
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
                    if name in {'pins', 'files'} and isinstance(child, dict):
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
        walk(value)
        for token in tokens:
            if not isinstance(token, str) or not token:
                continue
            token = token.split('=', 1)[-1]
            candidate = Path(token)
            candidates = [candidate] if candidate.is_absolute() else [ROOT / candidate, path.parent / candidate]
            require(not any(inside(item) or inside(item.resolve()) for item in candidates),
                    'remaining packet consumes retired scratch: ' + str(path) + ': ' + token)
        records[str(path.relative_to(ROOT))] = identity(path)
    return {'declaration_files': records, 'metadata_files_examined': len(records),
        'live_cleanup_root_path_found': False,
        'scope': 'Fresh selected out/ metadata frozen-input, pin/path and argv declarations; excludes the six retired packets and this cleanup apparatus. Not a general undeclared program-dependency proof.',
        'next_screen_rule': 'Use canonical/bin inputs and immutable focused/43ed/b220 archive paths. Retain the live finite packet until that screen is accepted if it still needs these paths.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-source', required=True)
    parser.add_argument('--focused-commit', required=True)
    parser.add_argument('--focused-manifest-sha256', required=True)
    parser.add_argument('--remove', action='store_true')
    parser.add_argument('--accepted-audit-sha256')
    args = parser.parse_args()
    for value, count in [(args.expected_source, 40), (args.focused_commit, 40), (args.focused_manifest_sha256, 64)]:
        require(re.fullmatch('[0-9a-f]{' + str(count) + '}', value) is not None, 'exact root-reviewed pins required')
    for name, expected in load(HERE / 'source-freeze.json')['pins'].items():
        require(identity(ROOT / name) == expected, 'cleanup preparation source changed')
    mature.ROOT, mature.SOURCE, mature.TARGETS = ROOT, args.expected_source, TARGETS
    mature.prefix_target = lambda path, targets=TARGETS: any(path == item or path.is_relative_to(item) for item in targets)
    source = source_join(args.expected_source)
    previous = None
    if args.remove:
        previous = load(HERE / 'audit.json')
        require(identity(HERE / 'audit.json')['sha256'] == args.accepted_audit_sha256
                and previous['status'] == 'verified-not-removed' and previous['source'] == source
                and previous['focused_commit'] == args.focused_commit
                and previous['focused_manifest_sha256'] == args.focused_manifest_sha256,
                'the exact independently reviewed audit is required')
        require(not RECEIPT.exists(), 'fresh fixed cleanup receipt root required')
    else:
        require(args.accepted_audit_sha256 is None and not (HERE / 'audit.json').exists(), 'fresh audit required')
    proof = preservation(args.focused_commit, args.focused_manifest_sha256)
    terminal = terminal_records()
    dependencies = declared_inputs()
    roots = snapshots()
    census = mature.consumer_census()
    document = {'schema': 'sirius-minimum-step-fixed-terminal-retirement-v1', 'status': 'verified-not-removed',
        'source': source, 'focused_commit': args.focused_commit, 'focused_manifest_sha256': args.focused_manifest_sha256,
        'preservation': proof, 'roots': roots, 'terminal_owners': terminal, 'declared_inputs': dependencies,
        'initial_census': census, 'removed': [], 'qualification_claimed': False,
        'global_consumer_absence_claimed': False, 'all_bin_build_roots_retained': True}
    if not args.remove:
        mature.write(HERE / 'audit.json', document)
        print(json.dumps({'audit': identity(HERE / 'audit.json'), 'fixed_roots': CONFIG['targets'],
            'terminal_owners': len(terminal), 'permission_gaps': len(census['permission_gaps'])}))
        return
    require(previous['preservation'] == proof and previous['roots'] == roots
            and previous['declared_inputs'] == dependencies, 'preserved roots/declarations changed after audit review')
    source_join(args.expected_source)
    document['final_census'] = mature.consumer_census()
    require(snapshots() == roots, 'owned bytes/modes/layout changed immediately before fixed removal')
    RECEIPT.mkdir(parents=True)
    for name in ['retire.py', 'mature_cleanup.py', 'configuration.json', 'source-freeze.json', 'commands.json', 'audit.json']:
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
        'Git_commit_tag_push_performed': False, 'all_bin_build_roots_retained': True}))


if __name__ == '__main__':
    main()
