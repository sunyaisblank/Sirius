"""Root-operated exact preservation of the completed original issue101 prerequisites.

This helper never builds, dispatches, edits a live prerequisite, or changes Git.
The final manifest is written only after acceptance and byte/birth rechecks.
"""
from pathlib import Path
import ast
import hashlib
import json
import re
import shutil
import stat
import struct
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
ROOT = Path.cwd().resolve()
WORK = ROOT / 'out/computed-carry-exit-original-prerequisites'
BUILD = ROOT / 'bin/linux-gcc/computed-carry-exit-original-prerequisites'
ROOTS = (WORK, BUILD)
ARCHIVE = ROOT / 'attestations/software-vulkan/4bcd987/computed-carry-exit-original-prerequisites'
PRODUCT = '4bcd9875cd9308deabed06a5fa2f17d2321b0cbc'


def identity(path):
    path = Path(path)
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': path.stat().st_size, 'sha256': digest}


def load(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def relative(path):
    return str(Path(path).relative_to(ROOT))


def no_symlink_ancestors(path):
    path = Path(path)
    assert path.is_relative_to(ROOT)
    for candidate in [path, *path.parents]:
        if candidate == ROOT:
            break
        assert not candidate.is_symlink(), candidate


def inventory():
    files, directories = {}, {}
    for owned in ROOTS:
        no_symlink_ancestors(owned)
        assert owned.is_dir(), owned
        for path in [owned, *sorted(owned.rglob('*'))]:
            info = path.lstat()
            assert not stat.S_ISLNK(info.st_mode), path
            assert stat.S_ISDIR(info.st_mode) or stat.S_ISREG(info.st_mode), path
            data = {'device': info.st_dev, 'inode': info.st_ino,
                    'mode': stat.S_IMODE(info.st_mode), 'mtime_ns': info.st_mtime_ns,
                    'ctime_ns': info.st_ctime_ns}
            name = relative(path)
            if stat.S_ISDIR(info.st_mode):
                directories[name] = data
            else:
                files[name] = {**identity(path), **data}
    return files, directories


def current_birth(pid):
    try:
        return Path('/proc/%d/stat' % pid).read_text().rsplit(')', 1)[1].split()[19]
    except (FileNotFoundError, ProcessLookupError):
        return None


def fresh_birth_checks(receipts):
    boot = Path('/proc/sys/kernel/random/boot_id').read_text().strip()
    checks = []
    for path in receipts:
        owner = load(path)
        assert owner['boot_id'] == boot, path
        births = [(b['pid'], b['saved_start_ticks']) for b in owner['observed_births']]
        births.append((owner['controller']['pid'], owner['controller']['start_ticks']))
        births.append((owner['child_pid'], owner['child_start_ticks']))
        for pid, saved in sorted(set(births)):
            now = current_birth(pid)
            assert now != saved, ('active observed birth', path, pid, saved)
            checks.append({'owner': relative(path), 'pid': pid, 'saved_start_ticks': saved,
                           'current_start_ticks': now, 'absent': True})
    preparation = load(WORK / 'preparation-owner.json')
    assert preparation['status'] == 'terminal' and preparation['boot_id'] == boot
    assert current_birth(preparation['pid']) != preparation['start_ticks']
    return checks


def acceptance_gate():
    R = ROOT
    old_work = 'out/computed-carry-exit-original-prerequisites/'
    old_stage = 'bin/linux-gcc/computed-carry-exit-original-prerequisites/controls/'

    def preserved(name):
        path = ROOT / name
        assert any(path.is_relative_to(owned) for owned in ROOTS), name
        no_symlink_ancestors(path)
        assert path.is_file() and stat.S_ISREG(path.lstat().st_mode), name
        return path

    producer = load(WORK / 'producer-acceptance.json')
    assert producer['passed'] and producer['unchanged_array_count'] == 27 and len(producer['changed_arrays']) == 13
    assert set(producer['actual_ELF_owner_arrays']) == {'baseline', 'candidate'}
    # The complete finite archive gate from the reviewed benefit preparation is inserted here.
    source = json.loads(preserved(old_work+'operator-context.json').read_text())
    def prerequisite_owner(label, provider=False):
        owner = json.loads(preserved(old_work+label+'/owner.json').read_text())
        assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0, label
        assert owner['source'] == {'revision':source['revision'], 'status':''}, label
        assert owner['source_unchanged'] and owner['whole_inputs_unchanged'], label
        assert owner['stop_reason'] is None and owner['cleanup_errors'] == [] and owner['remaining_owned_processes'] == [], label
        assert owner['observed_births'] and all(b['absent'] for b in owner['observed_births']), label
        if provider:
            assert owner['provider_mappings'], label
        assert preserved(old_work+label+'/inputs-before.json').read_bytes() == preserved(old_work+label+'/inputs-after.json').read_bytes(), label
        return owner

    controls = json.loads(preserved(old_work+'control-acceptance.json').read_text())
    routes = json.loads(preserved(old_work+'endpoint-route-acceptance.json').read_text())
    assert controls['passed'] and routes['passed'] and len(controls['original_controls']) == 8
    schema = controls['timestamp_schema']
    assert (schema['complete_readbacks'],schema['ordered_markers'],schema['actual_module_command_joins']) == (32,51,54)
    assert schema['all_readbacks_exact_to_repaired_original_baseline']
    assert routes['complete_readbacks'] == 56 and routes['complete_readbacks_exact_to_prior_original_control']
    assert len(routes['routes']) == 4 and {(v['declined_RTE32'],v['wide']) for v in routes['routes']} == {(d,w) for d in [False,True] for w in [False,True]}
    original_controls = json.loads((R/'attestations/software-vulkan/41fe5f2/portable-denominator-preparation/files/out/portable-denominator-preparation/original-controls/actions.json').read_text())
    archived_actions = json.loads(preserved(old_work+'actions.json').read_text())
    for index in range(1,9):
        label = 'control-'+str(index)
        owner = prerequisite_owner(label, True)
        assert owner['argv'] == archived_actions[label]['argv']
        assert owner['argv'][1] == original_controls[label]['argv'][1]
        assert owner['limits'] == original_controls[label]['limits']
        row = controls['original_controls'][index-1]
        assert row['label'] == label and row['owner'] == identity(preserved(old_work+label+'/owner.json'))
        assert row['xml'] == identity(preserved(old_work+label+'/gtest.xml'))
        xml = ET.parse(preserved(old_work+label+'/gtest.xml')).getroot()
        assert int(xml.get('tests','-1')) == 1
        assert not any(xml.findall('.//'+tag) for tag in ['failure','error','skipped'])
        for node in [xml,*xml.findall('.//testsuite')]:
            assert all(int(node.get(key,'0')) == 0 for key in ['failures','errors','skipped','disabled'])
        cases = xml.findall('.//testcase')
        assert len(cases) == 1 and cases[0].get('status') == 'run' and cases[0].get('result') == 'completed'
        assert owner['argv'][1] == '--gtest_filter='+cases[0].get('classname')+'.'+cases[0].get('name')
    protected = json.loads(preserved(old_work+'protected-before.json').read_text())
    assert len(protected) == controls['protected_payloads_unchanged'] == 114
    for name,pin in protected.items():
        assert identity(R/name) == pin, name
    products = json.loads(preserved(old_work+'products-acceptance.json').read_text())
    assert products['passed'] and len(products['modules']) == 35 and sum(row['changed'] for row in products['modules']) == 15
    for label in ['compile-products-stock','compile-products-candidate','verify-products','native-policy-selfcheck',
                  'compile-compute-sealed','compile-observer','link-candidate','arithmetic/validate-arithmetic']:
        prerequisite_owner(label)
    assert json.loads(preserved(old_work+'arithmetic/arithmetic-result.json').read_text())['passed']
    for mode in ['stock','candidate']:
        label = 'arithmetic/execute-gpu-'+mode
        prerequisite_owner(label, True)
        assert preserved(old_work+label+'/actual.bin').read_bytes() == preserved(old_work+'arithmetic/expected.bin').read_bytes()
    # Reuse the sealed original ELF reader's definitions, without its executable footer.
    reader_tree = ast.parse(preserved(old_work+'verify_producer.py').read_text())
    reader_defs = ast.Module(body=[n for n in reader_tree.body if isinstance(n,(ast.Import,ast.ImportFrom,ast.FunctionDef))],type_ignores=[])
    reader = {}
    exec(compile(reader_defs,'sealed_original_ELF_reader','exec'),reader)
    for mode in ['baseline','candidate']:
        directory = 'stock' if mode == 'baseline' else 'candidate'
        generated = 'bin/linux-gcc/computed-carry-exit-original-prerequisites/'+directory+'/'
        executable = preserved(old_stage+'original_'+mode+'_controls')
        assert identity(executable) == producer[mode+'_executable']
        actual,_ = reader['elf_arrays'](executable,reader['arrays'](preserved(generated+'retained_kernels.h')))
        assert actual == producer['actual_ELF_arrays'][mode]
        owner_arrays = {}
        for count,name,body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',preserved(generated+'retained_owner_kernels.h').read_text(),re.S):
            words = [int(v) for v in re.findall(r'(\d+)u',body)]
            assert len(words) == int(count)
            owner_arrays[name] = struct.pack('<'+str(len(words))+'I',*words)
        assert set(owner_arrays) == {'kTransportShader','kEndpointShader'}
        actual,_ = reader['elf_arrays'](executable,owner_arrays,'retained_owner')
        assert actual == producer['actual_ELF_owner_arrays'][mode]
    for label in ['compile-compute-sealed','compile-observer']:
        proof = producer['compile_provenance'][label]
        assert proof['owner'] == identity(preserved(old_work+label+'/owner.json'))
        assert proof['actual_dependencies'][old_stage+'headers/retained_kernels.h'] == identity(preserved(old_stage+'headers/retained_kernels.h'))
    for name in ['retained_owner_kernels.h','retained_owner_layout.h']:
        assert producer['compile_provenance']['compile-compute-sealed']['actual_dependencies'][old_stage+'headers/'+name] == identity(preserved(old_stage+'headers/'+name))
    assert producer['link_owner'] == identity(preserved(old_work+'link-candidate/owner.json'))
    assert producer['link_map'] == identity(preserved(old_stage+'candidate.map')) and producer['old_archive_compute_member_excluded']
    return {'controls': controls, 'routes': routes, 'products': products, 'producer': producer}


def receipts_and_inputs(files):
    context = load(WORK / 'operator-context.json')
    assert identity(WORK / 'source-freeze.json') == context['source_freeze']
    for name, pin in load(WORK / 'source-freeze.json').items():
        assert identity(name) == pin, name
    declarations = {}
    for folder in (WORK, WORK / 'arithmetic'):
        for label, action in load(folder / 'actions.json').items():
            declarations[folder / label / 'owner.json'] = (folder / label, action)
    receipts = sorted(p for owned in ROOTS for p in owned.rglob('owner.json'))
    assert len(receipts) == 32 and set(receipts) == set(declarations), ('unexpected receipt set', len(receipts))
    by_content = {}
    for name, pin in files.items():
        by_content.setdefault((pin['bytes'], pin['sha256']), []).append(name)
    for key in ['initial_actions', 'arithmetic_actions']:
        pin = context[key]
        assert (pin['bytes'], pin['sha256']) in by_content, ('missing historic actions', key)
    recovery, joins, current_cache = {}, {}, {}
    for receipt in receipts:
        folder, action = declarations[receipt]
        owner = load(receipt)
        assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0, receipt
        assert owner['stop_reason'] is None and owner['cleanup_errors'] == [] and owner['remaining_owned_processes'] == [], receipt
        assert owner['observed_births'] and all(b['absent'] for b in owner['observed_births']), receipt
        assert owner['whole_inputs_unchanged'] and owner['source_unchanged'], receipt
        assert owner['source'] == {'revision': context['revision'], 'status': ''}, receipt
        assert action['revision'] == context['revision'] and owner['limits'] == action['limits'], receipt
        assert owner['argv'] == [str(s).replace('{folder}', str(folder)) for s in action['argv']], receipt
        if action.get('software_provider'):
            assert owner['provider_mappings'], receipt
        before = load(folder / 'inputs-before.json')
        assert before == load(folder / 'inputs-after.json'), receipt
        for name, pin in before.items():
            expected = {k: pin[k] for k in ['bytes', 'sha256']}
            path = ROOT / name
            if name not in current_cache:
                current_cache[name] = (identity(path) if path.is_file() else None, str(path.resolve()))
            current, resolved = current_cache[name]
            if current == expected and resolved == pin['resolved_path']:
                joins[name] = {**expected, 'resolved_path': pin['resolved_path']}
                continue
            key = (expected['bytes'], expected['sha256'])
            if key in recovery:
                recovery[key]['uses'].append({'owner': relative(receipt), 'input': name})
                continue
            matches = by_content.get(key, [])
            if matches:
                origin = ROOT / matches[0]
                assert identity(origin) == expected
                source = {'owned_source': relative(origin)}
            else:
                # A changed tracked input can be recovered only from its exact owner revision.
                assert path.is_relative_to(ROOT), ('unrecoverable external input', name, expected)
                logical = relative(path)
                result = subprocess.run(['git', 'show', owner['source']['revision'] + ':' + logical], capture_output=True)
                assert result.returncode == 0, ('unrecoverable historic input', name, expected)
                assert {'bytes': len(result.stdout), 'sha256': hashlib.sha256(result.stdout).hexdigest()} == expected, name
                source = {'git_revision': owner['source']['revision'], 'git_path': logical}
            recovery[key] = {'expected': expected, 'current': current, 'source': source,
                             'uses': [{'owner': relative(receipt), 'input': name}]}
    return receipts, recovery, joins


def main():
    assert sys.argv[1:] == ['--seal']
    assert ROOT == Path('/home/astra/.project/Sirius')
    assert not ARCHIVE.exists(), ('archive already exists; preserve any partial attempt', ARCHIVE)
    no_symlink_ancestors(ARCHIVE)
    assert subprocess.run(['git', 'check-ignore', '-q', str(ARCHIVE)]).returncode == 0, 'archive must remain ignored'
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    assert not subprocess.check_output(['git', 'status', '--porcelain=v1'])
    files, directories = inventory()
    accepted = acceptance_gate()  # Stops before any archive write if a control is incomplete.
    receipts, recovery, current_joins = receipts_and_inputs(files)
    birth_checks = fresh_birth_checks(receipts)
    assert inventory() == (files, directories), 'live roots changed during acceptance'
    ARCHIVE.mkdir(parents=True)
    archive_files, modes = {}, {}
    for name, pin in files.items():
        source = ROOT / name
        destination = ARCHIVE / 'files' / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        destination.chmod(pin['mode'])
        expected = {key: pin[key] for key in ['bytes', 'sha256']}
        assert identity(source) == expected and identity(destination) == expected, name
        assert stat.S_IMODE(destination.stat().st_mode) == pin['mode'], name
        archive_files['files/' + name] = expected
        modes['files/' + name] = pin['mode']
    recovered = []
    for item in recovery.values():
        pin, source = item['expected'], item['source']
        target = ARCHIVE / 'input-versions' / pin['sha256']
        target.parent.mkdir(exist_ok=True)
        if 'owned_source' in source:
            shutil.copyfile(ROOT / source['owned_source'], target)
        else:
            result = subprocess.run(['git', 'show', source['git_revision'] + ':' + source['git_path']], capture_output=True, check=True)
            target.write_bytes(result.stdout)
        assert identity(target) == pin
        target.chmod(0o444)
        name = str(target.relative_to(ARCHIVE))
        archive_files[name] = pin
        modes[name] = 0o444
        recovered.append({**item, 'preserved_version': name})
    write(ARCHIVE / 'input-recovery.json', {'changed_inputs': recovered, 'unchanged_current_joins': current_joins,
          'meaning': 'Historic consumed bytes are preserved exactly; every current mismatch has a recorded immutable recovery. No live input was repaired or altered.'})
    archive_files['input-recovery.json'] = identity(ARCHIVE / 'input-recovery.json')
    modes['input-recovery.json'] = stat.S_IMODE((ARCHIVE / 'input-recovery.json').stat().st_mode)
    disposition = {'issue': 101, 'source_revision': PRODUCT, 'operator_revision': load(WORK / 'operator-context.json')['revision'],
                   'decision': 'PROCEED_TO_ORIGINAL_BENEFIT_GATES', 'original_eight_controls_passed': True,
                   'four_endpoint_routes_passed': True, 'complete_timestamp_readbacks': 32,
                   'ordered_timestamp_markers': 51, 'actual_command_module_joins': 54,
                   'complete_route_readbacks': 56, 'protected_payloads_unchanged': 114,
                   'fresh_matched_compile_header_seal': True, 'production_adopted': False,
                   'qualification_complete': False,
                   'scope': 'Original finite compiler/arithmetic/eight-control prerequisites only; no benefit or scientific/native/release verdict.'}
    disposition_name = 'files/' + relative(WORK / 'disposition.json')
    assert not (WORK / 'disposition.json').exists(), 'disposition must not alter or replace a live prerequisite'
    write(ARCHIVE / disposition_name, disposition)
    archive_files[disposition_name] = identity(ARCHIVE / disposition_name)
    modes[disposition_name] = stat.S_IMODE((ARCHIVE / disposition_name).stat().st_mode)
    for name, pin in sorted(directories.items(), key=lambda entry: len(Path(entry[0]).parts), reverse=True):
        destination = ARCHIVE / 'files' / name
        destination.mkdir(parents=True, exist_ok=True)
        shutil.copystat(ROOT / name, destination)
        assert stat.S_IMODE(destination.stat().st_mode) == pin['mode']
    assert inventory() == (files, directories), 'live roots changed while copying'
    final_birth_checks = fresh_birth_checks(receipts)
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == head
    assert not subprocess.check_output(['git', 'status', '--porcelain=v1'])
    for name, pin in archive_files.items():
        assert identity(ARCHIVE / name) == pin, name
    assert {str(p.relative_to(ARCHIVE)) for p in ARCHIVE.rglob('*') if p.is_file()} == set(archive_files)
    manifest = {'schema': 'sirius-raw-diagnostic-bundle-v1', 'qualification_claimed': False,
                'product_revision': PRODUCT, 'operator_revision': load(WORK / 'operator-context.json')['revision'],
                'seal_revision': head, 'roots': [relative(p) for p in ROOTS], 'files': archive_files,
                'preserved_file_modes': modes, 'live_file_count': len(files),
                'live_directory_identities': directories, 'file_count': len(archive_files),
                'terminal_owner_receipts': [relative(p) for p in receipts],
                'accepted_reports': {name: identity(WORK / name) for name in [
                    'products-acceptance.json', 'arithmetic/arithmetic-result.json',
                    'producer-acceptance.json', 'control-acceptance.json', 'endpoint-route-acceptance.json']},
                'initial_birth_recheck': birth_checks, 'final_birth_recheck': final_birth_checks,
                'generated_files': ['input-recovery.json', disposition_name] + [v['preserved_version'] for v in recovered],
                'manifest_files_excluded': ['manifest.json', 'whole-byte-manifest.json'],
                'sealing_helper': {'path': relative(Path(__file__).resolve()), **identity(__file__)},
                'scope': 'Exact two-root raw snapshot, all32 terminal owners, complete original finite acceptances and historic input recovery. Useful source/results require separate Git preservation. No live source, output or Git mutation.'}
    write(ARCHIVE / 'whole-byte-manifest.json', manifest)
    write(ARCHIVE / 'manifest.json', manifest)  # Final completion marker consumed by benefit preparation.
    assert (ARCHIVE / 'manifest.json').read_bytes() == (ARCHIVE / 'whole-byte-manifest.json').read_bytes()
    print(json.dumps({'archive': relative(ARCHIVE), 'live_files': len(files), 'terminal_owners': len(receipts),
                      'recovered_historic_inputs': len(recovered), 'manifest': identity(ARCHIVE / 'manifest.json'),
                      'decision': disposition['decision'], 'live_inputs_and_Git_unchanged': True}))


if __name__ == '__main__':
    main()
