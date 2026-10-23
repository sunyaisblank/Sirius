"""Passive one-attempt disposition. Exit zero means audit written, not frame accepted.

No subprocess, provider call, restart or PPM assertion is performed. The original
GTest alone checks rc==0, file existence and >1024 bytes, then removes its PPM.
"""
from pathlib import Path
import hashlib
import json
import math
import os
import re
import shlex
import xml.etree.ElementTree as ET

R = Path(__file__).resolve().parents[3]
W = Path(__file__).resolve().parent
P = W.parent
F = W / 'frame'
CASE = 'RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent'
LIMITS = {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
SEEN = {}


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size,
                'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


def doc(path):
    SEEN[str(path.relative_to(R)) if path.is_relative_to(R) else str(path)] = pin(path)
    return json.loads(path.read_text())


def need(value, message):
    if not value:
        raise ValueError(message)


def path_of(value):
    p = Path(value)
    return p if p.is_absolute() else R / p


def key(path):
    return str(path.relative_to(R)) if path.is_relative_to(R) else str(path)


def sealed(path, records):
    row = records[key(path)]
    actual = pin(path)
    need(actual == {k: row[k] for k in ('bytes', 'sha256')}, 'changed/unsealed input: ' + str(path))
    SEEN[key(path)] = actual
    return actual


def nonnegative_int(value):
    return type(value) is int and value >= 0


def finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def records(text, prefix):
    return [json.loads(line.split(prefix, 1)[1]) for line in text.splitlines() if prefix in line]


def check_owner(owner, before, after, source):
    need(owner['status'] == 'terminal' and owner['source'] == {'revision': source, 'status': ''}, 'owner source/terminal')
    need(owner['source_unchanged'] is True and owner['whole_inputs_unchanged'] is True and before == after,
         'source/frozen inputs changed')
    need(owner['cleanup_errors'] == [] and owner['remaining_owned_processes'] == [], 'owned cleanup incomplete')
    births = owner['observed_births']
    need(births and all(b['absent'] is True and b['current_start_ticks'] != b['saved_start_ticks'] for b in births),
         'saved owned births not retired')
    need(any(b['pid'] == owner['child_pid'] and b['saved_start_ticks'] == owner['child_start_ticks'] for b in births),
         'initial child birth missing')


def main():
    output = W / 'result.json'
    need(not output.exists(), 'refuse to overwrite prior disposition')
    result = {'schema': 'sirius-original-cold-frame-disposition-v1', 'attempted': (F / 'owner.json').is_file(),
              'accepted': False, 'audit_completed': False, 'failures': [],
              'scope': 'One original software Vulkan 128x128 one-SPP Kerr frame; no full science/Mandatory/release claim.',
              'limits_meaning': '4GiB sampled controller+owned-session RSS and 2GiB host reserve; not Vulkan residency or a GPU budget.',
              'output_assertions': 'Original completed zero-FESD GTest checks rc0/existence/>1024 bytes and deletes PPM; no new image/golden check.'}
    try:
        actions = doc(W / 'actions.json')
        action = actions['frame']
        owner = doc(F / 'owner.json')
        before, after = doc(F / 'inputs-before.json'), doc(F / 'inputs-after.json')
        source = action['revision']
        result.update(source=source, stop_reason=owner['stop_reason'], child_exit=owner['child_exit'],
                      wall_seconds=owner['wall_seconds'], peak_sampled_rss_kib=owner['peak_sampled_rss_kib'],
                      minimum_sampled_available_kib=owner['minimum_sampled_available_kib'])
        check_owner(owner, before, after, source)
        need(action['limits'] == owner['limits'] == LIMITS and action['kind'] == 'test'
             and type(action['expected_exit']) is int and action['expected_exit'] == 0
             and action['software_provider'] is True, 'original action/guard policy')
        need(owner['case'] == 'frame' and owner['argv'] == action['argv'], 'actual action argv')
        need(all(k in before for k in action['frozen_inputs']), 'declared input absent from frozen map')
        sealed(Path(__file__).resolve(), before)
        provider_path = P / 'providers.json'
        sealed(provider_path, before)
        provider = doc(provider_path)['providers']['candidate']
        library, icd = path_of(provider['library']), path_of(provider['ICD'])
        need(action['runtime_provider'] == str(library) and action['runtime_icd'] == str(icd), 'selected provider action')
        need(sealed(library, before) == provider['library_identity']
             and sealed(icd, before) == provider['ICD_identity'], 'actual dynamic provider identity')
        descriptor = doc(icd)
        described = Path(descriptor['ICD']['library_path'])
        need((described if described.is_absolute() else icd.parent / described).resolve() == library.resolve(), 'ICD library join')
        env = owner['declared_output_environment']
        need(env == {'GTEST_OUTPUT': 'xml:' + str(F / 'gtest.xml'), 'VK_DRIVER_FILES': str(icd),
                     'SIRIUS_VULKAN_DEVICE': '0', 'MESA_SHADER_CACHE_DISABLE': 'true'}, 'original declared runtime/cache environment')
        need(action.get('disable_software_shader_cache', True) is True, 'cold cache-disable policy')
        temp = path_of(action['build_environment']['TMPDIR'])
        need(temp == R / 'renders/lavapipe-owner-bulk-zero/cold-frame-tmp', 'workspace temporary output')
        need(set(action['build_environment']) == {'TMPDIR'}, 'extra environment override')
        births = {(b['pid'], b['saved_start_ticks']) for b in owner['observed_births']}
        mapped = []
        for m in owner['provider_mappings']:
            need((m['pid'], m['start_ticks']) in births and m['owned_session'] == owner['child_pid'], 'provider birth/session')
            for line in m['maps']:
                fields = line.split(None, 5)
                if len(fields) == 6 and fields[5] == str(library):
                    dev = tuple(int(x, 16) for x in fields[3].split(':'))
                    stat = library.stat()
                    need(int(fields[4]) == stat.st_ino and dev == (os.major(stat.st_dev), os.minor(stat.st_dev)), 'mapped provider inode/device')
                    mapped.append(m)
        need(mapped, 'actual same-session selected provider mapping missing')
        result['provider'] = {'library': str(library), 'identity': provider['library_identity'],
                              'ICD': str(icd), 'mapped_owned_births': sorted({(m['pid'], m['start_ticks']) for m in mapped}),
                              'environment_scope': 'Guard-declared inherited launch environment; no invented live-child environment readback.'}

        # The same current cohort authorizes this one frame, without copying historical provider identity.
        accepted_path = P / 'readers/controls-1-8-acceptance.json'
        sealed(accepted_path, before)
        accepted = doc(accepted_path)
        need(accepted['all8_original_controls_pass'] is True and accepted['all22_original_benefit_checks_pass'] is True,
             'original8/22 prerequisites')
        for name, field in [('control-8-acceptance.json', 'control8_acceptance'),
                            ('benefit.json', 'benefit'), ('decimal-agreement.json', 'decimal')]:
            p = P / 'readers' / name
            need(sealed(p, before) == accepted[field], 'cohort acceptance pin: ' + name)
        decimal = doc(P / 'readers/decimal-agreement.json')
        need(decimal['agrees_on_all22'] is True and decimal['all22_pass'] is True, 'Decimal22 prerequisite')

        recipe_path, build_path = W / 'link-recipe.json', W / 'build-result.json'
        sealed(recipe_path, before); sealed(build_path, before)
        recipe, build = doc(recipe_path), doc(build_path)
        link_owner_path = W / 'build/owner.json'
        link_before_path, link_after_path = W / 'build/inputs-before.json', W / 'build/inputs-after.json'
        for p in (link_owner_path, link_before_path, link_after_path):
            sealed(p, before)
        link_owner, link_before, link_after = doc(link_owner_path), doc(link_before_path), doc(link_after_path)
        check_owner(link_owner, link_before, link_after, source)
        need(link_owner['passed'] is True and link_owner['child_exit'] == 0 and link_owner['stop_reason'] is None, 'link not completed')
        sealed(recipe_path, link_before)
        need(recipe['source'] == build['source'] == source and build['passed'] is True
             and recipe['new_CPP_compilations'] == build['new_CPP_compilations'] == 0
             and recipe['changed_link_archive'] == 'sirius_backend' and build['original_objects_reused'] is True
             and build['original_render_binary_unchanged'] is True, 'actual original-object link authority')
        exe = path_of(build['executable'])
        need(sealed(exe, before) == build['executable_identity'], 'linked render executable')
        args, cwd = recipe['argv'], Path(recipe['cwd'])
        need(cwd == R / 'bin/linux-gcc' and args[0] == '/usr/bin/g++-14'
             and args[args.index('-o') + 1] == str(exe), 'original final link recipe/output')
        object_path = cwd / 'tests/render/CMakeFiles/sirius_render_tests.dir/app_render_execution_test.cpp.o'
        need(str(object_path.relative_to(cwd)) in args or str(object_path) in args, 'original testcase object not linked')
        sealed(object_path, link_before); sealed(object_path, before)
        backend = R / 'bin/linux-gcc/native-cooperative-transport-intergroup/matched/libsirius_backend.a'
        need(args.count(str(backend)) == 1 and sealed(backend, link_before) == recipe['backend'], 'mixed archive actual link join')
        need(recipe['resources'] == build['original_resources_copied_exactly'], 'original resource copy inventory')
        volume = exe.parent / 'resources'
        actual_names = {str(p.relative_to(volume)) for p in volume.rglob('*') if p.is_file()}
        need(actual_names == set(recipe['resources']), 'staged resource set')
        for name, identity in recipe['resources'].items():
            need(sealed(volume / name, before) == identity, 'staged original resource bytes: ' + name)
        compdb = R / 'bin/linux-gcc/compile_commands.json'
        sealed(compdb, before)
        rows = [r for r in doc(compdb) if Path(r['file']) == R / 'tests/render/app_render_execution_test.cpp']
        need(len(rows) == 1, 'selected testcase compile command cardinality')
        row = rows[0]
        tokens = row.get('arguments') or shlex.split(row['command'])
        compiled_output = Path(row['directory']) / row['output']
        need(compiled_output == object_path and '-c' in tokens and tokens[tokens.index('-c') + 1] == row['file'], 'compile object/source join')
        required = {'-DSIRIUS_CONTRACT_MODE=2', '-DSIRIUS_HAS_VULKAN_BACKEND=1', '-DSIRIUS_HAS_RETAINED_COMPUTE=1',
                    '-DSIRIUS_TEST_REQUIRE_VULKAN_RUNTIME=1', '-DSIRIUS_TEST_HAS_KERNELS=1', '-Werror'}
        includes = [tokens[i + 1] for i, t in enumerate(tokens[:-1]) if t == '-include']
        includes += [t[len('-include'):] for t in tokens if t.startswith('-include') and t != '-include']
        need(required.issubset(tokens) and str(R / 'tests/support/fail_on_skip.h') in includes,
             'required-runtime strict flags/skip refusal')
        for p in (R / row['file'], R / 'tests/support/fail_on_skip.h',
                  R / 'src/sirius/render/vulkan_renderer.cpp', R / 'src/sirius/render/render_evidence.cpp',
                  R / 'src/sirius/render/dispatch_governor.h', R / 'src/sirius/render/memory_governor.h',
                  R / 'src/sirius/render/session/render_session.cpp', R / 'src/sirius/app/cli/render_command.cpp'):
            sealed(p, before)
        runbook = R / 'tests/operational/vulkan_runtime_profile.cmake'
        sealed(runbook, before)
        need(action['argv'] == ['/usr/bin/cmake', '-DSIRIUS_RENDER_TESTS=' + str(exe), '-P', str(runbook)], 'original runbook argv')
        result['host_authority'] = {'executable': key(exe), 'identity': build['executable_identity'],
                                    'test_object': key(object_path), 'compile_command': row,
                                    'link_recipe': pin(recipe_path), 'link_owner': pin(link_owner_path),
                                    'scope': 'Current compile-command/source composition plus actual frozen reused object and link recipe; no new C++ compile or exhaustive historical tool-beforepin claim.'}
        result['audit_completed'] = True

        stdout, stderr = F / 'stdout.log', F / 'stderr.log'
        SEEN[key(stdout)], SEEN[key(stderr)] = pin(stdout), pin(stderr)
        text = stdout.read_text(errors='replace')
        errors = stderr.read_text(errors='replace')
        progress = re.findall(r'(\d+)/(16384)\b', text)
        snapshots = re.findall(r'\[Vulkan\] Retained progress: ([^\r\n]+)', errors)
        result['last_progress'] = list(map(int, progress[-1])) if progress else None
        result['last_complete_snapshot'] = snapshots[-1] if snapshots else None
        scene_authority = R / 'out/native-original-input-ownership/cold-frame/result.json'
        sealed(scene_authority, before)
        original_scene = doc(scene_authority)['scene']
        scene_rows = records(text, '[Session] Scene evidence: ') + records(text, '[Vulkan] Source scene evidence: ')
        need(scene_rows and all(s == original_scene for s in scene_rows), 'original08 complete scene identity')
        result['scene'] = scene_rows[-1]
        evidence_rows = records(text, '[Vulkan] Render evidence: ')
        result['XML_present'] = (F / 'gtest.xml').is_file()
        result['render_complete_marker'] = 'Render complete: ' in text
        result['completed_render_evidence_present'] = bool(evidence_rows)
        # A timeout is a valid negative disposition, never inferred to complete from progress.
        need(owner['passed'] is True and owner['stop_reason'] is None and owner['child_exit'] == 0
             and finite(owner['wall_seconds']) and 0 < owner['wall_seconds'] <= 300
             and 0 < owner['peak_sampled_rss_kib'] <= LIMITS['rss_kib']
             and owner['minimum_sampled_available_kib'] >= LIMITS['reserve_kib'], 'owned original frame did not pass its limits')
        xml_path = F / 'gtest.xml'
        SEEN[key(xml_path)] = pin(xml_path)
        xml = ET.parse(xml_path).getroot()
        cases = list(xml.iter('testcase'))
        need(len(cases) == 1 and cases[0].get('classname') + '.' + cases[0].get('name') == CASE
             and cases[0].get('status') == 'run' and cases[0].get('result') == 'completed'
             and not any(list(xml.iter(tag)) for tag in ('failure', 'error', 'skipped'))
             and all(int(node.get(k, '0')) == 0 for node in [xml, *xml.iter('testsuite')]
                     for k in ('failures', 'errors', 'skipped', 'disabled'))
             and int(xml.get('tests', '0')) == 1, 'exact completed one-case zero-FESD XML')
        need(result['render_complete_marker'] and '[Session] Vulkan render complete: ' in text
             and len(evidence_rows) == 1, 'missing/ambiguous completed renderer evidence')
        e = evidence_rows[0]
        result['render_evidence'] = e
        need(e['schema'] == 'sirius-vulkan-render-v1' and e['route'] == 'retained' and e['source_owner'] == 'host'
             and e['metric'] == 'Kerr' and e['width'] == e['height'] == 128 and e['precision'] == 'fp32'
             and e['device_index'] == 0 and isinstance(e['device_name'], str) and e['device_name'], 'original actual retained renderer route')
        need(all(nonnegative_int(e[k]) for k in ('budget_bytes', 'usable_bytes', 'allocated_bytes', 'work_items',
                                                'work_tile_edge', 'ray_capacity', 'maximum_dispatch_rays', 'dispatches',
                                                'camera_batches', 'accepted_intervals', 'target_overshoots',
                                                'batch_subdivisions', 'safety_fallbacks')), 'counter types')
        need(0 < e['allocated_bytes'] <= e['usable_bytes'] == int(float(e['budget_bytes']) * 0.5) < e['budget_bytes'], 'original residency half-budget')
        need(e['work_items'] == 16384 and e['work_tile_edge'] == 1 and e['ray_capacity'] == 64
             and 0 < e['maximum_dispatch_rays'] <= 64 and e['retained_timing']['projection_capacity'] == 128
             and e['camera_batches'] > 0 and e['accepted_intervals'] > 0, 'original full-frame work/capacity/counters')
        stages = e['retained_stage_dispatches']
        need(type(stages) is list and len(stages) == 7 and all(nonnegative_int(n) for n in stages)
             and sum(stages) == e['dispatches'] > 0, 'seven stage command counts/total')
        stage_rows = e['retained_stage_rows']['stages']
        need(len(stage_rows) == 7 and [r['submissions'] for r in stage_rows] == stages, 'producer stage row/count join')
        need(all(finite(e[k]) and e[k] >= 0 for k in ('wall_seconds', 'dispatch_seconds', 'maximum_dispatch_ms'))
             and 0 < e['wall_seconds'] <= owner['wall_seconds'] and e['maximum_dispatch_ms'] <= 1000, 'original governed physical work stop')
        header = re.findall(r'\[Vulkan\] Retained renderer: [^\r\n]*; budget (\d+) MiB, (\d+) ray rows, (\d+) projection slots, (\d+) host work tiles of (\d+)px', text)
        need(len(header) == 1 and tuple(map(int, header[0])) == (e['budget_bytes'] // (1024 * 1024), 64, 128, 16384, 1), 'actual planner header/evidence join')
        result['accepted'] = True
    except (AssertionError, ValueError, KeyError, TypeError, OSError, ET.ParseError) as error:
        result['failures'].append(type(error).__name__ + ': ' + str(error))
    result['inputs_read'] = SEEN
    result['accepted_requires'] = 'Owned clean completion, original8/22+Decimal, actual provider/host/resource/source joins, original scene, exact positive XML and completed governed render evidence.'
    result['limits'] = LIMITS
    result['no_raw_framebuffer_or_full_science_claim'] = True
    # Exclusive publication; readers never overwrite a prior result or touch the owned attempt.
    with output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({'accepted': result['accepted'], 'audit_completed': result['audit_completed'],
                      'stop_reason': result.get('stop_reason'), 'output': str(output)}))


if __name__ == '__main__':
    main()
