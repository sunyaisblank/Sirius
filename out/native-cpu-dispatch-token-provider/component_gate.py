"""Read sealed packed/raw and arithmetic evidence without restoring retired products."""
from pathlib import Path
import hashlib
import json
import re
import zipfile

R = Path.cwd()
MANIFEST = R / 'attestations/qualification/a47e219b/native-cpu-dispatch-admission-token-finite-compatibility/MANIFEST.json'
MANIFEST_SHA = '966a0d0013ab2363726d22654b1ac76dc58a23118972132481fe41c53d1221e2'
AS_GATE = R / 'out/native-cpu-inline-private-provider/component_gate.py'
AS_GATE_SHA = 'f8e57b5d3d78baad2891a44b460f29216172a92cdf6c422c1cc2d69ea559b045'
S = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8'
HEADERS = ('src/gallium/drivers/llvmpipe/lp_native_cpu_bridge.h',
           'src/gallium/drivers/llvmpipe/sirius_native_bridge.h')
STAGES = ('transport', 'endpoint', 'dopri_phase')

LAYOUT_ASSERTIONS = """/* AT passive-observer layout authority: compile-time only. */
#include <stddef.h>
_Static_assert(sizeof(struct lp_cs_job_info) == 112, "AT job size");
_Static_assert(offsetof(struct lp_cs_job_info, grid_size) == 0, "AT job grid");
_Static_assert(sizeof(((struct lp_cs_job_info *)0)->grid_size) == 12, "AT job grid width");
_Static_assert(sizeof(((struct lp_cs_job_info *)0)->grid_size[0]) == 4, "AT job grid element width");
_Static_assert(sizeof(float) == 4, "AT descriptor table element width");
_Static_assert(offsetof(struct lp_cs_job_info, native_preflight_passed) == 63, "AT positive job byte");
_Static_assert(sizeof(((struct lp_cs_job_info *)0)->native_preflight_passed) == 1, "AT positive job byte width");
_Static_assert(offsetof(struct lp_cs_job_info, native_cpu_stage) == 65, "AT job stage");
_Static_assert(sizeof(((struct lp_cs_job_info *)0)->native_cpu_stage) == 1, "AT job stage width");
_Static_assert(offsetof(struct lp_cs_job_info, current) == 72, "AT current views pointer");
_Static_assert(sizeof(((struct lp_cs_job_info *)0)->current) == 8, "AT current pointer width");
_Static_assert(offsetof(struct lp_cs_exec, jit_resources) + offsetof(struct lp_jit_resources, constants) + sizeof(struct lp_jit_buffer) == 24, "AT descriptor table view");
_Static_assert(sizeof(struct lp_jit_buffer) == 16, "AT buffer view size");
_Static_assert(offsetof(struct lp_jit_buffer, u) == 0, "AT buffer pointer");
_Static_assert(sizeof(((struct lp_jit_buffer *)0)->u) == 8, "AT buffer pointer width");
_Static_assert(offsetof(struct lp_jit_buffer, num_elements) == 8, "AT buffer extent");
_Static_assert(sizeof(((struct lp_jit_buffer *)0)->num_elements) == 4, "AT buffer extent width");
_Static_assert(sizeof(struct lp_descriptor) == 64, "AT descriptor stride");
_Static_assert(offsetof(struct lp_descriptor, buffer) == 0, "AT descriptor buffer origin");
/* End AT passive-observer layout authority. */
"""


def add_layout_assertions(headers):
    result=dict(headers);name=HEADERS[0]
    assert LAYOUT_ASSERTIONS not in result[name]
    result[name]=LAYOUT_ASSERTIONS + result[name]
    return result


def remove_layout_assertions(headers):
    result=dict(headers);name=HEADERS[0]
    assert result[name].startswith(LAYOUT_ASSERTIONS) and result[name].count(LAYOUT_ASSERTIONS)==1
    result[name]=result[name][len(LAYOUT_ASSERTIONS):]
    return result



def doc(path):
    return json.loads(path.read_text())


def data_pin(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def checked(path, row):
    assert pin(path) == {key: row[key] for key in ('bytes', 'sha256')}, path


def terminal(path):
    value = doc(path)
    assert value['status'] == 'terminal' and value['passed'] and value['child_exit'] == 0
    assert value['source_unchanged'] and value['whole_inputs_unchanged']
    assert not value['remaining_owned_processes'] and not value['cleanup_errors']
    return value


def contract_inputs(contract):
    for name, expected in contract['source_inputs'].items():
        checked(R / name, expected)


def bundle(path, sha):
    assert isinstance(sha, str) and re.fullmatch(r'[0-9a-f]{64}', sha), 'Dispatch-token component manifest is not yet bound'
    assert pin(path)['sha256'] == sha
    manifest = doc(path)
    assert manifest['finite_component_passed'] and not manifest['provider_adopted']
    archive = R / manifest['archive']['path']
    checked(archive, manifest['archive'])
    with zipfile.ZipFile(archive) as zipped:
        assert set(zipped.namelist()) == set(manifest['members'])
        files = {name: zipped.read(name) for name in zipped.namelist()}
    for name, expected in manifest['members'].items():
        assert data_pin(files[name]) == expected, name
    for row in manifest['owners']:
        checked_bytes = files[row['path']]
        assert data_pin(checked_bytes) == {key: row[key] for key in ('bytes', 'sha256')}
        value = json.loads(checked_bytes)
        assert value['status'] == 'terminal' and value['passed'] and value['child_exit'] == 0
        assert value['source'] == row['source'] and value['source_unchanged'] and value['whole_inputs_unchanged']
        assert not value['remaining_owned_processes'] and not value['cleanup_errors']
        assert value['observed_births'] and all(item['absent'] for item in value['observed_births'])
        assert row['all_owned_births_reaped']
        folder = row['path'].removesuffix('owner.json')
        assert files[folder + 'inputs-before.json'] == files[folder + 'inputs-after.json']
    return manifest, archive, files


def header_patch(texts, patch, reverse=False):
    """Apply only the two declared exact unified diffs, entirely in memory."""
    lines=patch.splitlines(keepends=True);i=0;changed={}
    while i<len(lines):
        assert lines[i].startswith('--- a/') and lines[i+1].startswith('+++ b/')
        name=lines[i][6:].rstrip('\n');assert lines[i+1][6:].rstrip('\n')==name and name in HEADERS
        original=texts[name].splitlines(keepends=True);out=[];cursor=0;i+=2
        while i<len(lines) and not lines[i].startswith('--- a/'):
            h=re.fullmatch(r'@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@[^\n]*\n',lines[i]);assert h,lines[i]
            start=int(h[3] if reverse else h[1])-1;assert cursor<=start
            out+=original[cursor:start];cursor=start;i+=1;before=[];after=[]
            while i<len(lines) and not lines[i].startswith(('@@ ','--- a/')):
                kind=lines[i][0];assert kind in ' +-'
                if kind in ' -':before.append(lines[i][1:])
                if kind in ' +':after.append(lines[i][1:])
                i+=1
            old,new=(after,before) if reverse else (before,after)
            assert original[cursor:cursor+len(old)]==old,(name,start)
            out+=new;cursor+=len(old)
        out+=original[cursor:];changed[name]=''.join(out)
    assert set(changed)==set(HEADERS)
    return changed


def component():
    manifest,archive,files=bundle(MANIFEST,MANIFEST_SHA)
    import importlib.util
    assert pin(AS_GATE)['sha256']==AS_GATE_SHA
    spec=importlib.util.spec_from_file_location('sealed_AS_components',AS_GATE)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    parent=module.component()
    prefix='out/native-cpu-dispatch-admission-token/'
    read=lambda name:json.loads(files[prefix+name])
    derivation,build,code,compatibility,linked=[read(n) for n in
        ('derivation.json','build-result.json','code-proof.json','compatibility.json','link-result.json')]
    assert len(manifest['members'])==72 and len(manifest['owners'])==3
    assert build['passed'] and code['passed'] and compatibility['passed'] and linked['passed']
    assert manifest['new_production_compilations']==build['new_production_compilations']==4
    assert manifest['new_diagnostic_compilations']==build['new_diagnostic_compilations']==1
    assert manifest['original_fixture_test_raw_math_and_admission_rebuilds']==build['original_fixture_test_raw_math_and_admission_rebuilds']==0
    assert manifest['fixture_invocations']==4 and manifest['fixture_modes']==['public','admitted']
    assert manifest['results']==compatibility['fixtures'] and manifest['code']==code
    assert manifest['source_derivation']['member']==prefix+'derivation.json'
    assert manifest['job_capability']==derivation['capability']
    assert manifest['original_public_wrapper_sources_preserved'] and manifest['original_pure_admission_object_reused']
    assert manifest['no_job_layout_change'] and manifest['no_new_state_buffers_allocations_or_caches']
    assert manifest['raw3_math_and_ISPC_guards_unchanged'] and not manifest['actual_Mesa_job_token_branch_observed']
    assert linked['fixture_GOT_targets_inverse_checked'] and linked['all_linked_functions_associated']
    assert linked['runtime_minimum_not_v3_or_v4']
    assert set(manifest['complete_link_association'])=={'transport','endpoint'}
    for stage,count in (('transport',18),('endpoint',17)):
        summary=manifest['complete_link_association'][stage]
        assert summary['objects']==count and summary['exact_preflight_public_and_private_target_joins']==15
        proof=read(stage+'-link-proof.json')
        assert proof['passed'] and len(proof['objects'])==count
        assert len(proof['exact_preflight_public_and_private_target_joins'])==15
        assert sum(len(row['functions']) for row in proof['objects'])==summary['functions']
    for row in manifest['reused_objects']:checked(R/row['path'],row)
    assert manifest['reused_objects']==derivation['reused_objects']
    assert derivation['new_compilations_if_selected']==4
    assert derivation['guard_inverse_restores_AP_header'] and derivation['default_false_projection_restores_original_strict_guard']
    assert derivation['Mesa_patch_inverse_restores_exact_two_headers']
    assert derivation['existing_public_wrapper_sources_preserved'] and derivation['original_pure_admission_object_reused']
    assert derivation['no_job_layout_change'] and derivation['no_new_state_buffers_allocations_or_caches']
    assert derivation['raw3_generated_sources_objects_and_all_ISPC_static_generic_guards_unchanged']
    assert compatibility['fixture_invocations']==4
    assert compatibility['original_assertions_counts_and_sources_unchanged_per_invocation']
    assert compatibility['actual_unchanged_pure_preflight_before_private_calls']
    assert compatibility['trusted_calls_after_failed_preflight']==0 and compatibility['public_and_private_routes_validated']
    assert set(compatibility['fixtures'])=={'public','admitted'}
    for mode in ('public','admitted'):
        assert set(compatibility['fixtures'][mode])=={'transport','endpoint'}
        for stage,row in compatibility['fixtures'][mode].items():
            result=row['result'];assert result['passed']
            if stage=='transport':
                for name,value in dict(case_count=38,active_rows=38,capacity_rows=40,
                        packets_per_implementation=5,raw_calls=10,compared_output_words=187960,
                        input_unchanged=True,inactive_sentinels_preserved=True,full_word_match=True,
                        mismatch_count=0,first_mismatch=None).items():assert result[name]==value,name
            else:
                for name,value in dict(selection_rows=23,capacity=32,header_values=6,layouts=2,
                        reference_rows=138,scalar_packets=360,raw_packets=24,pure_calls=21199,
                        other_program_word_refusals=21142,boundary_refusals=21,
                        compared_output_words=5189312960,compared_input_words=601713080,
                        whole3829_words_and_inactive_capacity_and_canaries=True,input_immutable=True,
                        MXCSR_controls_unchanged=True,caller_MXCSR_state_restored=True).items():assert result[name]==value,name
    assert code['public_guards_still_reference_actual_admission_or_memcmp']
    assert code['private_export_direct_edges_skip_program_comparisons_only']
    assert code['diagnostic_real_preflight_and_private_symbol_references']
    roles=('transport','endpoint','dopri_phase','rows_adapter')
    production=build['production_objects'];assert [row['role'] for row in production]==list(roles)
    assert production==manifest['production_objects']
    objects={row['role']:row['object'] for row in production}
    for role,row in zip(roles,production):
        obj=row['object'];checked(R/obj['path'],obj)
        assert {k:obj[k] for k in ('bytes','sha256')}==manifest['products'][obj['path']]
        code_row,=[item for item in code['objects'] if item['role']==role]
        assert code_row['object']==obj and code_row['baseline_legacy_encodings_only']
        assert code_row['no_FMA_x87_or_FP_control_writes']
    old=doc(R/'out/native-cpu-inline-private-provider/provider/source-contract.json')
    native_first7=[objects[stage] for stage in STAGES]+old['native_objects'][3:7]
    math=[row for row in parent['parts'] if Path(row['path']).name in ('dispatch.o','v3.o','optional.o')]
    assert len(math)==3
    parts=[objects['rows_adapter']]+parent['raw_objects']+math
    assert len(parts)==7
    for row in native_first7+parts:checked(R/row['path'],row)
    for stage in STAGES:
        assert derivation['scalar_source_inverse'][stage]['public_wrapper_exact_source_prefix']
        assert derivation['scalar_source_inverse'][stage]['private_body_inverse_restores_public_function']
        assert derivation['adapter_source_inverse'][stage]['private_inverse_restores_public_function']
        assert derivation['adapter_source_inverse'][stage]['full8_v3_v4_ISA_and_scalar_tail_selection_unchanged']
    patch=files[prefix+'mesa.patch'];assert data_pin(patch)==derivation['derived']['mesa.patch']
    header_basis={name:derivation['source_inputs'][str((S/name).relative_to(R))] for name in HEADERS}
    for name,row in derivation['derived'].items():
        assert data_pin(files[prefix+name])==row,name
    return dict(parts=parts,native_first7=native_first7,raw_objects=parent['raw_objects'],
        adapter=objects['rows_adapter'],scalar_objects=[objects[s] for s in STAGES],
        manifest=MANIFEST,archive=archive,math_manifest=parent['math_manifest'],
        math_archive=parent['math_archive'],math_code=parent['math_code'],math_code_proof=parent['math_code_proof'],
        code=parent['code'],admission_code=code,finite_results=compatibility['fixtures'],
        mesa_patch=patch.decode(),header_basis=header_basis,
        guard_shared_byte_authority='out/native-cpu-dispatch-admission-token/native_packet.hpp',
        source_derivation=derivation,parent_raw_manifest=parent['manifest'],parent_raw_archive=parent['archive'])


