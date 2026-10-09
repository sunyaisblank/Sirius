#!/usr/bin/env python3
"""Offline review of preserved capture; never launches a target or debugger."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ACTUAL=HERE/'actual'
def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()
pre=json.loads((HERE/'profile-preflight.json').read_text())
report=json.loads((ACTUAL/'report.json').read_text())
snapshot=json.loads((ACTUAL/'snapshot.json').read_text())
assert sha(ACTUAL/'snapshot.json')==report['snapshot_sha256']
assert sha(HERE/'controller.py')==pre['controller_sha256']==report['controller_sha256_after']
paths=[];deleted=[]
for line in snapshot['maps'].splitlines():
    fields=line.split(maxsplit=5)
    if len(fields)!=6 or not fields[5].startswith('/'):continue
    if fields[5].endswith(' (deleted)'):
        deleted.append(line)
        continue
    paths.append((Path(fields[5]).resolve(strict=True),line))
providers=[]
for item in pre['driver_binding_before']['libraries']:
    path=Path(item['path']).resolve(strict=True)
    lines=[line for mapped,line in paths if mapped==path]
    providers.append({**item,'snapshot_mapping_lines':lines,'post_capture_current_sha256':sha(path),
                      'mapped_path_matches_prebound_runtime_path':bool(lines),'post_capture_hash_matches_prebound_hash':sha(path)==item['sha256']})
frames_by_thread=[]
for entry in snapshot['stacks']:
    frames=[v['frame'] for v in entry['response']['fields'].get('stack',[])]
    frames_by_thread.append({'thread_id':entry['thread_id'],'target_id':entry['target_id'],'status':entry['response']['class'],
                             'frames':frames,'LLVM_frame_count':sum('libLLVM' in v.get('from','') for v in frames)})
active=[entry for entry in frames_by_thread if entry['frames'] and 'futex' not in entry['frames'][0].get('func','')]
summary={
    'scope':'Offline supplemental reading of the one saved compiler profiling snapshot. Does not modify or turn the original failed collector receipt into a complete receipt.',
    'original_report_sha256':sha(ACTUAL/'report.json'),'snapshot_sha256':sha(ACTUAL/'snapshot.json'),
    'preflight_sha256':sha(HERE/'profile-preflight.json'),'supplement_source_sha256':sha(__file__),
    'original_collector_disposition':report['disposition'],'original_error':report['error'],
    'error_boundary':'snapshot.json already written and snapshot_count=1 before mapped_paths rejected deleted mappings; runtime verification failed after capture, not interrupt or stack acquisition.',
    'unrelated_deleted_mappings':deleted,
    'deleted_mapping_exclusion':'Only the supplemental path matcher excludes deleted mappings. All are preserved verbatim here; the failed original collector and report remain unchanged.',
    'required_provider_mappings':providers,
    'all_required_provider_paths_and_current_post_hashes_match':all(v['mapped_path_matches_prebound_runtime_path'] and v['post_capture_hash_matches_prebound_hash'] for v in providers),
    'binding_time_limit':'Hashes were recorded from the actual provider files before launch and again after termination in the original receipt. Saved stopped-process maps name these files. This supplement adds offline map/path matching; it does not invent an in-snapshot hash or a completed original validation step.',
    'all100_original_frozen_source_artifact_preparation_and_link_bindings_unchanged':report['all_frozen_bindings_unchanged'] and len(pre['frozen_bindings'])==100,
    'driver_ICD_and_LLVM_before_after_unchanged':report['driver_bindings_unchanged'],
    'debugger_and_controller_before_after_unchanged':report['debugger_and_controller_unchanged'],
    'all_threads_listed':len(snapshot['threads']),'all_threads_stopped':all(v['state']=='stopped' for v in snapshot['threads']),
    'captured_stacks':len(snapshot['stacks']),'all_stacks_successful_and_nonempty':all(v['status']=='done' and v['frames'] for v in frames_by_thread),
    'omitted_stack_thread_ids':snapshot['omitted_stack_thread_ids'],
    'LLVM_module_frames_in_saved_stacks':sum(v['LLVM_frame_count'] for v in frames_by_thread),
    'non_futex_top_frame_threads':active,
    'main_thread_scientific_call_chain':[v['func'] for v in frames_by_thread[0]['frames'] if any(key in v.get('func','') for key in ['RetainedCompute','VulkanDevice','TestBody'])],
    'phase_finding':'At this one stopped instant, main thread waits in RetainedCompute::Step -> VulkanDevice::Dispatch; the only non-futex top-frame thread is LWP3265 in libc memset with stripped libvulkan_lvp.so caller frames. No LLVM-module frame appears in any captured stack. Exact Mesa compiler/pass phase remains unresolved; absence of LLVM frames at one instant does not establish that LLVM was never reached.',
    'smaps_rollup':snapshot['smaps_rollup'],
    'trigger_elapsed_seconds':report['trigger_elapsed_seconds'],'trigger_inferior_rss_kib':report['trigger_inferior_rss_kib'],
    'peak_sampled_debugger_and_inferior_tree_rss_kib':report['peak_sampled_debugger_and_inferior_tree_rss_kib'],
    'peak_sampled_inferior_rss_kib':report['peak_sampled_inferior_rss_kib'],'minimum_sampled_available_kib':report['minimum_sampled_available_kib'],
    'terminal_cleanup':report['cleanup'],'no_resource_guard_triggered':report['guard_reason'] is None,
    'qualification':'No numerical verdict, no readback or XML completion, no repeated observation, no isolated performance claim. Native Radeon P1 concurrency recorded in original receipt.',
}
xml=ACTUAL/'evidence/transport.xml'
summary['terminal_XML_present']=xml.exists()
summary['numerical_readback_files']=list(str(v.relative_to(ACTUAL)) for v in (ACTUAL/'evidence').rglob('*') if v.name in ['scientific-output.bin','production-readback.bin'])
assert not xml.exists() and not summary['numerical_readback_files']
assert summary['all_required_provider_paths_and_current_post_hashes_match']
assert summary['all_stacks_successful_and_nonempty'] and summary['all_threads_stopped'] and not summary['omitted_stack_thread_ids']
assert summary['terminal_cleanup']['inferior_identity_gone'] and not summary['terminal_cleanup']['remaining_owned_group_or_session_members'] and summary['terminal_cleanup']['all_streams_eof']
(ACTUAL/'supplemental-review.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({'supplement_sha256':sha(ACTUAL/'supplemental-review.json'),'original_report_unchanged_sha256':sha(ACTUAL/'report.json'),
                  'snapshot_unchanged_sha256':sha(ACTUAL/'snapshot.json'),'thread_count':summary['all_threads_listed'],'all_provider_mappings_match':True,
                  'exact_phase':'unresolved with installed stripped symbols','no_owned_processes':True}))
