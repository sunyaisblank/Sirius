"""Read sealed actual RGP event metadata and unavailable instruction samples."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, re, subprocess
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))
def seal(path):
    return {'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
bridge = read(WORK / 'profiler-interval-second-pair-bridge.json')
acquisition = read(WORK / 'native-capture-interval-console-bridge.json')
assert bridge['accepted'] and acquisition['accepted']
assert bridge['returncode'] == 0 and not bridge['outer_stop'] and not bridge['errors']
assert all(seal(Path(p) if Path(p).is_absolute() else ROOT / p) == s
           for p, s in bridge['input_seals'].items())
output = WORK / 'profiler-interval-second-pair'
owner = read(output / 'owner.json')
assert owner['owned_process_absent'] and not owner['cleanup_errors']
assert owner['intentional_owned_stop'] and owner['trace_unchanged'] and owner['tool_unchanged']
assert owner['developer_census_completed'] and not owner['accessible_remaining_developer_processes']
trace = WORK / 'native-capture-interval-console/original-legacy-interval.rgp'
expected = acquisition['original_control_and_acquired_trace']
assert seal(trace) == {'bytes': expected['trace_bytes'], 'sha256': expected['trace_sha256']}
assert len(owner['event_observations']) == 4
events = []
evidence = [trace, WORK / 'profiler-interval-second-pair-bridge.json', output / 'owner.json']
for label, event_name, unavailable in (
        ('transport', '16 vkCmdDispatch(3, 1, 1)', 'Event 16 \nNo instruction timing data'),
        ('endpoint', '22 vkCmdDispatch(6, 1, 1)', 'Event 22 \nNo instruction timing data')):
    observations = [v for v in owner['event_observations'] if v['event_name'] == event_name]
    assert [v['view'] for v in observations] == ['event', 'instructions']
    assert all(v['selected_event_confirmed'] and v['selection_state']['fresh_is_selected']
               for v in observations)
    event_path = output / (label + '-event-automation.json')
    instruction_path = output / (label + '-instructions-automation.json')
    event = read(event_path)
    instruction = read(instruction_path)
    def field(suffix):
        values = [v['name'] for v in event
                  if v['automation_id'].endswith(suffix) and not v['offscreen']]
        assert len(values) == 1, (suffix, values)
        return values[0]
    assert field('SingleEventSidePanel.event_name_') == event_name
    assert any(v['name'] == unavailable and not v['offscreen'] for v in instruction)
    shader = [v['name'] for v in event
              if re.fullmatch(r'CS: 0x[0-9A-F]+', v['name']) and not v['offscreen']]
    assert len(shader) == 1
    events.append({'label': label, 'actual_rgp_event': event_name,
                   'start_time_display': field('SingleEventSidePanel.start_time_'),
                   'end_time_display': field('SingleEventSidePanel.end_time_'),
                   'duration_display': field('SingleEventSidePanel.duration_'),
                   'work_duration_display': field('SingleEventSidePanel.work_duration_'),
                   'total_wavefronts_display': field('total_wavefronts_'),
                   'total_threads_display': field('total_threads_'),
                   'average_threads_per_wave_display': field('average_threads_'),
                   'api_shader_hash_display': shader[0],
                   'api_pso_hash_display': field('api_pso_hash_value_'),
                   'driver_pipeline_hash_display': field('internal_pipeline_hash_value_'),
                   'usable_instruction_samples': False,
                   'instruction_pane_message': unavailable})
    evidence.extend([event_path, instruction_path])
    for observation in observations:
        for window in observation['window_observation']:
            image = output / window['image_path'].replace('\\', '/').rsplit('/', 1)[-1]
            actual_windows_path = subprocess.check_output(
                ['wslpath', '-w', str(image.resolve())], text=True).strip()
            assert actual_windows_path.casefold() == window['image_path'].casefold()
            assert seal(image) == {'bytes': window['image_bytes'],
                                   'sha256': window['image_sha256']}
            evidence.append(image)
record = {'capture_acquired': True, 'gui_observation_accepted': True,
          'usable_executed_instruction_cost': False, 'events': events,
          'all_gui_frozen_inputs_rechecked_exact': len(bridge['input_seals']),
          'scope': 'Actual selected trace events and displayed identifiers are observed. RGP event indices16/22 are not process dispatch ordinals. Range and dimensions support intended legacy Transport/Endpoint but do not independently prove source ordinals or marker indices. Only this selected pair was observed; no whole-trace absence-of-samples claim. API/driver hashes are not SPIR-V SHA256. No analyzed instruction cost, exclusive arithmetic cost, ordinary timing benefit, frame or full qualification.',
          'evidence': {str(p.relative_to(ROOT)): seal(p) for p in evidence}}
destination = WORK / 'interval-second-pair-result.json'
assert not destination.exists()
destination.write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({k: record[k] for k in ('capture_acquired', 'gui_observation_accepted',
                                       'usable_executed_instruction_cost')}))
