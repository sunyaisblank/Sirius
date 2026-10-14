"""Compare durable JSON reports without changing consumed host/capture inputs."""
from pathlib import Path
import copy
import hashlib
import importlib.util
import json

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
HOST = ROOT / 'out/portable-divide64-jit-preparation'
spec = importlib.util.spec_from_file_location('consumed_host', HOST / 'host.py')
host = importlib.util.module_from_spec(spec); spec.loader.exec_module(host)

def canonical_json(value):
    # The original analyzer uses integer relocation keys; JSON object names
    # are strings. Normalize through its actual persistence format, retaining
    # every field, number type, sequence order and value.
    normalized = json.loads(json.dumps(value, allow_nan=False))
    return json.dumps(normalized, sort_keys=True, separators=(',', ':'), allow_nan=False)

def main():
    assert Path.cwd().resolve() == host.ROOT and not (WORK / 'pair-review.json').exists()
    host.verify_static()
    layout = json.loads((HOST / 'feature-layout.json').read_text())
    assert len(layout['features_fields']) == len(set(layout['features_fields'])) == 55
    original = host.module('original_pair_jit_analysis', HOST / 'analyze_capture.py')
    reader = host.module('paired_device_features', HOST / 'read_features.py')
    normalized, captures, controls, correction = {}, {}, [], {}
    for label in ('capture-stock-endpoint-normal', 'capture-candidate-endpoint-normal'):
        host.accepted(label)
        folder = HOST / label
        code = json.loads((folder / 'actual-code-review.json').read_text())
        recomputed = original.inspect(folder)
        assert canonical_json(recomputed) == canonical_json(code)
        correction[label] = {'native_python_equal': recomputed == code,
                             'complete_durable_json_equal': True,
                             'native_relocation_key_types': sorted({type(k).__name__ for k in recomputed['linked_section_review']['relocation_counts']}),
                             'json_relocation_key_types': sorted({type(k).__name__ for k in code['linked_section_review']['relocation_counts']})}
        for mutation in ('changed-count', 'missing-count', 'string-count'):
            altered = copy.deepcopy(code)
            counts = altered['linked_section_review']['relocation_counts']; key = next(iter(counts))
            if mutation == 'changed-count': counts[key] += 1
            elif mutation == 'missing-count': del counts[key]
            else: counts[key] = str(counts[key])
            assert canonical_json(recomputed) != canonical_json(altered)
            controls.append({'capture': label, 'mutation': mutation, 'rejected': True})
        features = json.loads((folder / 'feature-review.json').read_text())
        assert reader.inspect(folder) == features
        assert set(features['features']) == set(layout['features_fields'])
        features['features'] = {name: features['features'][name] for name in layout['features_fields']}
        assert features['enabled_extension_count'] == len(features['enabled_extensions'])
        normalized[label] = features
        files = sorted(path for path in folder.rglob('*') if path.is_file())
        assert files and all(not path.is_symlink() for path in files)
        captures[label] = {str(path.relative_to(folder)): host.record(path) for path in files}
    stock = normalized['capture-stock-endpoint-normal']
    assert stock == normalized['capture-candidate-endpoint-normal'] and len(controls) == 6
    host.verify_static()
    result = {'passed': True, 'normalized_feature_report': stock,
              'complete_55_enabled_fields_queue_extensions_order_and_feature_chain_equal': True,
              'captures': captures, 'serialization_correction': correction, 'negative_controls': controls,
              'consumed_host': host.record(HOST / 'host.py'), 'reader': host.record(Path(__file__).resolve()),
              'capture_reexecution_performed': False,
              'scope': 'Two original terminal short-input owners/readers with identical admitted features and complete durable JSON identity. No retained arithmetic, performance benefit, adoption or qualification verdict.'}
    (WORK / 'pair-review.json').write_text(json.dumps(result, indent=2) + '\n')
    print('PASS: both complete durable reports, identical feature policy, six serialization integrity negatives')

if __name__ == '__main__':
    main()
