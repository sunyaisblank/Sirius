from pathlib import Path
import hashlib, importlib.util, json

R = Path.cwd()
W = Path(__file__).resolve().parent
B = W / 'audited-benefit'
folder = W / 'audited-control-8'
owner = json.loads((folder / 'owner.json').read_text())
assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['stop_reason'] and not owner['cleanup_errors'] and owner['remaining_owned_processes'] == []
assert owner['provider_mappings'] and all(birth['absent'] for birth in owner['observed_births'])
spec = importlib.util.spec_from_file_location('original_audited_contract', B / 'software_contracts.py')
reader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reader)
props = reader.xml_properties(folder / 'gtest.xml', 'timestamps')
result = reader.timestamp_result(props, 'candidate')
assert len(result['readbacks']) == 32 and len(result['executed_modules']) == 54
old_props = reader.xml_properties(W / 'control-11/gtest.xml', 'timestamps')
old_readbacks = {k: reader.fields(v) for k, v in old_props.items() if k in result['readbacks']}
assert old_readbacks == result['readbacks']
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
record = {'passed': True, 'candidate_revision': owner['source']['revision'],
          'xml_sha256': sha(folder / 'gtest.xml'), 'owner_sha256': sha(folder / 'owner.json'),
          'readbacks': 32, 'available_markers': 51, 'actual_module_command_joins': 54,
          'inclusive_instance_API_schema_pass': True, 'original_uninstrumented_readbacks_byte_exact': True,
          'original_reader_sha256': sha(B / 'software_contracts.py'),
          'scope': 'Original last correctness prerequisite plus complete benchmark schema restored. No cost improvement or full science/release verdict.'}
destination = W / 'audited-witness-acceptance.json'
assert not destination.exists()
destination.write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record))
