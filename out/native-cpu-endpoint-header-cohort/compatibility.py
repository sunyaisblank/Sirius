from rebind import ao
exec(compile(ao('compatibility.py'), '/home/astra/.project/Sirius/out/native-cpu-v4-cohort/compatibility.py', "exec"))

_transport_check = check


def check(root):
    result = _transport_check(root)
    work = root / 'out/native-cpu-endpoint-header-domain'
    receipt = json.loads((work / 'compatibility.json').read_text())
    contract = json.loads((root / 'out/native-cpu-endpoint-header-provider/provider/source-contract.json').read_text())
    for name in ('compatibility.json', 'code-proof.json', 'build-result.json'):
        path = work / name
        identity = dict(bytes=path.stat().st_size, sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        assert identity == contract['source_inputs'][str(path.relative_to(root))]
    assert receipt['passed'] and receipt['source_unchanged'] and receipt['inputs_unchanged']
    assert receipt['objects'] == [contract['native_objects'][i] for i in (3, 1)]
    comparison = receipt['comparison']
    assert comparison['passed']
    for key, expected in dict(selection_rows=23, capacity=32, header_values=6, layouts=2,
            reference_rows=138, scalar_packets=360, raw_packets=24,
            other_program_word_refusals=21142, boundary_refusals=21).items():
        assert comparison[key] == expected, key
    for key in ('whole3829_words_and_inactive_capacity_and_canaries', 'input_immutable',
                'MXCSR_controls_unchanged', 'caller_MXCSR_state_restored'):
        assert comparison[key]
    result['inputs'][str((work / 'compatibility.json').relative_to(root))] = contract['source_inputs'][str((work / 'compatibility.json').relative_to(root))]
    result['scope']['focused_Endpoint_header_comparison'] = comparison
    return result
