"""Pin the exercised AP producer chain; bind only AT sources and raw symbols."""
from pathlib import Path
import hashlib

_AS_SOURCE = Path.cwd() / 'out/native-cpu-inline-private-cohort'
_AS_PINS = {'compatibility.py': 'dce4e1ab8a4eca768cec3dbf6b09b7cbf5e0de8fc10e7969bfccec085fd7e296', 'native_observer.py': '1268f40f3e085e9cff62cb9bdb9c3f6a10c5ec75f86c8b31888cb93fbc2f299c', 'optional_math.py': 'd3fd4d943fffd2f83534c7e3582bf045ead2c86aa6e44bc6b9f4f140669d041b', 'optional_observer.py': '6526a8455c67b3330093129394d8cf740c5caaeaaf7502f05b8afb4514e17523', 'owned_guard.py': 'ab5a3e58e47311e363842cfa04715f475483ee05f8bd63391a4b8a67787c0d92', 'prepare.py': '118cce14e369f5844c9435d5e85d5a63e67926df6ec38028a96c62884060d28d', 'read_result.py': '94e821db4f0e14db2daf09bdc5d5a8267b3247d40c779b04a8d7f33390c63d26', 'rebind.py': 'dc0372d2db801498d063723c7b5cca6e99f36f685dad835df69800676b4008f5', 'run_case.py': 'c8aa0b3f3622a278f1edc3839d4334e6c1dc61d814b68fa5425a6609dd7d20fd', 'run_owned.py': '7d832b5c9edc5518e0aca52e1325368744ae40215bf7978f4ead4ab2740672ef', 'runtime_transition.py': 'd3c3b1a3c6115c55e791132e45c793c64c7c0c26b2aa9980382d40ae880dde19'}
for _name, _digest in _AS_PINS.items():
    assert hashlib.sha256((_AS_SOURCE / _name).read_bytes()).hexdigest() == _digest, _name

_AR_SOURCE = Path.cwd() / 'out/native-cpu-packed-math-cohort'
_AR_PINS = {'compatibility.py': 'f2ae0fb85f8e6646fdee2a271bb30b0eb9b0ca38391ce4ca263989f6fd96b488', 'native_observer.py': 'd351a37aded9549738e05548ed1e3a983d3d50de6256b324dc0a7f298d559190', 'optional_math.py': 'd3fd4d943fffd2f83534c7e3582bf045ead2c86aa6e44bc6b9f4f140669d041b', 'optional_observer.py': '6526a8455c67b3330093129394d8cf740c5caaeaaf7502f05b8afb4514e17523', 'owned_guard.py': 'ab5a3e58e47311e363842cfa04715f475483ee05f8bd63391a4b8a67787c0d92', 'prepare.py': '118cce14e369f5844c9435d5e85d5a63e67926df6ec38028a96c62884060d28d', 'read_result.py': '94e821db4f0e14db2daf09bdc5d5a8267b3247d40c779b04a8d7f33390c63d26', 'rebind.py': 'c45fed40589d386f11721b931e155e291a6dae6567a6a758a3a174dda6dd78aa', 'run_case.py': 'c8aa0b3f3622a278f1edc3839d4334e6c1dc61d814b68fa5425a6609dd7d20fd', 'run_owned.py': '7d832b5c9edc5518e0aca52e1325368744ae40215bf7978f4ead4ab2740672ef', 'runtime_transition.py': 'd3c3b1a3c6115c55e791132e45c793c64c7c0c26b2aa9980382d40ae880dde19'}
for _name, _digest in _AR_PINS.items():
    assert hashlib.sha256((_AR_SOURCE / _name).read_bytes()).hexdigest() == _digest, _name

R = Path.cwd()
AP = R / 'out/native-cpu-endpoint-header-cohort'
AP_PINS = {'run_owned.py': '7d832b5c9edc5518e0aca52e1325368744ae40215bf7978f4ead4ab2740672ef', 'runtime_transition.py': 'd3c3b1a3c6115c55e791132e45c793c64c7c0c26b2aa9980382d40ae880dde19', 'optional_observer.py': '6526a8455c67b3330093129394d8cf740c5caaeaaf7502f05b8afb4514e17523', 'read_result.py': '94e821db4f0e14db2daf09bdc5d5a8267b3247d40c779b04a8d7f33390c63d26', 'optional_math.py': '873b7eed86d2234210d146d80369ac9cd622dc031b545dbc6daaf34f7f84ca5f', 'rebind.py': '420331da24d7e5768532d2466a6a1d7b80036579204c38aa2f82ef0af4371c95', 'run_case.py': 'c8aa0b3f3622a278f1edc3839d4334e6c1dc61d814b68fa5425a6609dd7d20fd', 'native_observer.py': 'e3819bb87b6133f30d98cb1a4bc4623a7bedf5d9e41de8d715383722c191fccf', 'compatibility.py': '2d8b51030d408bca810ce94b586f07258113fca12b45cfc0a9872a0116b0c224', 'prepare.py': '118cce14e369f5844c9435d5e85d5a63e67926df6ec38028a96c62884060d28d', 'owned_guard.py': 'ab5a3e58e47311e363842cfa04715f475483ee05f8bd63391a4b8a67787c0d92'}
for name, digest in AP_PINS.items():
    assert hashlib.sha256((AP / name).read_bytes()).hexdigest() == digest, name
_AP = {}
exec(compile((AP / 'rebind.py').read_text(), str(AP / 'rebind.py'), 'exec'), _AP)


def adapt(text, name):
    for before, after in (
        ('out/native-cpu-endpoint-header-cohort', 'out/native-cpu-dispatch-token-cohort'),
        ('out/native-cpu-endpoint-header-provider/provider', 'out/native-cpu-dispatch-token-provider/provider'),
        ('out/native-cpu-endpoint-header-compatibility', 'out/native-cpu-v4-raw-packets'),
        ('providerAP', 'providerAT'), ('current_AP', 'current_AT'),
        ('Same original H host/AP Endpoint admission provider', 'Same original H host/AT guarded raw provider')):
        text = text.replace(before, after)
    # Preserve actual observed function names; adjust only the new symbol prefix
    # and inverse stage/parent mapping, never normalize a fallback observation.
    text = text.replace("'sirius_simd_' + slug + '_packet'", "'sirius_simd_v4_' + slug + '_packet'")
    text = text.replace("len('sirius_simd_')", "len('sirius_simd_v4_')")
    text = text.replace(".replace('sirius_simd_', 'sirius_native_')", ".replace('sirius_simd_v4_', 'sirius_native_')")
    if name == 'optional_observer.py':
        before = "for index, slug in enumerate(('transport', 'endpoint', 'dopri_phase'), 4):\n    item = proof['objects'][index]"
        after = "for slug, item in zip(('transport', 'endpoint', 'dopri_phase'), proof['raw_v4_objects']):"
        assert text.count(before) == 1
        text = text.replace(before, after, 1)
    if name == 'optional_math.py':
        before = "raw = dict(zip(('transport', 'endpoint', 'dopri_phase'), proof['objects'][4:7]))"
        after = "raw = dict(zip(('transport', 'endpoint', 'dopri_phase'), proof['raw_v4_objects']))"
        assert text.count(before) == 1
        text = text.replace(before, after, 1)
    text = text.replace(
        'AP changes only Endpoint pure/scalar admission for relative header word11241 and provider identity; the unchanged shader retains its acceptance/refusal and all output behavior. Other six AO objects, actual v3 eligibility, FP/precision properties, raw3/scheduler/fallback and original physical contracts remain.',
        'AT replaces baseline scalar3/adapter and adds the exact existing-job-positive-preflight branch; private workers omit only repeated immutable-program comparisons. Public guards and all non-program guards remain. Pure admission, AS raw3, AVX2 fallback and math3 remain exact; physical/FP/scientific/count/scheduler contracts remain. C8 observes actual hidden admitted names and genuine job/views at first private entries.')
    return text


def ao(name):
    return adapt(_AP['ao'](name), name)


def load(name):
    return adapt(_AP['load'](name), name)


def authority_files():
    return _AP['authority_files']() + [AP / name for name in AP_PINS] + [_AR_SOURCE / name for name in _AR_PINS]


def seal_authorities(text):
    anchor = 'files = set(expected)'
    assert text.count(anchor) == 1
    return text.replace(anchor, 'for path in authority_files():\n    expected[str(path)] = pin(path)\n' + anchor, 1)


_previous_AS_authorities = authority_files
def authority_files():
    return _previous_AS_authorities() + [_AS_SOURCE / name for name in _AS_PINS]


_previous_adapt = adapt

def adapt(text, name):
    text = _previous_adapt(text, name)
    if name in ('native_observer.py', 'read_result.py', 'run_case.py'):
        # These are prospective symbol/expectation substitutions in source.
        # Captured inferior function names are never rewritten.
        text = text.replace('_rows_try_packet', '_rows_admitted_packet')
        text = text.replace("'_try_packet' if category == 'scalar_calls'", "'_admitted_packet' if category == 'scalar_calls'")
        text = text.replace('sirius_native_transport_try_packet', 'sirius_native_transport_admitted_packet')
    if name == 'native_observer.py' and "class Entered(gdb.Breakpoint):" in text:
        extension = (R / 'out/native-cpu-dispatch-token-cohort/dispatch_observer.py').read_text()
        anchor = 'class Entered(gdb.Breakpoint):'
        assert text.count(anchor) == 1
        text = text.replace(anchor, extension + '\n\n' + anchor, 1)
        anchor = "            observation[self.category].append(row)"
        assert text.count(anchor) == 1
        addition = "            if self.category in ('calls', 'simd_parent_calls', 'scalar_parent_calls'):\n                capture_dispatch_job(row, self.slug)\n"
        text = text.replace(anchor, addition + anchor, 1)
        anchor = "                    parent = parents[(row['thread'], self.slug, 'scalar')]"
        assert text.count(anchor) == 1
        text = text.replace(anchor, anchor + "\n                    row['dispatch_job'] = dict(parent['dispatch_job'])", 1)
    if name == 'read_result.py' and 'proof = ' in text:
        text = text.replace("endswith(('_try_packet', '_admit_packet'))", "endswith(('_try_packet', '_admit_packet', '_admitted_packet'))")
        anchor = '# Changed R Transport SIMD has fresh component execution and exact new-code joins.'
        assert text.count(anchor) == 1
        checks = '''
# Actual hidden worker identity and genuine positive job are independent of
# public preflight evidence and of final static branch/call-target proofs.
private_rows = observation['calls'] + observation['simd_parent_calls'] + observation['scalar_parent_calls']
assert private_rows
for row in private_rows:
    job = row['dispatch_job']
    slug = row['function'][len('sirius_native_'):-len('_rows_admitted_packet')]
    assert slug in slugs and row['function'] == 'sirius_native_' + slug + '_rows_admitted_packet'
    assert job['native_preflight_passed'] is True and job['preflight'] is False
    assert job['stage'] == slugs.index(slug) + 1 and job['job'] > 0
    assert job['caller_function'] == 'lp_native_cpu_packet_mode'
    assert job['job_type_bytes'] == 112 and job['descriptor_type_bytes'] == 64
    assert job['callback_function'] == 'cs_exec_fn'
    static_proof=json.loads((root/'out/native-cpu-dispatch-token-provider/provider/linked-native-proof.json').read_text())
    assert job['layout']==static_proof['job_observation_layout']
    assert job['register_evidence']==static_proof['callback_live_job_register']
    base=checked_entry(row)
    assert job['callback_pc']==base+job['register_evidence']['return_address']
    actual_call,=[entry for entry in static_proof['actual_hidden_admitted_worker_calls'] if entry['target']==row['function']]
    assert job['caller_pc']==base+actual_call['return_address']
    for key in ('thread','grid_x','grid_y','grid_z','input','input_bytes','output','output_bytes'):
        assert job[key] == row[key], key
for row in observation['scalar_calls']:
    parent = scalar_parents[row['parent_entry_sequence']]
    assert row['dispatch_job'] == parent['dispatch_job']
dispatch_admission_scope = dict(all_observed_private_rows_have_positive_live_job=True,
    private_rows_observed=len(private_rows), actual_functions=sorted({row['function'] for row in private_rows}),
    job_and_descriptor_views_match_actual_packet_arguments=True,
    source='Existing first-entry stops only; genuine unwound callee-saved job register and compile-time asserted layout, no extra workload or inferior calls/writes')
'''
        text = text.replace(anchor, checks + '\n' + anchor, 1)
        anchor = "'native_first_packets_handled': 3, 'pure_preflights_handled': 3,"
        assert text.count(anchor) == 1
        text = text.replace(anchor, "'native_first_packets_handled': 3, 'dispatch_admission': dispatch_admission_scope, 'pure_preflights_handled': 3,", 1)
    return text
