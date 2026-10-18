import hashlib, importlib.util, json, pathlib, sys

sys.dont_write_bytecode = True
ROOT = pathlib.Path.cwd()
BASE = ROOT / 'attestations/qualification/5b2df58b/native-constant-controls-benefit-pass/files/out'
reader = BASE / 'native-constant-pe-audit/audit_export_windows_paths.py'
spec = importlib.util.spec_from_file_location('saved_native_PE_reader', reader)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
candidate = BASE / 'native-constant-producer-inputs/candidate'
baseline = BASE / 'native-constant-producer-inputs/baseline'
expected = {f'main::{n}': b for n, b in module.arrays(candidate / 'retained_kernels.h').items()}
expected.update({f'owner::{n}': b for n, b in module.arrays(candidate / 'retained_owner_kernels.h').items()})
assert len(expected) == 42
old = module.arrays(baseline / 'retained_kernels.h')
changed = ['kTransportShader', 'kEndpointShader', 'kDenseShader', 'kInitializeShader', 'kDopriPhaseShader']
stage = pathlib.Path(sys.argv[1])
results = {}
for name in ['sirius.exe', 'sirius_render_tests.exe']:
    exe = stage / 'src/sirius/app/Release' / name
    image_base, sections = module.pe_sections(exe)
    joins = {}
    for array, payload in expected.items():
        locations = []
        for section, rva, raw, data in sections:
            start = 0
            while True:
                offset = data.find(payload, start)
                if offset < 0:
                    break
                locations.append(dict(section=section, file_offset=raw+offset, rva=rva+offset))
                start = offset + 1
        assert locations, (name, array, 'missing complete expected payload')
        joins[array] = dict(bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest(), locations=locations)
    for n in changed:
        assert all(old[n] not in data for _, _, _, data in sections), (name, n, 'old native payload also present')
    b = exe.read_bytes()
    results[name] = dict(path=str(exe), bytes=len(b), sha256=hashlib.sha256(b).hexdigest(),
                         main40_and_owner2_full_payloads=joins, old5_native_payloads_absent=True)
out = dict(schema='sirius-native-constant-canonical-frame-PE-payloads-v1', passed=True,
           source_revision='2c8ef08693d46f57d5dea14efdbddd4874dac684', stage=str(stage), PEs=results,
           admitted_candidate_evidence_commit='fd530d957d2e0737f5453d88d78d979f8ed17607',
           reused_reader=dict(path=str(reader.relative_to(ROOT)), sha256=hashlib.sha256(reader.read_bytes()).hexdigest()),
           scope='Complete expected main40/owner2 bytes in actual x64 PE initialized readable non-code sections; old5native modules absent. No MAP symbol/relocation proof or runtime execution. Combine actual CI source/gate/stage joins and unchanged selector source before frame release.',
           runtime_admitted=False, full_qualification_or_release=False)
pathlib.Path(sys.argv[2]).write_text(json.dumps(out, indent=2) + '\n')
print('PASS: canonical CLI/render PEs contain all42 admitted complete payloads; old5native absent')
