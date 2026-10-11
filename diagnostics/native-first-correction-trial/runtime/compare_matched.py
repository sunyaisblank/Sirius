"""Compare the original native controls; no retention threshold is applied."""

import hashlib
import json
import math
from pathlib import Path
import re
import sys
sys.dont_write_bytecode = True
import xml.etree.ElementTree as ET


WORKSPACE = Path(__file__).resolve().parent
ROOT = WORKSPACE.parent.parent.parent
LABELS = ("a1", "b1", "b2", "a2")
BASELINE_REVISION = "b81e061fa081cdaed37eca3add950bd254bf296d"
BASELINE_BACKEND = "5e01bad179c838022585d14bf6b7be6fa64f3fe2f070178970ca1aee40e494fa"
BASELINE_GATE = "96c1add302931977ced3168d9893a67c4a770e3a64054232d533a44104a03c64"
CASES = {
    "factored": "RetainedComputeTest.FactoredEndpointsPreserveIndependentRootsAndBoundaryAdmission",
    "timestamps": "RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults",
}
STAGES = ("camera", "transport", "endpoint", "dense", "initialize", "ray_camera", "dopri_phase")
FACTORED_ROWS = (29, 29, 29, 4, 8, 29, 27, 27, 23, 23, 23, 4, 1, 29)
ROW_WORDS = {"transport": 4699, "endpoint": 3829, "dopri_phase": 4149}
COUPLED_STREAM = (
    ("transport", 12), ("endpoint", 24), ("transport", 12), ("endpoint", 24),
    ("transport", 12), ("endpoint", 24), ("dopri_phase", 24), ("dopri_phase", 24),
    ("endpoint", 12), ("dopri_phase", 12), ("endpoint", 12),
)
MARKER_MS = ("device_ms", "host_submit_ms", "host_wait_ms", "host_completion_ms",
             "pipeline_ms", "setup_ms", "cleanup_ms", "total_ms")
COST_MS = ("device_ms", "host_submit_ms", "host_wait_ms", "host_completion_ms", "total_ms")
COUPLED_STAGES = ("transport", "endpoint", "dopri_phase")
CASES.update(json.loads((WORKSPACE/"native-cases.json").read_text())["cases"])
MAX_FILE_BYTES = 16 * 1024 * 1024


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def read(path):
    require(path.is_file() and path.stat().st_size <= MAX_FILE_BYTES, f"Missing/oversized input: {path}")
    payload = path.read_bytes()
    require(len(payload) <= MAX_FILE_BYTES, f"Input grew beyond bound: {path}")
    return payload


def identity(path):
    payload = read(path)
    return {"path": str(path.relative_to(ROOT)), "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest()}


def document(path):
    def invalid_constant(value):
        raise ValueError(f"Nonfinite JSON constant in {path}: {value}")
    return json.loads(read(path).decode("utf-8-sig"), parse_constant=invalid_constant)


def digest(value, size=64):
    require(isinstance(value, str) and re.fullmatch(f"[0-9a-f]{{{size}}}", value),
            f"Invalid {size}-digit identity: {value!r}")
    return value


def milliseconds(value, context):
    require(not isinstance(value, bool), f"Boolean timing: {context}")
    number = float(value)
    require(math.isfinite(number) and number >= 0, f"Invalid timing: {context}")
    return number


def total(values):
    result = math.fsum(values)
    require(math.isfinite(result), "Nonfinite aggregate")
    return result


def fields(value):
    result = {}
    for item in value.split(";"):
        key, separator, word = item.partition("=")
        require(separator and key not in result, f"Malformed/duplicate field: {value}")
        result[key] = word
    return result


def original_properties(folder, case):
    xml = ET.fromstring(read(folder / "gtest.xml"))
    require(int(xml.get("tests", "-1")) == 1, f"Wrong test count: {folder}")
    for node in [xml, *xml.findall(".//testsuite")]:
        for name in ("failures", "errors", "skipped", "disabled"):
            require(int(node.get(name, "0")) == 0, f"Nonzero {name}: {folder}")
    require(not any(xml.findall(".//" + tag) for tag in ("failure", "error", "skipped")),
            f"Failed/incomplete test: {folder}")
    tests = xml.findall(".//testcase")
    require(len(tests) == 1 and tests[0].get("classname") + "." + tests[0].get("name") == CASES[case],
            f"Wrong original filter: {folder}")
    require(tests[0].get("status") == "run" and tests[0].get("result") == "completed",
            f"Unexecuted original case: {folder}")
    properties = {}
    for item in xml.findall(".//property"):
        name = item.get("name")
        require(name is not None and name not in properties, f"Duplicate XML property: {folder}")
        properties[name] = item.get("value")
    return properties


def owner(folder, case):
    value = document(folder / "owner.json")
    require(value["status"] == "completed" and type(value["returncode"]) is int and value["returncode"] == 0,
            f"Owner did not complete: {folder}")
    require(value["outer_stop"] is False and not value.get("cleanup_errors"), f"Stopped/unclean owner: {folder}")
    for key in ("owned_process_absent", "numerical_control_passed", "test_executable_unchanged",
                "loaded_driver_modules_recorded", "source_build_gate_unchanged"):
        require(value[key] is True, f"Missing owner proof {key}: {folder}")
    require(value["expected_cases"] == [CASES[case]], f"Wrong owner case: {folder}")
    require([a for a in value["arguments"] if a.startswith("--gtest_filter=")] ==
            ["--gtest_filter=" + CASES[case]], f"Wrong actual owner filter: {folder}")
    for key in ("source_revision", "live_source_revision"):
        digest(value[key], 40)
    for key in ("source_build_gate_sha256", "test_executable_sha256"):
        digest(value[key])
    for key in ("elapsed_seconds", "peak_sampled_rss_bytes", "outer_guard_seconds", "rss_guard_bytes"):
        milliseconds(value[key], key)
    require(value["outer_guard_seconds"] > 0 and value["rss_guard_bytes"] > 0,
            f"Unbounded owner: {folder}")
    for key in ("owner_pid", "test_pid"):
        require(type(value[key]) is int and value[key] > 0, f"Invalid owner/process ID: {folder}")
    for key in ("owner_start_time_utc", "test_start_time_utc"):
        require(isinstance(value[key], str) and value[key], f"Missing owner/process birth: {folder}")
    for key, word in value["environment"].items():
        if key not in ("TEMP", "TMP"):
            require(word in (None, ""), f"Execution override {key}: {folder}")
    return value


def provider(folder):
    modules = document(folder / "driver-modules.json")
    require(isinstance(modules, list) and 2 <= len(modules) <= 16, f"Wrong provider module list: {folder}")
    normalized = []
    for module in modules:
        name = module["name"].casefold()
        require(name == "vulkan-1.dll" or "amdvlk" in name, f"Unexpected provider: {folder}")
        require(type(module["bytes"]) is int and module["bytes"] > 0, f"Invalid provider size: {folder}")
        normalized.append({"name": name, "path": module["path"].casefold(), "bytes": module["bytes"],
                           "sha256": digest(module["sha256"]), "version": module["version"]})
    normalized.sort(key=lambda item: item["name"])
    require(len({item["name"] for item in normalized}) == len(normalized), f"Duplicate provider: {folder}")
    require(sum(item["name"] == "vulkan-1.dll" for item in normalized) == 1 and
            any("amdvlk" in item["name"] for item in normalized), f"Loader/AMD provider absent: {folder}")
    return normalized


def readbacks(properties, case):
    prefix = "factored_endpoint_readback_" if case == "factored" else "coupled_readback_"
    expected = ({f"{prefix}0_{wide}_{call}" for wide in range(2) for call in range(14)}
                if case == "factored" else {f"{prefix}{i}" for i in range(11)})
    actual = {key for key in properties if key.startswith(prefix)}
    require(actual == expected, f"Wrong complete {case} readback keys/count")
    result = {}
    for key in sorted(actual):
        value = fields(properties[key])
        require(set(value) == ({"rows", "bytes", "sha256"} if case == "factored" else
                              {"stage", "bytes", "sha256"}), f"Wrong readback fields: {key}")
        if case == "factored":
            rows = FACTORED_ROWS[int(key.rsplit("_", 1)[1])]
            require(int(value["rows"]) == rows, f"Wrong original row count: {key}")
            stage = "endpoint"
        else:
            stage, rows = COUPLED_STREAM[int(key.rsplit("_", 1)[1])]
            require(value["stage"] == stage, f"Wrong readback stage: {key}")
        require(int(value["bytes"]) == rows * ROW_WORDS[stage] * 4, f"Partial complete buffer: {key}")
        digest(value["sha256"])
        result[key] = value
    return result


def factored_timings(folder):
    logs = "\n".join(read(folder / name).decode("utf-8-sig") for name in ("stdout.log", "stderr.log"))
    creates = re.findall(r"\[FactoredEndpoint\] create: portable=([01]);wide=([01]);decline_rte32=([01]);fma_enabled=([01])", logs)
    require(len(creates) == 2 and [(a, b, c) for a, b, c, _ in creates] ==
            [("0", "0", "0"), ("0", "1", "0")], "Wrong native Factored routes")
    starts = re.findall(r"\[FactoredEndpoint\] begin: call=(\d+);rows=(\d+);wide=([01]);decline_rte32=([01])", logs)
    ends = re.findall(r"\[FactoredEndpoint\] returned: call=(\d+);success=([01]);pipeline_ms=([^;\s]+);submit_wait_ms=([^;\s]+);total_ms=([^;\s]+)", logs)
    require(len(starts) == len(ends) == 28, "Incomplete original Factored call logs")
    selected = re.findall(r"\[FactoredEndpoint\] selected_stage=([^;\s]+);sha256=([0-9a-f]{64})", logs)
    endpoint = re.findall(r"\[FactoredEndpoint\] selected_endpoint_sha256=([0-9a-f]{64})", logs)
    require(len(selected) == 14 and len(endpoint) == 2, "Incomplete selected module logs")
    routes = []
    for wide in range(2):
        section = selected[7 * wide:7 * (wide + 1)]
        require(tuple(stage for stage, _ in section) == STAGES, "Wrong selected stage order")
        hashes = dict(section)
        require(endpoint[wide] == hashes["endpoint"], "Selected Endpoint mismatch")
        calls = []
        for call in range(14):
            i = wide * 14 + call
            require(starts[i] == (str(call), str(FACTORED_ROWS[call]), str(wide), "0"),
                    "Changed original Factored call order/rows")
            n, success, pipeline, wait, complete = ends[i]
            require(int(n) == call and success == "1", "Failed/out-of-order Factored call")
            calls.append({"call": call, "rows": FACTORED_ROWS[call],
                          "pipeline_ms": milliseconds(pipeline, "Factored pipeline"),
                          "wait_ms": milliseconds(wait, "Factored wait"),
                          "total_ms": milliseconds(complete, "Factored inclusive Dispatch")})
        routes.append({"route": "default" if wide == 0 else "wide_excluded_control", "wide": wide,
                       "fma_enabled": int(creates[wide][3]), "selected_stage_sha256": hashes,
                       "initial_pipeline_ms": calls[0]["pipeline_ms"],
                       "warm13_wait_ms": total(c["wait_ms"] for c in calls[1:]),
                       "warm13_total_ms": total(c["total_ms"] for c in calls[1:]), "all14_calls": calls})
    return routes


def timestamps(properties):
    require(properties["coupled_rows"] == "12" and properties["coupled_projection_budget"] == "24" and
            properties["product_mode"] == "default" and properties["retained_scalar_route"] == "native binary32",
            "Changed original modern cohort/mode")
    valid_bits = int(properties["timestamp_valid_bits"])
    require(1 <= valid_bits <= 64, "Unsupported timestamp queue")
    period = milliseconds(properties["timestamp_period_ns"], "timestamp period")
    require(period > 0, "Zero timestamp period")
    expected = {f"device_observation_{i}" for i in range(51)}
    require({key for key in properties if key.startswith("device_observation_")} == expected,
            "Incomplete/extra original markers")
    legacy = ("transport", "endpoint", "transport", "endpoint", "transport", "endpoint+dense")
    order = [(f"legacy_repeat{i}", stage) for i in range(3) for stage in legacy]
    order += [(f"coupled_repeat{i}", f"{stage}:rows={rows}") for i in range(3)
              for stage, rows in COUPLED_STREAM]
    markers = []
    for i, pair in enumerate(order):
        marker = fields(properties[f"device_observation_{i}"])
        require((marker["phase"], marker["stage"]) == pair and
                marker["available0"] == marker["available1"] == "1", "Changed/unavailable marker stream")
        require(set(marker) == {"phase", "stage", "begin", "end", "available0", "available1", *MARKER_MS},
                "Changed marker fields")
        begin, end = int(marker["begin"]), int(marker["end"])
        require(0 <= begin < 1 << 64 and 0 <= end < 1 << 64, "Invalid timestamp word")
        marker.update({key: milliseconds(marker[key], key) for key in MARKER_MS})
        markers.append(marker)
    values = {key: [milliseconds(properties[f"coupled_repeat{i}_{key}"], key) for i in range(3)]
              for key in ("attempt_ms", "sample_ms", "complete_call_ms")}
    for i in range(3):
        require(values["complete_call_ms"][i] == values["attempt_ms"][i] + values["sample_ms"][i],
                "Complete call does not equal original attempt plus sample")
    phases = {}
    for marker in markers:
        phase = phases.setdefault(marker["phase"], {"stages": {}, "dispatches": []})
        stage = marker["stage"].split(":", 1)[0]
        subtotal = phase["stages"].setdefault(stage, {"submissions": 0, **{key: [] for key in COST_MS}})
        subtotal["submissions"] += 1
        for key in COST_MS:
            subtotal[key].append(marker[key])
        phase["dispatches"].append(marker)
    for phase in phases.values():
        for stage in phase["stages"].values():
            for key in COST_MS:
                stage[key] = total(stage[key])
        phase["totals"] = {key: total(m[key] for m in phase["dispatches"])
                           for key in COST_MS}
    coupled = {key: total(phases[f"coupled_repeat{i}"]["totals"][key] for i in range(3))
               for key in COST_MS}
    coupled_stages = {
        stage: {"submissions": sum(phases[f"coupled_repeat{i}"]["stages"][stage]["submissions"] for i in range(3)),
                **{key: total(phases[f"coupled_repeat{i}"]["stages"][stage][key] for i in range(3))
                   for key in COST_MS}}
        for stage in COUPLED_STAGES}
    identity_fields = ("selected_device", "selected_driver", "retained_scalar_route", "product_mode",
                       "timestamp_valid_bits", "timestamp_period_ns", "coupled_rows", "coupled_projection_budget",
                       "timestamp_scope", "coupled_call_timing_scope")
    selected = {key: properties[key] for key in identity_fields}
    require(all(isinstance(value, str) and value for value in selected.values()), "Missing device/timing identity")
    return {"identity": selected, "values": values, "sums": {key: total(v) for key, v in values.items()},
            "all_phase_costs": phases, "modern_coupled_marker_sums": coupled,
            "modern_coupled_stage_sums": coupled_stages}


def change(a, b):
    difference = b - a
    require(math.isfinite(difference), "Nonfinite difference")
    percent = None if a == 0 else 100 * (b / a - 1)
    require(percent is None or math.isfinite(percent), "Nonfinite percent change")
    return {"baseline_ms": a, "candidate_ms": b, "difference_ms": difference,
            "candidate_percent_change": percent}


def whole_seal(path):
    with path.open('rb') as stream:
        result = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': path.stat().st_size, 'sha256': result}


def validate_execution(action):
    """Join a completed command and independent retirement to the current frozen plan."""
    plan = document(WORKSPACE/'execution-plan.json')
    case = plan['cases'][action]
    producer = plan['producers'][case['producer']]
    folder = WORKSPACE/action
    bridge_path = folder/'native-control-bridge.json'
    bridge = document(bridge_path)
    require(bridge['accepted'] is True and bridge['returncode'] == 0 and not bridge['errors'] and
            not bridge['outer_stop'] and not bridge['controller_signals'] and
            bridge['linux_child_birth_absent'] and bridge['linux_child_group_absent'], 'Incomplete bridge: '+action)
    require(bridge['deadline_seconds'] == 400 and bridge['action'] == action and
            bridge['source_revision'] == plan['source_revision'] and
            bridge['producer'] == case['producer'] and
            bridge['producer_source_revision'] == producer['source_revision'], 'Wrong bridge identity: '+action)
    for name, expected in plan['input_seals'].items():
        require(bridge['input_seals'][name] == expected, 'Plan/bridge seal differs: '+name)
    require(bridge['input_seals'][str((WORKSPACE/'execution-plan.json').relative_to(ROOT))] ==
            whole_seal(WORKSPACE/'execution-plan.json'), 'Plan changed after launch')
    paths = set()
    for records in (bridge['input_seals'], bridge['outputs']):
        for name, expected in records.items():
            path = Path(name) if Path(name).is_absolute() else ROOT/name
            require(whole_seal(path) == expected, 'Current bridge seal differs: '+name)
            paths.add(path)
    receipt = owner(folder/'run', action.split('-',1)[-1])
    require(bridge['Windows_owner'] == receipt, 'Owner/bridge content differs')
    require(receipt['source_revision'] == producer['source_revision'] and
            receipt['live_source_revision'] == plan['source_revision'] and
            receipt['source_build_gate_sha256'] == producer['gate_sha256'] and
            receipt['outer_guard_seconds'] == 90 and receipt['rss_guard_bytes'] == 4294967296 and
            receipt['peak_sampled_rss_bytes'] <= 4294967296, 'Wrong producer/guard: '+action)
    gate = document(ROOT/producer['build_root']/'generated/sirius/native_build_gate.json')
    require(receipt['test_executable_sha256'] == gate['tested_artifacts']['sirius_backend_tests']['sha256'],
            'Wrong sealed native executable')
    require(case['filter'] == CASES[action.split('-',1)[-1]], 'Changed original filter')
    actual_provider = provider(folder/'run')
    for module in actual_provider:
        expected = next(p for p in plan['providers'] if p['name'].casefold()==module['name'])
        require(all(module[key] == (expected[key].casefold() if key=='path' else expected[key])
                    for key in ('path','bytes','sha256')), 'Actual provider differs from sealed provider')
    retirement = document(folder/'retirement-run.json')
    require(retirement['accepted'] and retirement['inputs_unchanged'] and not retirement['errors'] and
            not retirement['signals'] and retirement['returncode']==0 and retirement['birth_absent'] and
            retirement['group_absent'], 'Retirement bridge incomplete: '+action)
    for name, expected in retirement['input_seals'].items():
        require(whole_seal(Path(name))==expected, 'Retirement input changed: '+name)
        paths.add(Path(name))
    census = document(folder/'native-retirement.json')
    require(census['pass'] and not census['errors'] and not census['matches'] and
            len(census['targets'])==2 and all(t['recorded_birth_absent'] for t in census['targets']),
            'Scoped native retirement incomplete: '+action)
    targets = {(t['kind'],t['pid'],t['birth']) for t in census['targets']}
    require(targets == {('owner',receipt['owner_pid'],receipt['owner_start_time_utc']),
                        ('child',receipt['test_pid'],receipt['test_start_time_utc'])}, 'Retirement targeted wrong births')
    paths.update(folder/name for name in ('native-control-bridge.json','retirement-run.json','native-retirement.json','retirement-bootstrap.json'))
    paths.update(folder/'run'/name for name in ('owner.json','driver-modules.json','gtest.xml','stdout.log','stderr.log'))
    return {str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path):whole_seal(path)
            for path in sorted(paths)}


def legacy_readbacks(properties):
    """Check all seven full original interval buffers for each of three repeats."""
    stream=(('transport',3),('endpoint',6),('transport',3),('endpoint',6),
            ('transport',3),('endpoint',6),('dense',3))
    expected={f'legacy_repeat{repeat}_readback_{index}' for repeat in range(3) for index in range(7)}
    actual={key for key in properties if re.match(r'legacy_repeat\d+_readback_',key)}
    require(actual==expected,'Incomplete/extra original legacy readbacks')
    result={}
    for key in sorted(actual):
        value=fields(properties[key]);stage,rows=stream[int(key.rsplit('_',1)[1])]
        require(set(value)=={'stage','bytes','sha256'} and value['stage']==stage,'Wrong legacy buffer stage')
        words={**ROW_WORDS,'dense':4279}[stage]
        require(int(value['bytes'])==rows*words*4,'Partial legacy buffer: '+key)
        digest(value['sha256']);result[key]=value
    return result


def selected_timestamp_modules(properties, hashes, arrays):
    expected={f'selected_module_{index}' for index in range(7)}
    require({key for key in properties if key.startswith('selected_module_')}==expected,'Incomplete/extra selected module records')
    result={}
    names=('Camera','Transport','Endpoint','Dense','Initialize','RayCamera','DopriPhase')
    for index,(stage,name) in enumerate(zip(STAGES,names)):
        value=fields(properties[f'selected_module_{index}']);record=arrays['k'+name+'Shader']
        require(set(value)=={'stage','bytes','sha256'} and value['stage']==stage,'Wrong selected module order')
        require(int(value['bytes'])==record['bytes'] and value['sha256']==record['sha256']==hashes[stage],
                'Selected timestamp module is not exact sealed original producer: '+stage)
        result[stage]=value
    return result


def main():
    require(Path.cwd().resolve() == ROOT, "Invoke from the canonical workspace root")
    output = WORKSPACE / "matched-native-comparison.json"
    require(not output.exists(), "Comparison output already exists")
    records, inputs, reference_reads, selected_provider, selected_device = {}, {}, {}, None, None
    producer_groups, live_revision, guards = {}, None, {}
    reference_legacy = None
    baseline_arrays=document(WORKSPACE/"baseline-native-binding.json")["arrays"]
    candidate_arrays=document(WORKSPACE/"product-inspection.json")["arrays"]
    owner_births, test_births = set(), set()
    validated_seals = {}
    for label in LABELS:
        record = {"owners": {}}
        for case in ("factored", "timestamps"):
            folder = WORKSPACE / (label+"-"+case) / "run"
            for path,expected in validate_execution(label+"-"+case).items():
                if path in validated_seals:require(validated_seals[path]==expected, "Input seal changed during comparison: "+path)
                validated_seals[path]=expected
            receipt, dlls = owner(folder, case), provider(folder)
            owner_birth = (receipt["owner_pid"], receipt["owner_start_time_utc"])
            test_birth = (receipt["test_pid"], receipt["test_start_time_utc"])
            require(owner_birth not in owner_births and test_birth not in test_births,
                    "A matched slot reuses an owner/process receipt")
            owner_births.add(owner_birth)
            test_births.add(test_birth)
            properties = original_properties(folder, case)
            buffers = readbacks(properties, case)
            if case not in reference_reads:
                reference_reads[case] = buffers
            require(reference_reads[case] == buffers, f"Complete {case} buffers differ: {label}")
            if selected_provider is None:
                selected_provider = dlls
            require(dlls == selected_provider, f"Loaded provider differs: {label}/{case}")
            if live_revision is None:
                live_revision = receipt["live_source_revision"]
            require(receipt["live_source_revision"] == live_revision, "Live candidate revision changed")
            group = "baseline" if label.startswith("a") else "candidate"
            producer = {key: receipt[key] for key in
                        ("source_revision", "source_build_gate_sha256", "test_executable_sha256")}
            if group not in producer_groups:
                producer_groups[group] = producer
            require(producer == producer_groups[group], f"Producer changed: {label}/{case}")
            if group == "baseline":
                require(producer == {"source_revision": BASELINE_REVISION, "source_build_gate_sha256": BASELINE_GATE,
                                     "test_executable_sha256": BASELINE_BACKEND}, "Wrong preserved b81 baseline")
            else:
                require(producer["source_revision"] == live_revision and live_revision != BASELINE_REVISION,
                        "Candidate producer is not the live trial revision")
            limits = (receipt["outer_guard_seconds"], receipt["rss_guard_bytes"])
            if case not in guards:
                guards[case] = limits
            require(limits == guards[case], "Matched owner limits changed within an original case")
            record["owners"][case] = {key: receipt[key] for key in
                                      ("source_revision", "live_source_revision", "source_build_gate_sha256",
                                       "test_executable_sha256", "owner_pid", "owner_start_time_utc", "test_pid",
                                       "test_start_time_utc", "elapsed_seconds", "peak_sampled_rss_bytes")}
            for name in ("owner.json", "driver-modules.json", "gtest.xml", "stdout.log", "stderr.log"):
                inputs[f"{label}/{case}/{name}"] = identity(folder / name)
            if case == "factored":
                record[case] = factored_timings(folder)
            else:
                record[case] = timestamps(properties)
                legacy = legacy_readbacks(properties)
                if reference_legacy is None: reference_legacy = legacy
                require(legacy == reference_legacy, "Complete legacy buffers differ: "+label)
                record["legacy_readbacks"] = legacy
                record["timestamp_selected_modules"] = selected_timestamp_modules(properties,record["factored"][0]["selected_stage_sha256"],baseline_arrays if label.startswith("a") else candidate_arrays)
                if selected_device is None:
                    selected_device = record[case]["identity"]
                require(selected_device == record[case]["identity"], "Device/driver/scalar/timing scope changed")
        records[label] = record
    for label in LABELS:
        for route in range(2):
            expected_label = "a1" if label.startswith("a") else "b1"
            require(records[label]["factored"][route]["selected_stage_sha256"] ==
                    records[expected_label]["factored"][route]["selected_stage_sha256"], "Selected module changed within producer")
            require(records[label]["factored"][route]["fma_enabled"] == records["a1"]["factored"][route]["fma_enabled"],
                    "Native FMA capability changed")
    for stage in STAGES:
        a, b = records["a1"]["factored"], records["b1"]["factored"]
        require(a[1]["selected_stage_sha256"][stage] == b[1]["selected_stage_sha256"][stage],
                "Excluded wide selected module changed")
        require((a[0]["selected_stage_sha256"][stage] != b[0]["selected_stage_sha256"][stage]) == (stage in ("transport", "endpoint")),
                "Default selected module change exceeds Transport/Endpoint scope")
    directions = {}
    for name, a, b in (("forward", "a1", "b1"), ("reverse", "a2", "b2")):
        left, right = records[a], records[b]
        directions[name] = {
            "baseline_label": a, "candidate_label": b,
            "factored_routes": [{"route": left["factored"][i]["route"],
                                 **{key: change(left["factored"][i][key], right["factored"][i][key])
                                    for key in ("warm13_wait_ms", "warm13_total_ms")}} for i in range(2)],
            "coupled_complete_calls": {key: change(left["timestamps"]["sums"][key], right["timestamps"]["sums"][key])
                                       for key in ("attempt_ms", "sample_ms", "complete_call_ms")},
            "coupled_per_repeat_calls": [{key: change(left["timestamps"]["values"][key][i], right["timestamps"]["values"][key][i])
                                          for key in ("attempt_ms", "sample_ms", "complete_call_ms")} for i in range(3)],
            "coupled_marker_sums": {key: change(left["timestamps"]["modern_coupled_marker_sums"][key],
                                                 right["timestamps"]["modern_coupled_marker_sums"][key])
                                    for key in COST_MS},
            "coupled_per_repeat_markers": {
                f"coupled_repeat{i}": {
                    key: change(left["timestamps"]["all_phase_costs"][f"coupled_repeat{i}"]["totals"][key],
                                right["timestamps"]["all_phase_costs"][f"coupled_repeat{i}"]["totals"][key])
                    for key in COST_MS} for i in range(3)},
            "coupled_stage_sums": {
                stage: {key: change(left["timestamps"]["modern_coupled_stage_sums"][stage][key],
                                    right["timestamps"]["modern_coupled_stage_sums"][stage][key])
                        for key in COST_MS} for stage in COUPLED_STAGES},
            "coupled_per_repeat_stages": {
                f"coupled_repeat{i}": {
                    stage: {key: change(left["timestamps"]["all_phase_costs"][f"coupled_repeat{i}"]["stages"][stage][key],
                                        right["timestamps"]["all_phase_costs"][f"coupled_repeat{i}"]["stages"][stage][key])
                            for key in COST_MS}
                    for stage in COUPLED_STAGES} for i in range(3)},
        }
    report = {"kind": "original_native_endpoint_and_complete_coupled_abba", "sequence": list(LABELS),
              "observer": identity(Path(__file__).resolve()), "input_files": inputs,
              "producer_groups": producer_groups, "live_source_revision": live_revision,
              "selected_provider_modules": selected_provider, "selected_device": selected_device,
              "owner_limits": {case: {"seconds": limits[0], "rss_bytes": limits[1]}
                               for case, limits in guards.items()},
              "records": records, "directions": directions, "factored_readbacks": reference_reads["factored"],
              "coupled_readbacks": reference_reads["timestamps"], "all_28_factored_complete_readbacks_equal": True,
              "all_11_coupled_complete_readbacks_equal": True,
              "all_21_legacy_complete_readbacks_equal": True, "legacy_readbacks": reference_legacy,
              "all_seven_timestamp_selected_modules_joined_to_actual_producer": True,
              "all_51_original_marker_order_and_availability_equal": True,
              "all_eight_original_owners_completed": True, "retention_decision": "not evaluated",
              "scope": selected_device["coupled_call_timing_scope"], "cold_cache_state_unverified": True,
              "timing_labels": {
                  "warm13_wait_ms": "Factored legacy submit_wait_ms: submit plus wait; excludes the first call",
                  "host_submit_ms": "Marker host submission cost",
                  "host_wait_ms": "Marker host wait cost only",
                  "host_completion_ms": "Marker host submission plus wait cost",
                  "total_ms": "Inclusive synchronous Dispatch cost, not the complete Endpoint API call",
                  "complete_call_ms": "Original attempt_ms plus sample_ms; scope recorded by the test"},
              "bridge_pre_post_source_and_product_seals_checked_here": True,
              "limits": "Recorded loaded DLL identity and completed owner/XML proofs only; this parser verifies bridge source/product seals; root separately joins all eight original controls. Factored total is inclusive Dispatch, not the complete Endpoint API call. Markers include original synchronization/scheduling; no exclusive arithmetic, frame or release qualification.",
              "full_qualification_claimed": False}
    for path,expected in validated_seals.items():
        require(whole_seal(Path(path) if Path(path).is_absolute() else ROOT/path)==expected, "Input changed before comparison completion: "+path)
    report['complete_input_seals']=validated_seals
    payload = json.dumps(report, indent=2, allow_nan=False) + "\n"
    require(len(payload.encode()) <= 4 * 1024 * 1024, "Comparison report exceeds bound")
    with output.open("x", encoding="utf-8") as stream:
        stream.write(payload)
    print(json.dumps({"output": str(output.relative_to(ROOT)), "directions": directions,
                      "retention_decision": "not evaluated"}, allow_nan=False))


if __name__ == "__main__":
    try:
        main()
    except (KeyError, TypeError, ValueError, OSError, OverflowError, ET.ParseError) as error:
        print(f"Matched native comparison invalid: {error}", file=sys.stderr)
        sys.exit(1)
