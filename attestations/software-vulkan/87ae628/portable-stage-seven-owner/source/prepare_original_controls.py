from pathlib import Path
import hashlib, json, shlex, shutil, subprocess

ROOT = Path.cwd()
WORK = ROOT / 'out/portable-stage-seven-owner-construction'
PHASE = WORK / 'original-controls-preparation'
BIN = ROOT / 'bin/linux-gcc/portable-stage-seven-owner'
PHASE.mkdir()

def sha(data):
    return hashlib.sha256(data).hexdigest()

def saved(path, data):
    assert not path.exists(), path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)

def original(relative):
    return subprocess.run(['git', 'show', '9fbd0e8:' + relative], check=True, capture_output=True).stdout

def transform(text, old, new, edits):
    assert text.count(old) == 1, (old, text.count(old))
    edits.append({'old': old, 'new': new})
    return text.replace(old, new)

def region(text, start, end, callback, edits):
    a, b = text.index(start), text.index(end, text.index(start))
    old = text[a:b]
    new = callback(old)
    return transform(text, old, new, edits)

def compute_source(base, candidate):
    text = base.decode()
    edits = []
    owner = 'true' if candidate else 'false'
    text = transform(text, 'namespace {\nusing namespace sirius::backend;',
                     'namespace {\nusing namespace sirius::backend;\n'
                     f'inline constexpr bool kDiagnosticOwner = {owner};', edits)
    def probe(s):
        s = s.replace('    Failure fail_next = Failure::None;',
            '    std::vector<KernelHandle> loaded_handles;\n'
            '    KernelHandle last_kernel{};\n'
            '    std::array<KernelHandle, 2> last_pair{};\n'
            '    std::size_t last_bindings = 0;\n'
            '    std::array<std::size_t, 2> pair_bindings{};\n'
            '    Failure fail_next = Failure::None;')
        s = s.replace('        if (loaded && after_load) after_load(code);',
            '        if (loaded) loaded_handles.push_back(*loaded);\n'
            '        if (loaded && after_load) after_load(code);')
        s = s.replace('        ++pair_calls;',
            '        last_pair = {commands[0].kernel, commands[1].kernel};\n'
            '        pair_bindings = {commands[0].buffers.size(), commands[1].buffers.size()};\n'
            '        ++pair_calls;')
        s = s.replace('        ++dispatch_calls;',
            '        last_kernel = kernel; last_bindings = buffers.size();\n'
            '        ++dispatch_calls;')
        return s
    text = region(text, 'class TransferProbeDevice final', '\nvoid CheckStickyError', probe, edits)
    def timestamps(s):
        s = s.replace('    std::size_t modules = 0;',
            '    std::size_t modules = 0;\n    std::vector<std::string> module_descriptions;')
        s = s.replace('        ASSERT_LT(modules, RetainedCompute::kStageCount);',
            '        ASSERT_LT(modules, RetainedCompute::kStageCount + (kDiagnosticOwner ? 2 : 0));')
        old = '''        RecordProperty("selected_module_" + std::to_string(modules),
                       std::format("stage={};bytes={};sha256={}",
                                   RetainedCompute::StageName(
                                       static_cast<RetainedCompute::KernelStage>(modules)),
                                   code.size_bytes(), *digest));'''
        new = '''        const auto role = modules < RetainedCompute::kStageCount
            ? RetainedCompute::StageName(static_cast<RetainedCompute::KernelStage>(modules))
            : (modules == RetainedCompute::kStageCount ? "transport" : "endpoint");
        const auto description = std::format("stage={};bytes={};sha256={}", role, code.size_bytes(), *digest);
        module_descriptions.push_back(description);
        RecordProperty((modules < RetainedCompute::kStageCount ? "selected_module_" : "optional_module_")
                       + std::to_string(modules), description);'''
        assert s.count(old) == 1
        s = s.replace(old, new)
        call = 'RetainedCompute::Create(probe, 24, fp64_products)'
        if candidate:
            s = s.replace(call, 'RetainedCompute::Create(probe, 24, fp64_products, 250, true)')
            s = s.replace('    ASSERT_TRUE(observed) << observed.error().Description();',
                '    ASSERT_TRUE(observed) << observed.error().Description();\n'
                '    ASSERT_TRUE((*observed)->StageSevenOwnershipEnabled());\n'
                '    RecordProperty("owner_provisioning_ms", std::format("{:.17g}", (*observed)->StageSevenProvisioningMs()));')
        else:
            s = s.replace('    ASSERT_TRUE(observed) << observed.error().Description();',
                '    ASSERT_TRUE(observed) << observed.error().Description();\n'
                '    RecordProperty("owner_provisioning_ms", "0");')
        s = s.replace('    auto observed = ',
            '    const auto create_started = std::chrono::steady_clock::now();\n    auto observed = ', 1)
        s = s.replace('    ASSERT_TRUE(observed)',
            '    const auto create_finished = std::chrono::steady_clock::now();\n'
            '    RecordProperty("observed_create_ms", std::format("{:.17g}", std::chrono::duration<double, std::milli>(create_finished - create_started).count()));\n'
            '    ASSERT_TRUE(observed)', 1)
        s = s.replace('    ASSERT_EQ(modules, RetainedCompute::kStageCount);',
            '    ASSERT_EQ(modules, RetainedCompute::kStageCount + (kDiagnosticOwner ? 2 : 0));')
        s = s.replace('    ASSERT_EQ(probe.allocations.size(), 2 * RetainedCompute::kStageCount);',
            '    ASSERT_EQ(probe.allocations.size(), 2 * RetainedCompute::kStageCount + (kDiagnosticOwner ? 1 : 0));')
        s = s.replace('i < probe.allocations.size()', 'i < 2 * RetainedCompute::kStageCount')
        s = s.replace('    const auto capture =', '''    const auto execution = [&](std::size_t marker, std::size_t command, KernelHandle kernel, std::size_t bindings) {
        const auto found = std::find_if(probe.loaded_handles.begin(), probe.loaded_handles.end(),
            [&](KernelHandle loaded) { return loaded.value == kernel.value; });
        ASSERT_NE(found, probe.loaded_handles.end());
        const auto index = static_cast<std::size_t>(found - probe.loaded_handles.begin());
        ASSERT_LT(index, module_descriptions.size());
        RecordProperty("executed_module_" + std::to_string(marker) + "_" + std::to_string(command),
            module_descriptions[index] + std::format(";loaded_index={};handle={};bindings={}", index, kernel.value, bindings));
    };
    const auto capture =''')
        s = s.replace('host->command_setup_ms, host->cleanup_ms, host->total_ms));',
            'host->command_setup_ms, host->cleanup_ms, host->total_ms));\n'
            '        if (stage == "endpoint+dense") {\n'
            '            execution(samples - 1, 0, probe.last_pair[0], probe.pair_bindings[0]);\n'
            '            execution(samples - 1, 1, probe.last_pair[1], probe.pair_bindings[1]);\n'
            '        } else execution(samples - 1, 0, probe.last_kernel, probe.last_bindings);')
        s = s.replace('    const auto coupled_baseline = AttemptRetainedDopriIntervals(**observed, coupled, 24);',
            '    const auto disabled_attempt_started = Clock::now();\n'
            '    const auto coupled_baseline = AttemptRetainedDopriIntervals(**observed, coupled, 24);\n'
            '    const auto disabled_attempt_finished = Clock::now();\n'
            '    RecordProperty("disabled_attempt_ms", std::format("{:.17g}", milliseconds(disabled_attempt_started, disabled_attempt_finished)));')
        s = s.replace('    const auto sample_baseline =\n',
            '    const auto disabled_sample_started = Clock::now();\n    const auto sample_baseline =\n')
        s = s.replace('    ASSERT_TRUE(sample_baseline)',
            '    const auto disabled_sample_finished = Clock::now();\n'
            '    RecordProperty("disabled_sample_ms", std::format("{:.17g}", milliseconds(disabled_sample_started, disabled_sample_finished)));\n'
            '    ASSERT_TRUE(sample_baseline)', 1)
        return s
    text = region(text, 'void RetainedComputeTest::CheckOriginalIntervalTimestamps(',
                  '\nTEST_F(RetainedComputeTest, DeviceTimestamps', timestamps, edits)
    if candidate:
        def mixed(s):
            s = s.replace('RetainedCompute::Create(probe, 24, wide)', 'RetainedCompute::Create(probe, 24, wide, 250, true)')
            s = s.replace('            ASSERT_EQ(probe.allocations.size(), 2 * RetainedCompute::kStageCount);',
                '            const bool owner = !decline_rte32 && RetainedUsesPortableArithmetic(probe.Info()) && probe.Info().rounds_fp32_to_nearest;\n'
                '            ASSERT_EQ(instance.StageSevenOwnershipEnabled(), owner);\n'
                '            ASSERT_EQ(probe.allocations.size(), 2 * RetainedCompute::kStageCount + (owner ? 1 : 0));')
            s = s.replace('            EXPECT_EQ(allocation, *required);',
                '            const auto extra = RetainedCompute::RequiredStageSevenOwnerAllocationBytes(probe, 24);\n'
                '            ASSERT_TRUE(extra);\n'
                '            EXPECT_EQ(allocation, *required + (owner ? *extra : 0));')
            return s
        text = region(text, 'TEST_F(RetainedComputeTest, MixedIndependentTransportLayersPreserveWordsAndRefusal)',
                      '\nTEST_F(RetainedComputeTest, FactoredEndpoints', mixed, edits)
        def fixture(s):
            s = s.replace('        auto created = RetainedCompute::Create(*device, 24);',
                '        const bool owner = std::string_view(::testing::UnitTest::GetInstance()->current_test_info()->name()) == "CoupledIntervalsRequireEmbeddedAndIndependentDenseAgreement";\n'
                '        auto created = RetainedCompute::Create(*device, 24, false, 250, owner);')
            s = s.replace('        EXPECT_EQ(device->BufferAllocationBytes(), *required);',
                '        const auto extra = RetainedCompute::RequiredStageSevenOwnerAllocationBytes(*device, 24);\n'
                '        ASSERT_TRUE(extra);\n'
                '        EXPECT_EQ(device->BufferAllocationBytes(), *required + (owner ? *extra : 0));')
            return s
        text = region(text, 'class RetainedComputeTest :', '\n#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE\n// Observe', fixture, edits)
    return text.encode(), edits

for candidate, name in [(False, 'baseline'), (True, 'candidate')]:
    base = (ROOT/'tests/backend/retained_compute_test.cpp').read_bytes() if candidate else original('tests/backend/retained_compute_test.cpp')
    saved(PHASE/name/'original_retained_compute_test.cpp', base)
    source, edits = compute_source(base, candidate)
    saved(PHASE/name/'retained_compute_test.cpp', source)
    saved(PHASE/name/'edits.json', (json.dumps({'base_sha256':sha(base), 'source_sha256':sha(source), 'edits':edits},indent=2)+'\n').encode())
saved(PHASE/'baseline/include/sirius/backend/retained_compute.h', original('src/sirius/backend/retained_compute.h'))
dopri = (ROOT/'tests/backend/retained_dopri_test.cpp').read_text()
assert dopri.encode() == original('tests/backend/retained_dopri_test.cpp')
dopri_edits = []
dopri = transform(dopri, 'auto created = RetainedCompute::Create(*device, 24, wide);',
                  'auto created = RetainedCompute::Create(*device, 24, wide, 250, true);', dopri_edits)
dopri = transform(dopri, '            EXPECT_EQ(device->BufferAllocationBytes(), *required);',
    '            const auto extra = RetainedCompute::RequiredStageSevenOwnerAllocationBytes(*device, 24);\n'
    '            ASSERT_TRUE(extra);\n'
    '            ASSERT_TRUE((*created)->StageSevenOwnershipEnabled());\n'
    '            EXPECT_EQ(device->BufferAllocationBytes(), *required + *extra);', dopri_edits)
saved(PHASE/'candidate/retained_dopri_test.cpp', dopri.encode())
saved(PHASE/'candidate/dopri-edits.json', (json.dumps(dopri_edits,indent=2)+'\n').encode())

limits = {'seconds':300, 'rss_kib':4194304, 'reserve_kib':2097152}
commands = json.loads((ROOT/'bin/linux-gcc/compile_commands.json').read_text())
actions = {}
for name, relative, source, output, includes in [
    ('compile-baseline-observer', 'tests/backend/retained_compute_test.cpp', PHASE/'baseline/retained_compute_test.cpp', BIN/'original_baseline_test.o', [PHASE/'baseline/include']),
    ('compile-candidate-observer', 'tests/backend/retained_compute_test.cpp', PHASE/'candidate/retained_compute_test.cpp', BIN/'original_candidate_test.o', [WORK/'full-builder/products']),
    ('compile-candidate-dopri', 'tests/backend/retained_dopri_test.cpp', PHASE/'candidate/retained_dopri_test.cpp', BIN/'original_candidate_dopri.o', [WORK/'full-builder/products']),
    ('compile-candidate-executor', 'src/sirius/backend/retained_trace_executor.cpp', ROOT/'src/sirius/backend/retained_trace_executor.cpp', BIN/'executor.o', [WORK/'full-builder/products']),
]:
    preset = next(c for c in commands if c['file']==str(ROOT/relative))
    argv = shlex.split(preset['command'])
    argv[argv.index('-o')+1] = str(output)
    argv[argv.index('-c')+1] = str(source)
    argv[1:1] = ['-I'+str(i) for i in includes]
    actions[name] = {'revision':subprocess.run(['git','rev-parse','HEAD'],capture_output=True,text=True,check=True).stdout.strip(),
        'kind':'diagnostic','expected_exit':0,'argv':argv,'limits':limits,
        'frozen_inputs':[str(source.relative_to(ROOT)), str((ROOT/'bin/linux-gcc/compile_commands.json').relative_to(ROOT)),
                         str((WORK/'prepare_original_controls.py').relative_to(ROOT)), '/usr/bin/g++-14', '/usr/libexec/gcc/x86_64-linux-gnu/14/cc1plus', '/usr/bin/as', '/usr/bin/ld',
                         'attestations/software-vulkan/dc714b9/guarded-upmultiply-trial/diagnostic/actions.json']
            + [str(p.relative_to(ROOT)) for p in PHASE.rglob('*') if p.is_file()]
            + [str((WORK/'full-builder/products/retained_owner_layout.h').relative_to(ROOT))],
        'scope':'Original-control source preparation with configured C++ preset flags; no Vulkan or qualification execution.'}
saved(PHASE/'actions.json',(json.dumps(actions,indent=2)+'\n').encode())
for file in ['run_owned.py','owned_guard.py']:
    shutil.copyfile(WORK/'runtime-stream-v2'/file, PHASE/file)
print(json.dumps({'phase':str(PHASE.relative_to(ROOT)), 'actions':list(actions)},indent=2))
