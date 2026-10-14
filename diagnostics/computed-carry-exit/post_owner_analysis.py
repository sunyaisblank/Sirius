"""Analyze only the exact issue-101 postprocessing-order failure, without repair.

Original owner, capture files and frozen helpers remain read-only. The original
analyzer runs in memory with its owner precondition replaced by the stricter
incident-specific terminal-failure check below; every binding/byte check remains.
"""
from pathlib import Path
import hashlib
import json
import os
import resource
import sys

sys.dont_write_bytecode = True
ROOT = Path.cwd().resolve()
WORK = Path(__file__).resolve().parent
CAPTURE = WORK / 'capture-candidate-transport-normal'
OUTPUT = WORK / 'post-owner-analysis'
ANALYZER = WORK / 'analyze_capture.py'
PINS = {
    ANALYZER: '58d484832210971e160c836a98f33e22bb9767d1e2c7a6883a1332c6acecf708',
    WORK / 'capture_review.py': 'f34f19c89ea07faf005e11c1001038f2b7bcab7aa84e3b0faab09fcbfe1f8cf1',
    CAPTURE / 'owner.json': 'd13d06580e22bd899f65b31a707997d2c7bb7198acea1dacbc275b7f8428462a',
    CAPTURE / 'stderr.log': '6b326ebd1dd44a54198596e4309966737e9e0406333a220bc70181884c7885dc',
    CAPTURE / 'capture-verification.json': 'c42083d19638220d9b0cc90913bda8d569c29102be3b88389e551b05c1335b1a',
}


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def identity(path):
    path = Path(path)
    return {'resolved_path': str(path.resolve(strict=True)),
            'bytes': path.stat().st_size, 'sha256': sha(path)}


def exact_failed_owner(owner, folder):
    assert folder == CAPTURE
    assert sha(folder / 'owner.json') == PINS[folder / 'owner.json']
    assert owner['passed'] is False and owner['child_exit'] == 1
    assert owner['status'] == 'terminal' and owner['stop_reason'] is None
    assert owner['cleanup_errors'] == [] and owner['remaining_owned_processes'] == []
    assert owner['source_unchanged'] is True and owner['whole_inputs_unchanged'] is True
    assert 'observer_error' not in owner
    assert owner['observed_births'] and all(x['absent'] is True for x in owner['observed_births'])
    assert owner['provider_mappings']
    assert owner['limits'] == {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152}
    assert owner['wall_seconds'] < 300 and owner['peak_sampled_rss_kib'] <= 4194304
    assert owner['minimum_sampled_available_kib'] >= 2097152
    assert owner['source'] == {'revision': '6e13108d5e77b053ff2b0b391e42984e243733c7', 'status': ''}
    actions = json.loads((WORK / 'actions.json').read_text())
    expected = [x.replace('{folder}', str(CAPTURE)) for x in
                actions['capture-candidate-transport-normal']['argv']]
    assert owner['argv'] == expected
    assert json.loads((folder / 'inferior-exit.json').read_text()) == {'exit_code': 0}
    assert json.loads((folder / 'capture-verification.json').read_text())['passed'] is True
    assert not (folder / 'actual-code-review.json').exists()
    stderr = (folder / 'stderr.log').read_text()
    assert 'analyze_capture.py", line 113, in inspect' in stderr
    assert "assert owner['passed'] and owner['child_exit'] == 0 and not owner['stop_reason']" in stderr
    assert 'capture_review.py", line 12' in stderr
    assert stderr.count('AssertionError') == 1 and stderr.count('CalledProcessError') == 2


def main():
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    resource.setrlimit(resource.RLIMIT_AS, (4 * 1024**3, 4 * 1024**3))
    resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024**2, 64 * 1024**2))
    assert not OUTPUT.exists()
    for path, expected in PINS.items():
        assert sha(path) == expected, path
    before_raw = {str(p.relative_to(CAPTURE)): identity(p)
                  for p in CAPTURE.rglob('*') if p.is_file()}
    before = json.loads((CAPTURE / 'inputs-before.json').read_text())
    assert before == json.loads((CAPTURE / 'inputs-after.json').read_text())
    for name, expected in before.items():
        path = Path(name) if Path(name).is_absolute() else ROOT / name
        assert identity(path) == expected, name
    owner = json.loads((CAPTURE / 'owner.json').read_text())
    exact_failed_owner(owner, CAPTURE)
    inputs = json.loads((CAPTURE / 'capture-inputs.json').read_text())
    loaded = json.loads((CAPTURE / 'actual-loaded-module.json').read_text())
    verified = json.loads((CAPTURE / 'capture-verification.json').read_text())
    assert loaded == verified['actual_loaded_module']
    assert loaded['actual_span_matches_frozen_input'] is True
    assert sha(inputs['module']) == loaded['sha256'] == inputs['module_sha256']
    assert loaded['sha256'] == '0d8e5784035d2870c66420dd91607d4f5c35ed4e1360c402501ea1792a0708fb'
    original = ANALYZER.read_text()
    needle = "    assert owner['passed'] and owner['child_exit'] == 0 and not owner['stop_reason']"
    assert original.count(needle) == 1
    adapted = original.replace(needle, '    exact_failed_owner(owner, folder)', 1)
    namespace = {'__name__': 'issue101_exact_failure_analyzer', '__file__': str(ANALYZER),
                 'exact_failed_owner': exact_failed_owner}
    exec(compile(adapted, str(ANALYZER) + ':exact-terminal-failure-adapter', 'exec'), namespace)
    review = namespace['inspect'](CAPTURE)
    after_raw = {str(p.relative_to(CAPTURE)): identity(p)
                 for p in CAPTURE.rglob('*') if p.is_file()}
    assert after_raw == before_raw
    assert sha(CAPTURE / 'owner.json') == PINS[CAPTURE / 'owner.json']
    OUTPUT.mkdir()
    (OUTPUT / 'actual-code-review.json').write_text(json.dumps(review, indent=2) + '\n')
    receipt = {'analysis_passed': True, 'original_owner_passed': owner['passed'],
        'original_child_exit': owner['child_exit'], 'original_owner_preserved': True,
        'original_analyzer_sha256': PINS[ANALYZER],
        'original_capture_wrapper_sha256': PINS[WORK / 'capture_review.py'],
        'adapter_sha256': sha(Path(__file__)),
        'in_memory_analyzer_sha256': hashlib.sha256(adapted.encode()).hexdigest(),
        'only_analyzer_change': {'original_statement': needle.strip(),
                                 'replacement': 'exact_failed_owner(owner, folder)'},
        'phase_ordering_diagnosis': 'capture_review.py invoked analyze_capture.py before its enclosing run_owned.py owner became terminal/passed. The canary and capture verification succeeded first. The analyzer then failed only at its initial owner precondition, propagating exit1. This analysis accepts only that exact immutable terminal failure; it does not promote or rewrite the failed owner.',
        'unchanged_raw_files': before_raw, 'exact_original_input_records_rehashed': len(before),
        'postprocessing_process_pid': os.getpid(),
        'scope': 'Post-owner code-binding/relocation analysis only. Original command remains failed; no GPU rerun, retained arithmetic, speed or qualification verdict.'}
    (OUTPUT / 'analysis-receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'analysis_passed': True, 'original_owner_passed': False,
        'original_child_exit': 1, 'input_records_rehashed': len(before),
        'raw_files_unchanged': len(before_raw), 'jit_executable_section_bytes': review['jit_executable_section_bytes'],
        'coroutine_frames': review['llvm_resume_frame_dereferenceable_bytes'],
        'output': str(OUTPUT)}, indent=2))


if __name__ == '__main__':
    main()
