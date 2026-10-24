"""Reproduce/check only the three private-owner context statement moves.

Default is read-only. --write may recreate missing derived headers and receipts;
existing bytes must already equal this derivation and are never overwritten.
"""
from pathlib import Path
import difflib
import hashlib
import json
import re
import sys

root = Path.cwd()
work = Path(__file__).resolve().parent
basis = root / 'bin/linux-gcc/native-cpu-slang-bridge/typed-batch'
destination = root / 'bin/linux-gcc/native-cpu-owner-context/typed-batch'
assert sys.argv[1:] in ([], ['--write'])
write = sys.argv[1:] == ['--write']
records = {}
patch = ''


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def publish(path, data):
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        assert write, 'Missing derived input: ' + str(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)


for slug in ('transport', 'endpoint', 'dopri_phase'):
    source = basis / (slug + '.native.hpp')
    text = source.read_text()
    start = text.index('void sirius_' + slug + '_ComputeMain(void* ')
    end = text.index('// [numthreads(64, 1, 1)]', start)
    body = text[start:end]
    declaration, = re.findall(r'^    KernelContext_0 kernelContext_\d+\{\};\n', body, re.M)
    assignment, = re.findall(r'^    \(&kernelContext_\d+\)->globalParams_0 = \(slang_bit_cast<GlobalParams_0\*>\(globalParams_1\)\);\n', body, re.M)
    assert re.search(r'kernelContext_\d+', declaration)[0] == re.search(r'kernelContext_\d+', assignment)[0]
    assert body.count('ComputeOwnedRow_0(') == 1 and body.count('== 0U)') == 1
    stripped = body.replace(declaration, '', 1).replace(assignment, '', 1)
    location = stripped.index('        ComputeOwnedRow_0(')
    moved = stripped[:location] + '    ' + declaration + '    ' + assignment + stripped[location:]
    assert moved.replace('    ' + declaration, '', 1).replace('    ' + assignment, '', 1) == stripped
    derived = text[:start] + moved + text[end:]
    assert derived[:start] == text[:start]
    assert derived[start + len(moved):] == text[end:]
    path = destination / (slug + '.native.hpp')
    publish(path, derived.encode())
    records[slug] = {'basis_path': str(source.relative_to(root)), 'basis': pin(source.read_bytes()),
        'derived_path': str(path.relative_to(root)), 'derived': pin(derived.encode()),
        'entrypoint': 'sirius_' + slug + '_ComputeMain',
        'moves': [declaration.strip(), assignment.strip()], 'all_other_header_bytes_equal': True,
        'original_Group_text_equal': derived[end + len(moved) - len(body):] == text[end:]}
    patch += ''.join(difflib.unified_diff(text.splitlines(True), derived.splitlines(True),
        fromfile='a/' + str(source.relative_to(root)), tofile='b/' + str(path.relative_to(root))))
record = {'scope': 'Source preparation only: move only private KernelContext{} declaration and globalParams assignment inside the preexisting lane-zero branch. All64Group calls, integer lanecheck/groupIDreads, ownerbody/arithmetic/tables/order remain unchanged. No compile/build/runtime or gain claim.', 'headers': records}
publish(work / 'owner-context.patch', patch.encode())
publish(work / 'header-preparation.json', (json.dumps(record, indent=2) + '\n').encode())
print(json.dumps({'derived_headers_checked': 3, 'statement_moves': 6,
                  'other_generated_bytes_changed': False, 'compiled_or_executed': False}))
