"""Inspect actual ordinary IR-to-called-JIT joins; no arithmetic/cost verdict.

ELF reader and defined-relocation/uncovered-byte checks derive from the sealed
b5d2abf observer. Undefined callback addresses are checked for consistency;
this diagnostic does not independently resolve them.
"""
from pathlib import Path
import collections
import hashlib
import json
import re
import struct
import sys
sys.dont_write_bytecode = True

def elf_sections(data):
    assert data[:6] == b'\x7fELF\x02\x01'
    assert struct.unpack_from('<HH', data, 16) == (1, 62)
    offset = struct.unpack_from('<Q', data, 40)[0]
    size, count, names_index = struct.unpack_from('<HHH', data, 58)
    assert size == 64 and count < 100
    sections = [struct.unpack_from('<IIQQQQIIQQ', data, offset + i * size)
                for i in range(count)]
    strings = sections[names_index]
    names = data[strings[4]:strings[4] + strings[5]]
    named = []
    for i, row in enumerate(sections):
        name = names[row[0]:].split(b'\0', 1)[0].decode()
        named.append({'index': i, 'name': name, 'type': row[1],
                      'flags': row[2], 'address': row[3], 'offset': row[4],
                      'bytes': row[5], 'link': row[6], 'info': row[7],
                      'entry_bytes': row[9]})
    symbols = []
    for table in named:
        if table['type'] != 2:
            continue
        strings = named[table['link']]
        names = data[strings['offset']:strings['offset'] + strings['bytes']]
        assert table['entry_bytes'] == 24
        for start in range(table['offset'], table['offset'] + table['bytes'], 24):
            name, info, other, section, value, size = struct.unpack_from('<IBBHQQ', data, start)
            if info & 15 == 2 and section < len(named) and size:
                symbols.append({'name': names[name:].split(b'\0', 1)[0].decode(),
                                'section': section, 'offset': value, 'bytes': size,
                                'address': named[section]['address'] + value})
    return named, symbols

def linked_review(folder, event, variant):
    data = (folder / event['file']).read_bytes()
    assert hashlib.sha256(data).hexdigest() == event['sha256']
    sections, functions = elf_sections(data)
    assert sections == event['sections'] and functions == event['functions']
    live = {s['index']: s for s in event['live_sections_after_pointer_lookup']}
    assert set(live) == {s['index'] for s in sections if s['flags'] & 2 and s['address'] and s['bytes']}
    decoded_symbols = {}
    for table in sections:
        if table['type'] != 2:
            continue
        names = sections[table['link']]
        strings = data[names['offset']:names['offset'] + names['bytes']]
        syms = []
        for start in range(table['offset'], table['offset'] + table['bytes'], 24):
            name, info, other, section, value, size = struct.unpack_from('<IBBHQQ', data, start)
            syms.append({'name': strings[name:].split(b'\0', 1)[0].decode(), 'section': section, 'value': value})
        decoded_symbols[table['index']] = syms
    external = collections.defaultdict(set)
    relocation_counts = collections.Counter()
    reviewed = []
    for index, capture in live.items():
        section = sections[index]
        live_data = (folder / capture['path']).read_bytes()
        assert len(live_data) == section['bytes']
        assert hashlib.sha256(live_data).hexdigest() == capture['sha256']
        unlinked = data[section['offset']:section['offset'] + section['bytes']]
        expected = bytearray(unlinked)
        mask = bytearray(len(expected))
        for table in sections:
            if table['type'] != 4 or table['info'] != index:
                continue
            assert table['entry_bytes'] == 24
            for start in range(table['offset'], table['offset'] + table['bytes'], 24):
                offset, info, addend = struct.unpack_from('<QQq', data, start)
                kind, symbol = info & 0xffffffff, decoded_symbols[table['link']][info >> 32]
                assert kind in (1, 24) and offset + 8 <= len(expected)
                relocation_counts[kind] += 1
                observed = struct.unpack_from('<Q', live_data, offset)[0]
                place = section['address'] + offset
                if symbol['section'] == 0:
                    # Undefined callback addresses are checked for consistency,
                    # not independently resolved from a public exported name.
                    target = (observed - addend + (place if kind == 24 else 0)) & ((1 << 64) - 1)
                    external[symbol['name']].add(target)
                else:
                    target = (sections[symbol['section']]['address'] + symbol['value']
                              if symbol['section'] < len(sections) else symbol['value'])
                relocated = (target + addend - (place if kind == 24 else 0)) & ((1 << 64) - 1)
                struct.pack_into('<Q', expected, offset, relocated)
                mask[offset:offset + 8] = b'\1' * 8
        assert bytes(expected) == live_data, 'unexpected non-relocation or defined-target byte difference'
        reviewed.append({'name': section['name'], 'bytes': len(expected),
                         'changed_bytes': sum(a != b for a, b in zip(unlinked, live_data)),
                         'relocation_field_bytes': sum(mask),
                         'all_defined_relocations_and_uncovered_bytes_exact': True})
    assert all(len(values) == 1 for values in external.values())
    return {'file': event['file'], 'variant': variant,
            'live_sections': reviewed, 'relocation_counts': dict(relocation_counts),
            'undefined_callback_targets_consistent': {k: list(v)[0] for k, v in external.items()},
            'undefined_callback_targets_independently_resolved': False}


def inspect(folder):
    owner = json.loads((folder / 'owner.json').read_text())
    assert owner['passed'] and owner['child_exit'] == 0 and not owner['stop_reason']
    assert json.loads((folder / 'inputs-before.json').read_text()) == json.loads(
        (folder / 'inputs-after.json').read_text())
    limited = json.loads((folder / 'capture-verification.json').read_text())
    assert limited['passed']
    bindings = json.loads((folder / 'bindings.json').read_text())
    dispatches = bindings['dispatches']
    assert len(dispatches) == 1
    dispatch = dispatches[0]
    assert dispatch['id'] == 0 and dispatch['completed']
    assert dispatch['kernel'] == 0 and dispatch['groups_x'] == 1
    calls = [c for c in bindings['entry_calls'] if c['dispatch_id'] == 0]
    assert len(calls) == 1
    call = calls[0]
    assert call['buffers'] == dispatch['buffers'] and call['kernel'] == dispatch['kernel']
    events = [e for e in bindings['jit_events'] if e.get('file') == call['file']]
    assert len(events) == 1
    event = events[0]
    assert call['function']['name'] == 'cs_variant'
    assert call['function'] in event['functions']
    assert call['actual_entry_pc'] == call['function']['address'] == event['actual_callable_pointer']
    lookups = [item for item in bindings['pointer_lookups'] if call['file'] in item['matched_files']]
    assert lookups and all(item['pointer'] == call['actual_entry_pc'] for item in lookups)
    generations = {item['module_generation'] for item in lookups}
    assert len(generations) == 1
    generation = generations.pop()
    created = bindings['module_creations'][generation]
    assert created['id'] == generation
    assert all(item['module'] == created['module'] for item in lookups)
    exports = [e for e in bindings['ir_exports'] if e['module_generation'] == generation]
    assert len(exports) == 1 and exports[0] == limited['ordinary_ir_export']
    export = exports[0]
    assert export['module'] == created['module']
    raw = (folder / export['path']).read_bytes()
    assert len(raw) == export['bytes'] and hashlib.sha256(raw).hexdigest() == export['sha256']
    pipelines = [bindings['optimizer_passes'][i] for i in export['successful_optimizer_pass_ids']]
    assert all(p['completed'] and p['error_pointer'] == 0 and p['module_generation'] == generation for p in pipelines)
    assert all(p['id'] == i and p['module'] == created['module']
        for i, p in zip(export['successful_optimizer_pass_ids'], pipelines))
    assert [p['pipeline'] for p in pipelines] == ['default<O0>',
        'sroa,early-cse,simplifycfg,reassociate,mem2reg,instsimplify,instcombine<no-verify-fixpoint>']
    text = raw.decode()
    functions = re.findall(r'^define[^\n]* @([^ (]+)\(', text, re.M)
    assert 'cs_variant' in functions
    name = folder.name
    producer, variant = name.removeprefix('capture-').split('-', 1)
    assert producer in {'stock', 'candidate'}
    constants = re.findall(r'@(?:coro_malloc|malloc|aligned_alloc)\([^\n]*', text)
    frame_sizes = [int(n) for n in re.findall(r'@coro_malloc\(i32 (\d+)\)', text)]
    resume_frames = [int(n) for n in re.findall(
        r'^define[^\n]*@cs_co_variant\.resume\([^\n]*dereferenceable\((\d+)\)', text, re.M)]
    linked = linked_review(folder, event, variant)
    executable_sections = [s for s in event['sections'] if s['flags'] & 6 == 6]
    return {'case': name, 'producer': producer, 'variant': variant,
        'actual_loaded_module': limited['actual_loaded_module'],
        'single_completed_dispatch': dispatch, 'actual_entry': call,
        'ordinary_ir_export': export, 'ordinary_pipelines': [p['pipeline'] for p in pipelines],
        'linked_section_review': linked,
        'jit_object_bytes': event['object_bytes'],
        'jit_defined_functions': event['functions'],
        'jit_executable_section_bytes': sum(s['bytes'] for s in executable_sections),
        'llvm_defined_functions': functions,
        'llvm_float_multiply_instruction_count': len(re.findall(r' = fmul ', text)),
        'llvm_alloca_instruction_count': len(re.findall(r' = alloca ', text)),
        'llvm_allocator_call_text': constants,
        'llvm_coroutine_allocator_constant_bytes': frame_sizes,
        'llvm_resume_frame_dereferenceable_bytes': resume_frames,
        'scope': 'Actual loaded complete module -> ordinary successful LLVM passes -> lookup generation -> called JIT entry and allocated-section byte/defined-relocation joins. Short-input canary only; no retained arithmetic or runtime benefit verdict. Undefined callbacks are consistency checked only.'}


if __name__ == '__main__':
    work = Path(__file__).resolve().parent
    names = sys.argv[1:]
    assert names and all(n.startswith('capture-') and Path(n).name == n for n in names)
    reports = [inspect(work / name) for name in names]
    for report in reports:
        target = work / report['case'] / 'actual-code-review.json'
        assert not target.exists()
        target.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps([{'case': r['case'], 'entry_function_bytes': r['actual_entry']['function']['bytes'],
        'executable_section_bytes': r['jit_executable_section_bytes'],
        'defined_functions': [(f['name'], f['bytes']) for f in r['jit_defined_functions']]} for r in reports], indent=2))
