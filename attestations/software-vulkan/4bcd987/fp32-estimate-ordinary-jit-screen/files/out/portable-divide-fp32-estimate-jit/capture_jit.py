import gdb
import hashlib
import json
import os
import pathlib
import struct
import time

folder = pathlib.Path(os.environ['SIRIUS_GEOMETRY_CAPTURE'])
dispatches, events, captures, entries, lookups = [], [], [], [], []
active = None
last_completed = None
passes, exports, creations = [], [], []
module_generation = {}
optimizer_installed = False


def memory(address, size):
    return bytes(gdb.selected_inferior().read_memory(address, size))


def reg(name):
    return int(gdb.parse_and_eval('$' + name))


def write():
    (folder / 'bindings.json').write_text(json.dumps({
        'dispatches': dispatches, 'jit_events': events, 'entry_calls': entries,
        'pointer_lookups': lookups, 'optimizer_passes': passes, 'ir_exports': exports, 'module_creations': creations,
    }, indent=2) + '\n')


def failure(error):
    (folder / 'capture-error.txt').write_text(type(error).__name__ + ': ' + str(error))
    return True


class Returned(gdb.FinishBreakpoint):
    def __init__(self, item):
        super().__init__(gdb.newest_frame(), internal=True)
        self.item = item

    def stop(self):
        global active, last_completed
        try:
            assert active is self.item
            self.item['completed'] = True
            last_completed, active = self.item, None
            write()
            return False
        except Exception as error:
            return failure(error)


class Dispatch(gdb.Breakpoint):
    def stop(self):
        global active
        try:
            assert active is None, 'overlapping synchronous application dispatch'
            count = reg('r8')
            assert count == 2
            buffers = list(struct.unpack('<II', memory(reg('rcx'), 8)))
            active = {'id': len(dispatches), 'kernel': reg('edx'),
                      'buffers': buffers, 'groups_x': reg('r9d'),
                      'completed': False, 'captured_ns': time.time_ns(),
                      'host_thread': gdb.selected_thread().global_num}
            dispatches.append(active)
            Returned(active)
            write()
            return False
        except Exception as error:
            return failure(error)


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


class JitEntry(gdb.Breakpoint):
    def __init__(self, record, sections, symbol):
        super().__init__('*' + hex(symbol['address']), internal=True)
        self.record, self.sections, self.symbol = record, sections, symbol
        self.seen = set()

    def stop(self):
        try:
            # The original control also exercises independent Endpoint/Dense
            # pairs before the modern coupled capture. They are deliberately
            # outside this observer's single-command dispatch binding scope.
            key = None if active is None else active['id']
            if key not in self.seen:
                self.seen.add(key)
                entries.append({'file': self.record['file'], 'function': self.symbol,
                                'dispatch_id': key,
                                'buffers': None if active is None else active['buffers'],
                                'kernel': None if active is None else active['kernel'],
                                'worker_thread': gdb.selected_thread().global_num,
                                'actual_entry_pc': reg('rip')})
            if active is not None and 'first_actual_dispatch_id' not in self.record:
                self.record['first_actual_dispatch_id'] = key
                self.record['maps_at_first_call'] = pathlib.Path(
                    '/proc/' + str(gdb.selected_inferior().pid) + '/maps').read_text()
            write()
            return False
        except Exception as error:
            return failure(error)


class JitCapture(gdb.Breakpoint):
    def stop(self):
        try:
            descriptor = int(gdb.parse_and_eval('(unsigned long)&__jit_debug_descriptor'))
            version, action, entry, first = struct.unpack('<IIQQ', memory(descriptor, 24))
            assert version == 1 and action in (1, 2)
            event = {'action': action, 'entry': entry, 'captured_ns': time.time_ns()}
            if entry:
                following, previous, address, size = struct.unpack('<QQQQ', memory(entry, 32))
                event.update(object_address=address, object_bytes=size)
                if action == 1:
                    assert 0 < size < 32 * 1024**2
                    raw = memory(address, size)
                    name = 'jit-%02d.elf' % len(events)
                    (folder / name).write_bytes(raw)
                    event.update(file=name, sha256=hashlib.sha256(raw).hexdigest())
                    sections, symbols = elf_sections(raw)
                    event.update(functions=symbols, sections=sections)
            events.append(event)
            write()
            return False
        except Exception as error:
            return failure(error)


class LookupReturned(gdb.FinishBreakpoint):
    def __init__(self, module, generation):
        super().__init__(gdb.newest_frame(), internal=True)
        self.module, self.generation = module, generation

    def stop(self):
        try:
            pointer = reg('rax')
            matched = [(event, symbol) for event in events if event.get('file')
                       for symbol in event['functions']
                       if symbol['address'] == pointer and symbol['name'] == 'cs_variant']
            lookups.append({'module': self.module, 'module_generation': self.generation, 'pointer': pointer, 'matched_files': [e['file'] for e, s in matched],
                            'active_dispatch_id': None if active is None else active['id']})
            assert len(matched) <= 1, 'ambiguous actual callable JIT entry'
            for event, symbol in matched:
                if 'live_sections_after_pointer_lookup' in event:
                    assert event['actual_callable_pointer'] == pointer
                    continue
                live = []
                for section in event['sections']:
                    if section['flags'] & 2 and section['address'] and section['bytes']:
                        assert section['bytes'] < 16 * 1024**2
                        raw = memory(section['address'], section['bytes'])
                        name = event['file'] + '.live-' + str(section['index']) + '.bin'
                        (folder / name).write_bytes(raw)
                        live.append({**section, 'path': name,
                                     'sha256': hashlib.sha256(raw).hexdigest()})
                event['live_sections_after_pointer_lookup'] = live
                event['actual_callable_pointer'] = pointer
                # Capture linked bytes before inserting the entry breakpoint.
                JitEntry(event, event['sections'], symbol)
            write()
            return False
        except Exception as error:
            return failure(error)


class ModuleCreated(gdb.FinishBreakpoint):
    def stop(self):
        try:
            module = reg('rax')
            assert module
            generation = len(creations)
            creations.append({'id': generation, 'module': module})
            module_generation[module] = generation
            write()
            return False
        except Exception as error:
            return failure(error)


class ModuleCreate(gdb.Breakpoint):
    def stop(self):
        global optimizer_installed
        try:
            if not optimizer_installed:
                parent_api = int(gdb.parse_and_eval('(unsigned long)&LLVMGetGlobalParent'))
                assert memory(parent_api, 16).hex() == 'f30f1efa554889e5488b47285dc36690'
                Optimizer('*' + hex(parent_api - 0x10cc070 + 0x473aad0), internal=True)
                optimizer_installed = True
            ModuleCreated(gdb.newest_frame(), internal=True)
            return False
        except Exception as error:
            return failure(error)


class PassReturned(gdb.FinishBreakpoint):
    def __init__(self, item):
        super().__init__(gdb.newest_frame(), internal=True)
        self.item = item

    def stop(self):
        try:
            self.item.update(completed=True, error_pointer=reg('rax'))
            write()
            return False
        except Exception as error:
            return failure(error)


class Optimizer(gdb.Breakpoint):
    def stop(self):
        try:
            if pathlib.Path(gdb.solib_name(reg('rip')) or '').name != 'libLLVM.so.20.1':
                return False
            pipeline = bytearray()
            for i in range(4096):
                c = memory(reg('rsi') + i, 1)
                if c == b'\0':
                    break
                pipeline.extend(c)
            else:
                raise AssertionError('unterminated optimizer pipeline')
            item = {'id': len(passes), 'module': reg('rdi'),
                    'module_generation': module_generation[reg('rdi')],
                    'export_rows_before': len((folder / 'optimizer-exports.tsv').read_text().splitlines()) if (folder / 'optimizer-exports.tsv').exists() else 0,
                    'pipeline': pipeline.decode(), 'completed': False}
            passes.append(item)
            PassReturned(item)
            write()
            return False
        except Exception as error:
            return failure(error)


class Lookup(gdb.Breakpoint):
    def stop(self):
        try:
            # Provider-specific passive query: installed LLVMGetGlobalParent
            # reads GlobalValue::Parent at +0x28. Exact implementation and
            # provider file identity are preserved with this diagnostic.
            function = reg('rsi')
            parent_api = int(gdb.parse_and_eval('(unsigned long)&LLVMGetGlobalParent'))
            assert memory(parent_api, 16).hex() == 'f30f1efa554889e5488b47285dc36690'
            module = struct.unpack('<Q', memory(function + 0x28, 8))[0]
            generation = module_generation[module]
            LookupReturned(module, generation)
            if any(e['module_generation'] == generation for e in exports):
                return False
            successful = [p['id'] for p in passes if p['module'] == module and p['module_generation'] == generation
                          and p['completed'] and p['error_pointer'] == 0]
            meta = (folder / 'optimizer-exports.tsv').read_text().splitlines()
            final_passes = [p for p in passes if p['id'] in successful and 'instcombine' in p['pipeline']]
            assert final_passes
            matching = [line.split('\t') for i, line in enumerate(meta) if i >= final_passes[-1]['export_rows_before'] and int(line.split('\t')[0],16) == module]
            assert successful and matching, 'no completed ordinary optimizer serialization for current module'
            address, name, size = matching[-1]
            raw = (folder / name).read_bytes()
            assert len(raw) == int(size) and raw
            exports.append({'module': module, 'module_generation': generation, 'function': function,
                            'successful_optimizer_pass_ids': successful, 'path': name,
                            'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest(),
                            'optimized_module_observed': True})
            write()
            return False
        except Exception as error:
            return failure(error)


class Readback(gdb.Breakpoint):
    def stop(self):
        try:
            handle = struct.unpack('<I', memory(reg('rsi'), 4))[0]
            data, size = struct.unpack('<QQ', memory(reg('rdx'), 16))
            assert last_completed is not None and last_completed['completed']
            assert handle == last_completed['buffers'][1]
            assert 0 < size < 16 * 1024**2 and size % 4 == 0
            raw = memory(data, size)
            name = 'readback-%02d.bin' % len(captures)
            (folder / name).write_bytes(raw)
            captures.append({'index': len(captures), 'handle': handle, 'path': name,
                             'bytes': size, 'sha256': hashlib.sha256(raw).hexdigest(),
                             'dispatch_id': last_completed['id']})
            (folder / 'captures.json').write_text(json.dumps(captures, indent=2) + '\n')
            assert len(captures) <= 44
            return False
        except Exception as error:
            return failure(error)


class Started(gdb.Breakpoint):
    def stop(self):
        (folder / 'inferior.json').write_text(json.dumps({'pid': gdb.selected_inferior().pid}) + '\n')
        return False


Started('main', temporary=True, internal=True)
Dispatch('*' + os.environ['SIRIUS_DISPATCH_SYMBOL'], internal=True)
# No readback callback: this direct canary returns before retained arithmetic.
JitCapture('__jit_debug_register_code', internal=True)
Lookup('LLVMGetPointerToGlobal', internal=True)

ModuleCreate('LLVMModuleCreateWithNameInContext', internal=True)

# One direct probe loads exactly one complete module. Observe the actual span,
# then join its single completed dispatch to successful optimizer/lookup/code.
class KernelLoad(gdb.Breakpoint):
    def stop(self):
        try:
            # Current x86_64 GCC ABI: expected<KernelHandle,Error> uses sret;
            # this is rsi and span pointer/count are rdx/rcx.
            count = reg('rcx')
            assert 5 < count <= 2 * 1024 * 1024
            raw = memory(reg('rdx'), count * 4)
            assert raw[:4] == bytes.fromhex('03022307')
            expected = pathlib.Path(os.environ['SIRIUS_MODULE_INPUT']).read_bytes()
            assert raw == expected
            prior = folder / 'actual-loaded-module.json'
            assert not prior.exists(), 'probe must load exactly one module'
            prior.write_text(json.dumps({'words': count, 'bytes': len(raw),
                'sha256': hashlib.sha256(raw).hexdigest(),
                'actual_span_matches_frozen_input': True}) + '\n')
            return False
        except Exception as error:
            return failure(error)

def ended(event):
    (folder / 'inferior-exit.json').write_text(json.dumps({'exit_code':getattr(event,'exit_code',None)})+'\n')
KernelLoad('*' + os.environ['SIRIUS_LOAD_SYMBOL'], internal=True)
gdb.events.exited.connect(ended)
