"""GDB-only passive engine/factory observer; no MachineFunction callbacks."""
# Sourced after the original capture_jit/features/properties scripts. All
# register reads use the existing x86_64 SysV observer; no inferior calls.
cg_rows = {'engines': [], 'bridges': [], 'lookups': [], 'emissions': [], 'factories': [], 'clear_bindings': [], 'clears': []}
cg_clock = 0
cg_active_emission = {}

def cg_tick():
    global cg_clock
    cg_clock += 1
    return cg_clock

def cg_write():
    (folder/'machine-policy.json').write_text(json.dumps(cg_rows, indent=2)+'\n')

def cg_boundary(expected):
    library = gdb.solib_name(reg('rip'))
    assert library and pathlib.Path(library).resolve(strict=True) == pathlib.Path(expected).resolve(strict=True)
    return {'pc': reg('rip'), 'library': library,
            'thread': gdb.selected_thread().global_num, 'sequence': cg_tick()}

class CgEngineReturned(gdb.FinishBreakpoint):
    def __init__(self, row, output):
        super().__init__(gdb.newest_frame(), internal=True)
        self.row, self.output = row, output
    def stop(self):
        try:
            assert reg('eax') & 0xffffffff == 0
            self.row.update(engine=struct.unpack('<Q', memory(self.output, 8))[0],
                            completed=True, return_sequence=cg_tick())
            assert self.row['engine']
            cg_write(); return False
        except Exception as error: return failure(error)

class CgEngine(gdb.Breakpoint):
    def stop(self):
        try:
            # Actual Mesa helper: OutJIT, OutCode, cache, M, MM, OptLevel.
            row = {'id': len(cg_rows['engines']), **cg_boundary(os.environ['SIRIUS_ACTUAL_PROVIDER']),
                   'module': reg('rcx'), 'module_generation': module_generation[reg('rcx')],
                   'initial_OptLevel': reg('r9d') & 0xffffffff, 'completed': False}
            assert row['initial_OptLevel'] == 2
            cg_rows['engines'].append(row); CgEngineReturned(row, reg('rdi'))
            cg_write(); return False
        except Exception as error: return failure(error)

class CgLookupReturned(gdb.FinishBreakpoint):
    def __init__(self, row):
        super().__init__(gdb.newest_frame(), internal=True); self.row = row
    def stop(self):
        try:
            self.row.update(completed=True, pointer=reg('rax'), return_sequence=cg_tick())
            cg_attach_clear_breakpoints(self.row)
            cg_write(); return False
        except Exception as error: return failure(error)

class CgLookup(gdb.Breakpoint):
    def stop(self):
        try:
            module = struct.unpack('<Q', memory(reg('rsi')+0x28, 8))[0]
            row = {'id': len(cg_rows['lookups']), **cg_boundary(os.environ['SIRIUS_LLVM_LIBRARY']),
                   'engine': reg('rdi'), 'function': reg('rsi'), 'module': module,
                   'module_generation': module_generation[module], 'completed': False}
            # The unchanged Lookup observer authenticates GlobalValue::Parent
            # offset against the installed LLVMGetGlobalParent implementation.
            cg_rows['lookups'].append(row); CgLookupReturned(row)
            cg_write(); return False
        except Exception as error: return failure(error)

class CgEmissionReturned(gdb.FinishBreakpoint):
    def __init__(self, row):
        super().__init__(gdb.newest_frame(), internal=True); self.row = row
    def stop(self):
        try:
            stack = cg_active_emission[self.row['thread']]
            assert stack.pop() == self.row['id']
            self.row.update(completed=True, return_sequence=cg_tick())
            cg_write(); return False
        except Exception as error: return failure(error)

class CgEmission(gdb.Breakpoint):
    def stop(self):
        try:
            row = {'id': len(cg_rows['emissions']), **cg_boundary(os.environ['SIRIUS_LLVM_LIBRARY']),
                   'engine': reg('rdi'), 'module': reg('rsi'),
                   'module_generation': module_generation[reg('rsi')], 'completed': False}
            cg_rows['emissions'].append(row)
            cg_active_emission.setdefault(row['thread'], []).append(row['id'])
            CgEmissionReturned(row); cg_write(); return False
        except Exception as error: return failure(error)

class CgFactory(gdb.Breakpoint):
    def __init__(self, symbol, kind):
        super().__init__(symbol, internal=True); self.kind = kind
    def stop(self):
        try:
            boundary = cg_boundary(os.environ['SIRIUS_LLVM_LIBRARY'])
            stack = cg_active_emission.get(boundary['thread'], [])
            emission = None if not stack else cg_rows['emissions'][stack[-1]]
            cg_rows['factories'].append({**boundary, 'allocator': self.kind,
                'emission_id': None if emission is None else emission['id'],
                'engine': None if emission is None else emission['engine'],
                'module_generation': None if emission is None else emission['module_generation']})
            cg_write(); return False
        except Exception as error: return failure(error)

cg_clear_addresses = set()

class CgClearReturned(gdb.FinishBreakpoint):
    def __init__(self, row):
        super().__init__(gdb.newest_frame(), internal=True); self.row = row
    def stop(self):
        try:
            raw = memory(self.row['destination'], self.row['bytes'])
            assert raw == b'\0' * self.row['bytes']
            assert reg('rax') == self.row['destination']
            self.row.update(completed=True, complete_zero_bytes_observed=True,
                            resulting_sha256=hashlib.sha256(raw).hexdigest(), return_sequence=cg_tick())
            cg_write(); return False
        except Exception as error: return failure(error)

class CgClear(gdb.Breakpoint):
    def __init__(self, binding):
        super().__init__('*'+hex(binding['address']), internal=True); self.binding = binding
    def stop(self):
        try:
            caller = struct.unpack('<Q', memory(reg('rsp'),8))[0]
            if not any(x['address'] <= caller < x['address']+x['bytes'] for x in self.binding['functions']):
                return False
            size = int(os.environ['SIRIUS_BULK_ZERO_BYTES'])
            if reg('rdx') != size: return False
            boundary = cg_boundary(os.environ['SIRIUS_LIBC_LIBRARY'])
            if self.binding['symbol'] == 'memset': assert reg('esi') & 0xff == 0
            else: assert memory(reg('rsi'),size) == b'\0'*size
            assert active is not None and active['id'] == 0
            row = {**boundary, 'symbol':self.binding['symbol'], 'address':self.binding['address'],
                   'destination':reg('rdi'), 'bytes':size, 'caller_return_pc':caller,
                   'dispatch_id':active['id'], 'jit_file':self.binding['jit_file'],
                   'module_generation':self.binding['module_generation'], 'completed':False}
            cg_rows['clears'].append(row); CgClearReturned(row); self.enabled = False
            cg_write(); return False
        except Exception as error: return failure(error)

def cg_attach_clear_breakpoints(lookup):
    # Read resolved absolute relocations from the actual registered JIT ELF.
    # Do not confuse libc's IFUNC resolver st_value or a provider PLT with the
    # implementation to which MCJIT actually linked the emitted call.
    for event in events:
        if not event.get('file'): continue
        selected = [x for x in event['functions'] if x['name']=='cs_variant' and x['address']==lookup['pointer']]
        if not selected: continue
        raw = (folder/event['file']).read_bytes()
        sections = event['sections']
        for relocation in sections:
            if relocation['type'] != 4: continue
            assert relocation['entry_bytes'] == 24
            symbols = sections[relocation['link']]; strings = sections[symbols['link']]
            names = raw[strings['offset']:strings['offset']+strings['bytes']]
            target = sections[relocation['info']]
            for position in range(relocation['offset'],relocation['offset']+relocation['bytes'],24):
                offset, info, addend = struct.unpack_from('<QQq',raw,position)
                index, kind = info>>32, info & 0xffffffff
                entry = symbols['offset']+index*symbols['entry_bytes']
                name, sym_info, other, shndx, value, length = struct.unpack_from('<IBBHQQ',raw,entry)
                symbol = names[name:].split(b'\0',1)[0].decode()
                if symbol not in {'memset','memcpy'}: continue
                assert shndx == 0 and kind == 1 and addend == 0, (symbol,kind,addend)
                site = target['address']+offset
                address = struct.unpack('<Q',memory(site,8))[0]
                library = gdb.solib_name(address)
                assert library and pathlib.Path(library).resolve(strict=True) == pathlib.Path(os.environ['SIRIUS_LIBC_LIBRARY']).resolve(strict=True)
                binding = {'symbol':symbol, 'address':address, 'relocation_site':site,
                    'relocation_kind':kind, 'addend':addend, 'library':library,
                    'jit_file':event['file'], 'module_generation':lookup['module_generation'],
                    'functions':event['functions']}
                key = (address, event['file'], lookup['module_generation'])
                if key not in cg_clear_addresses:
                    cg_rows['clear_bindings'].append(binding)
                    CgClear(binding); cg_clear_addresses.add(key)
    cg_write()

cg_installed = set()
cg_offsets = json.loads(os.environ['SIRIUS_MACHINE_SYMBOL_OFFSETS'])
def cg_install(event=None):
    # Bind actual frozen ELF symbol values to that object's zero-offset mapping.
    # GDB's plain &symbol lookup can select a provider PLT stub for LLVM APIs.
    # Retry on ordinary objfile loads; never call the inferior to resolve names.
    maps = pathlib.Path('/proc/'+str(gdb.selected_inferior().pid)+'/maps')
    if not maps.exists(): return
    bases = {}; executable = {}
    for line in maps.read_text().splitlines():
        fields = line.split(maxsplit=5)
        if len(fields) != 6 or not fields[5].startswith('/'): continue
        path = str(pathlib.Path(fields[5]).resolve())
        if path not in {x['library'] for x in cg_offsets.values()}: continue
        assert pathlib.Path(path).is_file()
        start, end = [int(x, 16) for x in fields[0].split('-')]
        if int(fields[2], 16) == 0: bases.setdefault(path, set()).add(start)
        if 'x' in fields[1]: executable.setdefault(path, []).append((start, end))
    definitions = [
        (os.environ['SIRIUS_ENGINE_CREATOR_SYMBOL'], CgEngine, None),
        ('LLVMGetPointerToGlobal', CgLookup, None),
        (os.environ['SIRIUS_EMITTER_SYMBOL'], CgEmission, None),
        (os.environ['SIRIUS_FAST_ALLOCATOR_SYMBOL'], CgFactory, 'fast'),
        (os.environ['SIRIUS_GREEDY_ALLOCATOR_SYMBOL'], CgFactory, 'greedy'),
    ]
    for symbol, cls, kind in definitions:
        if symbol in cg_installed: continue
        binding = cg_offsets[symbol]
        candidates = bases.get(binding['library'], set())
        if not candidates: continue
        assert len(candidates) == 1, (symbol, candidates)
        address = next(iter(candidates)) + binding['offset']
        assert any(start <= address < end for start, end in executable.get(binding['library'], [])), symbol
        library = gdb.solib_name(address)
        if not library: continue
        assert str(pathlib.Path(library).resolve(strict=True)) == binding['library'], (symbol, library)
        if kind is None: cls('*'+hex(address), internal=True)
        else: cls('*'+hex(address), kind)
        cg_installed.add(symbol)
gdb.events.new_objfile.connect(cg_install)
cg_install()
cg_write()
