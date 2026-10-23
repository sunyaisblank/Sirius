"""GDB-only passive engine/factory observer; no MachineFunction callbacks."""
# Sourced after the original capture_jit/features/properties scripts. All
# register reads use the existing x86_64 SysV observer; no inferior calls.
cg_rows = {'engines': [], 'bridges': [], 'lookups': [], 'emissions': [], 'factories': []}
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

class CgBridgeReturned(gdb.FinishBreakpoint):
    def __init__(self, row):
        super().__init__(gdb.newest_frame(), internal=True); self.row = row
    def stop(self):
        try:
            self.row.update(completed=True, returned_bool=reg('rax') & 0xff,
                            return_sequence=cg_tick())
            assert self.row['returned_bool'] == 1
            cg_write(); return False
        except Exception as error: return failure(error)

class CgBridge(gdb.Breakpoint):
    def stop(self):
        try:
            engine, = [x for x in cg_rows['engines'] if x.get('completed') and x.get('engine') == reg('rdi')]
            completed = [p['id'] for p in passes if p['module'] == engine['module']
                         and p['module_generation'] == engine['module_generation']
                         and p['completed'] and p['error_pointer'] == 0]
            assert completed and any('instcombine' in passes[i]['pipeline'] for i in completed)
            row = {'id': len(cg_rows['bridges']), **cg_boundary(os.environ['SIRIUS_ACTUAL_PROVIDER']),
                   'engine': engine['engine'], 'module': engine['module'],
                   'module_generation': engine['module_generation'],
                   'completed_optimizer_pass_ids_at_bridge': completed, 'completed': False}
            cg_rows['bridges'].append(row); CgBridgeReturned(row)
            cg_write(); return False
        except Exception as error: return failure(error)

class CgLookupReturned(gdb.FinishBreakpoint):
    def __init__(self, row):
        super().__init__(gdb.newest_frame(), internal=True); self.row = row
    def stop(self):
        try:
            self.row.update(completed=True, pointer=reg('rax'), return_sequence=cg_tick())
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
        (os.environ['SIRIUS_POLICY_BRIDGE_SYMBOL'], CgBridge, None),
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
