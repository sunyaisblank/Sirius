"""Run only inside the future owned capture, after nm and its DSOs are sealed."""
llvm_library = Path('/usr/lib/x86_64-linux-gnu/libLLVM.so.21.1').resolve(strict=True)
llvm_symbols = subprocess.check_output(
    ['/usr/bin/nm', '-D', '--defined-only', str(llvm_library)], text=True).splitlines()
provider_symbols = subprocess.check_output(
    ['/usr/bin/nm', '--defined-only', str(driver)], text=True).splitlines()

def exact_symbol(rows, expected):
    found = {line.split()[-1].split('@', 1)[0] for line in rows
             if line.split() and line.split()[-1].split('@', 1)[0] == expected}
    assert found == {expected}, (expected, sorted(found))
    return expected

codegen_environment = {
    'SIRIUS_FAST_ALLOCATOR_SYMBOL': exact_symbol(llvm_symbols, '_ZN4llvm27createFastRegisterAllocatorEv'),
    'SIRIUS_GREEDY_ALLOCATOR_SYMBOL': exact_symbol(llvm_symbols, '_ZN4llvm29createGreedyRegisterAllocatorEv'),
    # void MCJIT::generateCodeForModule(Module*): this=rdi, module=rsi.
    # Avoid the hidden-return-pointer ABI of emitObject's unique_ptr result.
    'SIRIUS_EMITTER_SYMBOL': exact_symbol(llvm_symbols, '_ZN4llvm5MCJIT21generateCodeForModuleEPNS_6ModuleE'),
    'SIRIUS_POLICY_BRIDGE_SYMBOL': exact_symbol(provider_symbols, 'lp_set_mcjit_codegen_none_after_ir'),
    'SIRIUS_LLVM_LIBRARY': str(llvm_library),
}
(folder/'codegen-symbols.json').write_text(json.dumps({
    'provider': {'path': str(driver), 'sha256': hashlib.sha256(driver.read_bytes()).hexdigest()},
    'LLVM': {'path': str(llvm_library), 'sha256': hashlib.sha256(llvm_library.read_bytes()).hexdigest()},
    'symbols': codegen_environment,
}, indent=2)+'\n')
