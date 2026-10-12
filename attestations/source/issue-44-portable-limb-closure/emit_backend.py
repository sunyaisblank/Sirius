"""Re-emit protected optimized IR offline; output equality remains mandatory."""
from pathlib import Path
import ctypes as C
import hashlib
import json
import resource
import sys

sys.dont_write_bytecode = True
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
resource.setrlimit(resource.RLIMIT_FSIZE, (64 * 1024 * 1024, 64 * 1024 * 1024))
root = Path.cwd()
work = Path(__file__).resolve().parent
output = root / 'bin/linux-gcc/portable-execution-attribution'
assert sys.argv[1:] == ['endpoint']
source = root / 'attestations/software-vulkan/b5d2abf/optimized-ir-stage-code-binding/capture/optimized-05.ll'
library = Path('/usr/lib/x86_64-linux-gnu/libLLVM.so.20.1')
assert hashlib.sha256(library.read_bytes()).hexdigest() == '4b93481bd1c9d42c0951275dd275d4a9d30db78c309ae2cc0405669bc6dbbf8a'
data = source.read_bytes()
assert hashlib.sha256(data).hexdigest() == '3271b6d61996590a4c74f822659bbe8c33cd3cb0b91f4bdf836ae71a95a290b7'
llvm = C.CDLL(str(library))
pointer, integer, string = C.c_void_p, C.c_int, C.c_char_p

def api(name, result, *arguments):
    function = getattr(llvm, name)
    function.restype = result
    function.argtypes = arguments
    return function

dispose_message = api('LLVMDisposeMessage', None, pointer)
def text(name):
    value = api(name, pointer)()
    assert value
    result = C.string_at(value).decode()
    dispose_message(value)
    return result

for name in ('LLVMInitializeX86TargetInfo', 'LLVMInitializeX86Target',
             'LLVMInitializeX86TargetMC', 'LLVMInitializeX86AsmPrinter'):
    api(name, None)()
triple, cpu = text('LLVMGetDefaultTargetTriple'), text('LLVMGetHostCPUName')
host_features = text('LLVMGetHostCPUFeatures')
available = {v[1:]: v[0] == '+' for v in host_features.split(',')}
# Mesa explicitly selects this feature subset from util_get_cpu_caps(). Here
# its current LLVM host projection is a reconstruction, accepted only if all
# complete emitted function/constant/relocation bytes match the protected JIT.
names = ('sse', 'sse2', 'sse3', 'ssse3', 'sse4.1', 'sse4.2', 'avx', 'f16c',
         'fma', 'avx2', 'avx512f', 'avx512cd', 'avx512bw', 'avx512dq',
         'avx512vl', 'avx512vbmi')
assert all(name in available for name in names)
features = ','.join(['+64bit', *(('+' if available[name] else '-') + name for name in names)])
target, error = pointer(), pointer()
assert api('LLVMGetTargetFromTriple', integer, string, C.POINTER(pointer), C.POINTER(pointer))(
    triple.encode(), C.byref(target), C.byref(error)) == 0
# LLVM20 enum values: default backend level2, static relocation1, JIT-default
# code model1. No additional IR optimizer or executable code is run.
machine = api('LLVMCreateTargetMachine', pointer, pointer, string, string, string,
              integer, integer, integer)(target, triple.encode(), cpu.encode(), features.encode(), 2, 1, 1)
assert machine
api('LLVMSetTargetMachineAsmVerbosity', None, pointer, integer)(machine, 1)
context = api('LLVMContextCreate', pointer)()
buffer = api('LLVMCreateMemoryBufferWithMemoryRangeCopy', pointer, string, C.c_size_t, string)
parse = api('LLVMParseIRInContext', integer, pointer, pointer, C.POINTER(pointer), C.POINTER(pointer))
emit = api('LLVMTargetMachineEmitToFile', integer, pointer, pointer, string, integer, C.POINTER(pointer))
settings = {'stage': 'endpoint', 'library_sha256': hashlib.sha256(library.read_bytes()).hexdigest(),
            'input_sha256': hashlib.sha256(data).hexdigest(), 'triple': triple, 'cpu': cpu,
            'features': features, 'LLVM_host_feature_projection': host_features,
            'backend_optlevel': 2, 'reloc_enum': 1, 'code_model_enum': 1,
            'original_producer_options_uniquely_recovered': False,
            'scope': 'Offline reconstruction only. A complete artifact equality check is required '
                     'before offset correspondence; no GPU execution, timing or source operation attribution.'}
(work / 'emission-settings.json').write_text(json.dumps(settings, indent=2) + '\n')
try:
    for suffix, kind in (('s', 0), ('o', 1)):
        destination = output / ('endpoint.' + suffix)
        assert not destination.exists()
        module, error = pointer(), pointer()
        memory = buffer(data, len(data), b'protected-endpoint')
        assert parse(context, memory, C.byref(module), C.byref(error)) == 0, (
            C.string_at(error).decode() if error else 'IR parse failure')
        try:
            api('LLVMSetTarget', None, pointer, string)(module, triple.encode())
            assert emit(machine, module, str(destination).encode(), kind, C.byref(error)) == 0, (
                C.string_at(error).decode() if error else 'code generation failure')
            assert destination.is_file() and 0 < destination.stat().st_size < 64 * 1024 * 1024
        finally:
            api('LLVMDisposeModule', None, pointer)(module)
finally:
    api('LLVMContextDispose', None, pointer)(context)
    api('LLVMDisposeTargetMachine', None, pointer)(machine)
print(json.dumps({'emission_complete': True, 'settings': settings}))
