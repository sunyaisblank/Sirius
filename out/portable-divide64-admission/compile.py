"""Finish the old stopped guard on current source; reuse exact current stock bytes."""
from pathlib import Path
from collections import Counter
import hashlib,importlib.util,json,re,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/portable-divide64-admission';B.mkdir()
def load(name,p):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
builder=load('candidate_builder',W/'candidate/scripts/build-retained-kernels.py')
old=load('original_structural_reader',W/'historical-screen.py')
r=json.loads((W/'endpoint-arguments.json').read_text());t=r['tools'];dest=B/'candidate-endpoint-normal.spv'
builder.compile_shader(W/'candidate/src/sirius/kernels/retained_endpoint.slang',dest,t['compiler'],t['assembler'],t['disassembler'],t['validator'],r['registers'],5,r['layers'],r['prefix'],portable=True,optimizer=t['optimizer'],normal_sum32=True,endpoint_specialized=tuple(r['endpoint_specialized']),diagnostic_divide64=True)
assembly=dest.with_suffix('.final.spvasm');subprocess.run([t['disassembler'],str(dest),'-o',str(assembly)],check=True)
candidate=assembly.read_text();stock=(W/'reused/stock-endpoint-normal.final.spvasm').read_text()
assert (W/'reused/stock-endpoint-normal.spv').read_bytes()==(R/'bin/linux-gcc/src/sirius/backend/retained/retained_endpoint_portable_normal_sum.spv').read_bytes()
assert old.interface_signature(stock)==old.interface_signature(candidate)
modes=lambda s:sorted(re.findall(r'OpExecutionMode %\S+ (.*)$',s,re.M))
assert modes(stock)==modes(candidate)
assert all(stock.count(op)==candidate.count(op) for op in ['OpControlBarrier','OpMemoryBarrier'])
modules={}
for arm,text in [('stock',stock),('candidate',candidate)]:
 caps=re.findall(r'OpCapability (\S+)',text);expected=['Shader','RoundingModeRTE']+(['Int64'] if arm=='candidate' else [])
 assert Counter(caps)==Counter(expected)
 widths=dict(re.findall(r'(%\S+) = OpTypeInt (\d+) [01]',text));body=old.function_body(text,'PB32Divide')
 divisions=[int(widths[kind]) for kind in re.findall(r'= OpUDiv (%\S+) ',body)]
 assert (divisions==[64] and body.count('OpLoopMerge')==0) if arm=='candidate' else (divisions and set(divisions)=={32})
 if arm=='candidate':
  assert re.findall(r'OpTypeInt 64 ([01])',text)==['0']
  wide=next(identifier for identifier,width in widths.items() if width=='64')
  assert re.search(r'= OpIMul '+re.escape(wide)+r' ',body) and re.search(r'= OpISub '+re.escape(wide)+r' ',body)
 (B/(arm+'-PB32Divide.spvasm')).write_text(body)
 modules[arm]={'capabilities':caps,'integer_widths':sorted(set(widths.values())),'divide_udiv_widths':divisions,'divide_loop_count':body.count('OpLoopMerge'),'control_barriers':text.count('OpControlBarrier'),'execution_modes':modes(text)}
result={'schema':'sirius-divide64-finished-compiler-screen-v1','passed':True,'stock_compilation_reused_exactly':True,
 'candidate_module':{'path':str(dest),'bytes':dest.stat().st_size,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest()},'modules':modules,
 'complete_descriptor_and_workgroup_type_decorations_exact':True,'float_controls_modes_and_barriers_exact':True,'capability_multiset_exact_order_irrelevant':True,
 'source_argument_not_machine_checked_proof':'After unchanged special-value exits and normalization,D<=R<2D with0<D<2^24. N=R*2^26<2^51. One unsigned64 quotient q=floor(N/D) is in[2^26,2^27), and residual N-qD is in[0,D). The original8+8+8+2 radix chunks compute the same q and residual. Exact sticky-bit insertion and PB32Round are unchanged. All64-bit intermediates fit; all stored/public words remain32-bit.',
 'Vulkan_devices_pipelines_dispatches':0,'numerical_correctness_JIT_cost_or_benefit_established':False,
 'next_dependency':'Actual advertised shaderInt64 support and explicit enablement, then ordinary called64-bit divide/remainder and code/frame/JIT costs. No dispatch with an unenabled optional capability.'}
(W/'compiler-result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'passed':True,'candidate_bytes':dest.stat().st_size,'divide_udiv_widths':modules['candidate']['divide_udiv_widths'],'Vulkan_calls':0}))
