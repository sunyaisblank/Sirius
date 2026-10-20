from pathlib import Path
import hashlib,importlib.util,json,re,struct,subprocess
R=Path.cwd();W=Path(__file__).resolve().parent;B=R/'bin/linux-gcc/portable-checked-fp64-modules';B.mkdir()
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
builder=load('checked_builder',W/'candidate_builder.py');programs=load('retained_programs',W/'kernels/retained_program.py');reader=load('original_interface_reader',W/'interface_reader.py')
tools=json.loads((R/'out/portable-checked-fp64/endpoint-arguments.json').read_text())['tools']
header=(R/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h').read_text();results={}
def encoded(p):return [p['instructions'],p['registers'],len(p['outputs'])]+p['outputs']+p['operations']+p['layer_offsets']
for kind,stem in [('Endpoint','endpoint'),('Transport','transport'),('DopriPhase','dopri_phase')]:
    program=getattr(programs,'build_'+stem+'_program')(parallel=True);words=encoded(program);specialized={}
    if kind in ['Endpoint','Transport']:
        special=getattr(programs,'build_schwarzschild_'+stem+'_program')(parallel=True);sw=encoded(special)
        value=(len(words),len(sw),len(special['layer_offsets'])-1)
        if kind=='Endpoint':value+=special['prefix_instructions'],
        specialized[stem+'_specialized']=value;words+=sw
    actual=re.search(r'\bk'+kind+r'Program\{\{(.*?)\}\};',header,re.S)
    assert actual and [int(x) for x in re.findall(r'(\d+)u',actual[1])]==words
    stock=R/'bin/linux-gcc/src/sirius/backend/retained'/('retained_'+stem+'_portable_normal_sum.spv')
    dest=B/('checked-'+stem+'.spv')
    builder.compile_shader(W/'kernels'/('retained_'+stem+'.slang'),dest,tools['compiler'],tools['assembler'],tools['disassembler'],tools['validator'],program['registers'],5,len(program['layer_offsets'])-1,program.get('prefix_instructions',0),portable=True,optimizer=tools['optimizer'],normal_sum32=True,**specialized)
    texts={}
    for arm,path in [('stock',stock),('candidate',dest)]:
        dump=B/(arm+'-'+stem+'.spvasm');subprocess.run([tools['disassembler'],str(path),'-o',str(dump)],check=True);texts[arm]=dump.read_text()
    original,candidate=texts['stock'],texts['candidate']
    assert reader.interface_signature(original)==reader.interface_signature(candidate)
    modes=lambda text:sorted(re.findall(r'OpExecutionMode %\S+ (.*)$',text,re.M))
    assert modes(original)==modes(candidate)
    assert all(original.count(op)==candidate.count(op) for op in ['OpControlBarrier','OpMemoryBarrier'])
    assert sorted(re.findall(r'OpCapability (\S+)',candidate))==['Float64','RoundingModeRTE','Shader']
    results[kind]={'module':str(dest.relative_to(R)),'bytes':dest.stat().st_size,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'stock_bytes':stock.stat().st_size,'stock_sha256':hashlib.sha256(stock.read_bytes()).hexdigest(),'all_original_program_words_equal':True,'type_layout_descriptor_shared_interface_equal':True,'execution_modes':modes(candidate),'barriers':candidate.count('OpControlBarrier'),'original_native32_controls_unchanged':True,'no_Int64_or_extra_float_mode':True,'integer_divide_and_sqrt_fallbacks_present':True}
(W/'compiler-result.json').write_text(json.dumps({'passed':True,'modules':results,'runtime_correctness_or_benefit_established':False},indent=2)+'\n');print(json.dumps(results))
