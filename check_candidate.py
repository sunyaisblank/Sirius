"""Bind all generated payloads and challenge the compiled floating whitelist."""
from pathlib import Path
import hashlib, importlib.util, json, re, struct, subprocess, sys
sys.dont_write_bytecode=True
R=Path.cwd(); W=Path(__file__).resolve().parent
header=R/'bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
arrays={}
for count,name,body in re.findall(r'inline constexpr std::array<std::uint32_t, (\d+)> (\w+)\{\{(.*?)\}\};',header.read_text(),re.S):
    words=[int(word[:-1],0) for word in re.findall(r'(?:0x[0-9a-fA-F]+|\d+)u',body)]
    assert len(words)==int(count)
    data=struct.pack('<'+'I'*len(words),*words)
    arrays[name]={'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
baseline=json.loads((W/'baseline.json').read_text())
assert len(arrays)==40 and arrays.keys()==baseline['arrays'].keys()
changed={name for name in arrays if arrays[name]!=baseline['arrays'][name]}
assert changed=={'kTransportPortableNormalSumShader','kEndpointPortableNormalSumShader'},changed
rev=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
assert not subprocess.check_output(['git','status','--porcelain'])
record={'source_revision':rev,'whole_header_sha256':hashlib.sha256(header.read_bytes()).hexdigest(),
    'arrays':arrays,'changed_arrays':sorted(changed),'unchanged_complete_arrays':38,
    'scope':'Exact generated complete payload comparison with clean 7f baseline; all 38 excluded arrays unchanged. Compiled modules validated; runtime/benefit pending.'}
(W/'build-and-arrays.json').write_text(json.dumps(record,indent=2)+'\n')
spec=importlib.util.spec_from_file_location('strict_builder',R/'scripts/build-retained-kernels.py')
builder=importlib.util.module_from_spec(spec); spec.loader.exec_module(builder)
cache={}
for line in (R/'bin/linux-gcc/CMakeCache.txt').read_text().splitlines():
    match=re.match(r'([^#/:][^:]*):[^=]*=(.*)$',line)
    if match: cache[match[1]]=match[2]
results={}
for stage in ['transport','endpoint','dopri_phase']:
    module=R/f'bin/linux-gcc/src/sirius/backend/retained/retained_{stage}_portable_normal_sum.spv'
    assembly=W/(stage+'-actual.spvasm')
    subprocess.run([cache['SIRIUS_SPIRV_DIS'],str(module),'-o',str(assembly)],check=True)
    text=assembly.read_text()
    assert re.findall(r'OpCapability (\S+)',text)==['Shader','RoundingModeRTE']
    entry=re.search(r'OpEntryPoint GLCompute (%\S+)',text)[1]
    assert len(re.findall(r'OpExecutionMode '+re.escape(entry)+r' RoundingModeRTE 32',text))==1
    raw=re.sub(r'^.*OpCapability RoundingModeRTE.*\n','',text,flags=re.M)
    raw=re.sub(r'^.*OpExecutionMode '+re.escape(entry)+r' RoundingModeRTE 32.*\n','',raw,flags=re.M)
    multiply=stage!='dopri_phase'
    builder.portable_normal_sum_controls(raw,upmultiply=multiply)
    rejected=[]
    def reject(label,mutant,admit=multiply):
        try: builder.portable_normal_sum_controls(mutant,upmultiply=admit)
        except ValueError as error: rejected.append({'control':label,'reason':str(error)})
        else: raise AssertionError('Invalid compiled control accepted: '+label)
    reject('preexisting_RTE_mode',text)
    reject('relaxed_arithmetic',raw+'\n OpDecorate %bad RelaxedPrecision\n')
    reject('wide_float_type',raw.replace('OpTypeFloat 32','OpTypeFloat 64'))
    if multiply:
        body=re.search(r'%RPUpMultiply(?:_\d+)? = OpFunction .*?OpFunctionEnd',raw,re.S)[0]
        owner=body.split()[0]
        operation=re.search(r'(%\S+) = OpFMul (%\S+) (%\S+) (%\S+)',body)
        assert operation
        value=operation[1]
        reject('lost_multiply_NoContraction',re.sub(r'^.*OpDecorate '+re.escape(value)+r' NoContraction.*\n','',raw,flags=re.M))
        reject('wrong_floating_owner',raw.replace(owner,'%UnexpectedProduct'))
        reject('DP_admission_of_multiply',raw,False)
        reject('wrong_multiply_operation',raw.replace(operation[0],operation[0].replace('OpFMul','OpFAdd')))
        reject('duplicate_multiply_operation',raw.replace(operation[0],operation[0]+'\n '+operation[0].replace(value,'%extraProduct',1)))
        (W/(stage+'-upmultiply-body.spvasm')).write_text(body+'\n')
    else: reject('missing_required_multiply',raw,True)
    # Keep the original FastTwoSum dependency requirement live.
    sum_body=re.search(r'%RetainedSum32(?:_\d+)? = OpFunction .*?OpFunctionEnd',raw,re.S)[0]
    links=list(re.finditer(r'(%\S+) = OpF(?:Add|Sub) %\S+ (%\S+) (%\S+)',sum_body))
    assert len(links)==3
    last=links[2]
    mutant=sum_body[:last.start()]+last[0].replace(last[3],links[0][2])+sum_body[last.end():]
    reject('broken_FastTwoSum_dependency',raw.replace(sum_body,mutant))
    results[stage]={'prequalification_checker_pass':True,'single_upmultiply_required':multiply,
        'compiled_FMul_count':len(re.findall(r'= OpFMul ',raw)), 'negative_controls':rejected,
        'module_bytes':module.stat().st_size,'module_sha256':hashlib.sha256(module.read_bytes()).hexdigest()}
(W/'compiled-controls.json').write_text(json.dumps({'source_revision':rev,'stages':results,
    'scope':'Actual emitted module whitelist/type/mode/dependency checks and finite rejecting mutations; not universal semantic equivalence or runtime qualification.'},indent=2)+'\n')
subprocess.run([sys.executable,'-B',str(W/'bind_linux.py')],check=True,cwd=R)
print(json.dumps({'changed_arrays':sorted(changed),'other_complete_arrays_unchanged':38,
    'actual_readonly_payload_joins':160,'compiled_controls':True,'rejected_mutations':sum(len(v['negative_controls']) for v in results.values())}))
