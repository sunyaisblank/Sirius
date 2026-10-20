from pathlib import Path
import collections,hashlib,json,re,subprocess
W=Path(__file__).resolve().parent;F=W/'capture-candidate'
verify=json.loads((F/'capture-verification.json').read_text());bindings=json.loads((F/'bindings.json').read_text())
event,=[x for x in bindings['jit_events'] if x.get('file')==verify['actual_entry']['file']]
section,=[x for x in event['live_sections_after_pointer_lookup'] if x['name']=='.ltext']
data=(F/section['path']).read_bytes();assert hashlib.sha256(data).hexdigest()==section['sha256']
functions=event['functions'];regions=[(f['address'],f['address']+f['bytes']) for f in functions]
dump=subprocess.check_output(['/usr/bin/objdump','-D','-b','binary','-m','i386:x86-64','-M','intel','--adjust-vma='+str(section['address']),str(F/section['path'])],text=True)
lines=[];instructions=[]
for line in dump.splitlines():
    match=re.match(r'\s*([0-9a-f]+):\s+(?:[0-9a-f]{2}\s+)+\s*([a-z][a-z0-9]*)\s*(.*)',line)
    if match and any(a<=int(match[1],16)<b for a,b in regions):lines.append(line);instructions.append((match[2],match[3]))
(W/'called-functions.asm').write_text('\n'.join(lines)+'\n')
counts=collections.Counter(x for x,_ in instructions)
ir=(F/verify['ordinary_ir_export']['path']).read_text()
assert hashlib.sha256(ir.encode()).hexdigest()==verify['ordinary_ir_export']['sha256']
assert counts['vdivpd']>0 and counts['vsqrtpd']>0
assert len(re.findall(r'\bfdiv <8 x double>',ir))>0 and re.search(r'@llvm.sqrt.v8f64',ir)
features=json.loads((F/'feature-verification.json').read_text())['features'];assert features['shaderFloat64']==1
baseline=json.loads((W/'baseline-ordinary-divide-review.json').read_text())['rows']['baseline']
assert counts['mfence']==baseline['fences']
result={'structural_premise_observed':True,'actual_complete_module':verify['actual_loaded_module'],'actual_enabled_features':features,'actual_called_entry':verify['actual_entry'],'baseline_capture_reused':baseline,'candidate':{'functions':functions,'total_function_bytes':sum(x['bytes'] for x in functions),'ordinary_IR_fdiv8double_sites':len(re.findall(r'\bfdiv <8 x double>',ir)),'ordinary_IR_sqrt8double_calls':len(re.findall(r'\bcall <8 x double> @llvm.sqrt.v8f64',ir)),'vdivpd_instructions':counts['vdivpd'],'vsqrtpd_instructions':counts['vsqrtpd'],'integer_divide_instructions':counts['div'],'fences':counts['mfence'],'stack_adjustments':[line for line in lines if re.search(r'\bsub\s+rsp,',line)]},'scope':'Actual ordinary called code includes packed hardware estimates and all integer fallbacks. Short-input dispatch exits before retained arithmetic; scalar corpus supplies numerical evidence. No fast-path frequency, integrated correctness or benefit claim.'}
(W/'ordinary-code-review.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result['candidate']))
