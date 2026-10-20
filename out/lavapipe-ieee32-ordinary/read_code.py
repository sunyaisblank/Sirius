from pathlib import Path
import collections,hashlib,io,json,re,subprocess,zipfile
R=Path.cwd();W=Path(__file__).resolve().parent;F=W/'capture'
verification=json.loads((F/'capture-verification.json').read_text());bindings=json.loads((F/'bindings.json').read_text());event,=[x for x in bindings['jit_events'] if x.get('file')==verification['actual_entry']['file']]
section,=[x for x in event['live_sections_after_pointer_lookup'] if x['name']=='.ltext'];raw=(F/section['path']).read_bytes();assert hashlib.sha256(raw).hexdigest()==section['sha256']
regions=[(x['address'],x['address']+x['bytes']) for x in event['functions']]
dump=subprocess.check_output(['/usr/bin/objdump','-D','-b','binary','-m','i386:x86-64','-M','intel','--adjust-vma='+str(section['address']),str(F/section['path'])],text=True);lines=[];ops=[]
for line in dump.splitlines():
    m=re.match(r'\s*([0-9a-f]+):\s+(?:[0-9a-f]{2}\s+)+\s*([a-z][a-z0-9]*)\s*(.*)',line)
    if m and any(a<=int(m[1],16)<b for a,b in regions):lines.append(line);ops.append(m[2])
(W/'called-native-functions.asm').write_text('\n'.join(lines)+'\n');counts=collections.Counter(ops)
ir=(F/verification['ordinary_ir_export']['path']).read_text();assert hashlib.sha256(ir.encode()).hexdigest()==verification['ordinary_ir_export']['sha256']
oldprefix='attestations/qualification/8c2be67d/lavapipe-vector-workgroup-rejected/'
z=zipfile.ZipFile(io.BytesIO(subprocess.check_output(['git','show','43389d91:'+oldprefix+'diagnostics-and-trial.zip'])))
baseline=json.loads(z.read('out/lavapipe-vector-capture/ordinary-code-review.json'))['rows']['baseline']
stock=R/'bin/linux-gcc/src/sirius/backend/retained/retained_transport_portable_normal_sum.spv';assert hashlib.sha256(stock.read_bytes()).hexdigest()==baseline['loaded_module']['sha256']
locked_fences=sum(bool(re.search(r'\block or\s+DWORD PTR \[rsp[^]]*\],0x0\b',line)) for line in lines)
ir_fences=len(re.findall(r'^\s*fence seq_cst\b',ir,re.M))
shader=(W/'native-transport.spvasm').read_text();assert shader.count('OpControlBarrier')==9
assert counts['mfence']+locked_fences==ir_fences==baseline['memory_fences']==9 and counts['vdivps']>0 and counts['vsqrtps']>0
result={'structural_premise_observed':True,'baseline_reused_from':'43389d91','baseline':baseline,'candidate':{'module':verification['actual_complete_module'],'functions':event['functions'],'total_function_bytes':sum(x['bytes'] for x in event['functions']),'seq_cst_IR_fences':ir_fences,'mfence_instructions':counts['mfence'],'locked_stack_OR_fences':locked_fences,'integer_divide_instructions':counts['div'],'packed_binary32_divisions':counts['vdivps'],'packed_binary32_square_roots':counts['vsqrtps'],'packed_binary32_additions':counts['vaddps'],'packed_binary32_multiplications':counts['vmulps'],'stack_adjustments':[line for line in lines if re.search(r'\bsub\s+rsp,',line)],'actual_entry_MXCSR':verification['actual_entry']['actual_mxcsr']},'reuse_scope':'Same exact unmodified baseline provider/LLVM21 and complete portable Transport bytes. Historical baseline factory did not enable unused Int64; the module has no Int64 types/capability. This is structural code evidence; all performance comparisons require fresh matched current public-factory slots.','scope':'Actual complete existing native route and ordinary called code, explicit real IEEE32 support; no scientific/benefit/frame acceptance.'}
(W/'ordinary-code-review.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['candidate']))
