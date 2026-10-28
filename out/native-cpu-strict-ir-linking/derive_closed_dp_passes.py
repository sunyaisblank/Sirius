"""Emit fixed closed-DP pass sources only; never compile or invoke a pipeline.

scope: after the owned ordinary emit, derive the original exact FP converter.
value: after owned scope conversion, derive its original local-bank lowering.
"""
from collections import Counter
from pathlib import Path
import hashlib
import json
import re
import sys

R=Path(__file__).resolve().parents[2]
V=Path(__file__).resolve().parent
W=R/'out/native-cpu-canonical-closed-dp'
B=R/'bin/linux-gcc/native-cpu-canonical-closed-dp'
RAW=B/'dopri_phase.ll'
SOURCE=B/'ispc-source/dopri_phase.ispc'
STRICT=B/'dp.strict.bc'
MASKED='___un_3C_unu_3E_unuun_3C_unu_3E_unuunuunu'
OLD_EXPORT='sirius_simd_closed_dp_packet'
NEW_EXPORT='sirius_closed_dp_packet'
POINTER_PREFIX='sirius_closed_scope_v4_'
VALUE_PREFIX='sirius_local_value_'
RECIPES=('dense_add8','dense_multiply8','dense_divide8','dense_sqrt8',
         'tagged_add8','tagged_subtract8','tagged_multiply8','tagged_divide8','tagged_sqrt8')
ORIGINAL_MATH=set(RECIPES)-{'tagged_sqrt8'}
PINS={
 'out/native-cpu-strict-ir-linking/constrain_scope.cpp':(14483,'fa6b225b21edb880054e006504ec1c401155bbd2448c68c250be953ababf8fab'),
 'out/native-cpu-strict-ir-linking/shared_value_abi.cpp':(28405,'4cbefca757c267222e0bcde0e6818776ec058ab70fa6913d470785716abab1e6'),
 'out/native-cpu-strict-ir-linking/scope_contract.json':(9197,'adef7709ec5c94c55442840a30023a5db00972ab4d078ea625950bf8b2391cf9'),
 'out/native-cpu-strict-ir-linking/canonical_private_dp.py':(9467,'a9fbfbd69cbb4ebc1aad739fda95c39352a7d3281f1e03ca5e3ef4c96dd12508'),
 'out/native-cpu-strict-ir-linking/canonical_private_dp_contract.json':(11905,'242e1e101927473daea1c0cb8a9238d438e5f569a27d6218e01b56c6b46a674e'),
 'bin/linux-gcc/native-cpu-local-effects/math.bc':(296980,'93e076866edeb28c329a65a81dad1485a2717a8104ec4b7879002c78a9e6a6a1'),
 'bin/linux-gcc/native-cpu-local-effects/math.o':(259216,'69f9d771c5cf2d0de63e32ac70ac534961e9283a02753047784902e3ee54cdbe')}


def pin(p):
    with p.open('rb') as stream:
        return dict(bytes=p.stat().st_size,sha256=hashlib.file_digest(stream,'sha256').hexdigest())


def doc(p):return json.loads(p.read_text())


def write(p,text):
    with p.open('x') as stream:stream.write(text)


def owner(name):
    row=doc(W/name/'owner.json')
    assert row['passed'] and row['status']=='terminal' and row['child_exit']==0
    assert row['source_unchanged'] and row['whole_inputs_unchanged']
    assert not row['remaining_owned_processes'] and not row['cleanup_errors']
    return row


def authorities():
    inputs={str(Path(__file__).resolve().relative_to(R)):pin(Path(__file__).resolve())}
    for name,(size,sha) in PINS.items():
        expected=dict(bytes=size,sha256=sha)
        assert pin(R/name)==expected,name
        inputs[name]=expected
    contract=doc(V/'canonical_private_dp_contract.json')
    proof=contract['source_proof']
    for key in ('source_inverse_exact','unchanged_outside_export_binding_and_EvaluateProgram',
        'canonical_payload_unwrapped_only','canonical_tuple_records_match_exact_admitted_table',
        'canonical_root_validation_and_publication_match_table','all_validation_before_first_publication',
        'complete_private_context_memset_preserved','original_owned_row_shape_IO_guards_and_callers_unchanged',
        'sole_supported_caller_required200_lane0','full8_physical_payloads_and_tags_unchanged',
        'generic_program_fallback_removed_from_new_private_body_only'):
        assert proof[key] is True,key
    assert proof['ordered_layer_numbers']==list(range(45))
    assert (proof['original_tuple_records'],proof['layer_calls'],proof['layer_status_checks'],
            proof['validated_roots'],proof['published_roots'],proof['publication_words'])==(2088,45,45,200,200,1000)
    assert contract['new_export']==OLD_EXPORT
    assert contract['new_masked_export']==OLD_EXPORT+MASKED
    assert contract['desired_strict_export_prefix']==NEW_EXPORT
    assert pin(SOURCE)==contract['derived_source']
    for name,expected in contract['source_inputs'].items():
        assert pin(R/name)==expected,name
        inputs[name]=expected
    inputs[str(SOURCE.relative_to(R))]=pin(SOURCE)
    return contract,inputs


def inventory():
    receipt=doc(W/'inventory.json');assert receipt['output']==dict(path=str(RAW.relative_to(R)),**pin(RAW))
    assert receipt['fp_instruction_flags']=={}
    text=RAW.read_text()
    assert re.search(r'^source_filename = '+re.escape(json.dumps(str(SOURCE)))+r'$',text,re.M)
    headers=re.findall(r'^define .*$',text,re.M)
    ordinary=Counter(re.findall(r'^\s+%[^=]+ = (fadd|fsub|fmul|fdiv|frem|fneg|fcmp|fpext|fptrunc|sitofp|uitofp|fptosi|fptoui)\b',text,re.M))
    predicates=Counter(re.findall(r'= fcmp (\w+) ',text))
    calls=Counter(re.findall(r'\bcall\s+.*?@([^ (]+)\(',text))
    assert set(ordinary)<={'fadd','fsub','fcmp','fneg'},ordinary
    assert set(predicates)<={'oeq','une','uno','oge','ole','ogt','ugt'},predicates
    assert dict(ordinary)==receipt['ordinary_fp_operations']
    assert dict(predicates)==receipt['comparison_predicates']
    assert dict(calls)==receipt['call_inventory']
    assert headers==receipt['definition_headers'] and len(headers)==receipt['definitions']
    assert not re.search(r'^\s+.*\b(?:fneg|fadd|fsub|fmul|fdiv|fcmp|call)\s+(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b',text,re.M)
    assert 'strictfp' not in text and 'llvm.experimental.constrained.' not in text
    math={k:v for k,v in calls.items() if k.startswith('sirius_native_v4_')}
    assert math==receipt['math_calls'] and math
    assert {n.removeprefix('sirius_native_v4_') for n in math}<=ORIGINAL_MATH
    declarations=re.findall(r'^declare .*?@([^ (]+)\(',text,re.M)
    math_decls=[n for n in declarations if n.startswith('sirius_native_v4_')]
    assert {n.removeprefix('sirius_native_v4_') for n in math_decls}<=ORIGINAL_MATH
    assert set(math)==set(math_decls), 'Unexpected unused math declaration requires review.'
    exports=[re.search(r'@([^ (]+)',h)[1] for h in headers if not h.startswith('define internal ')]
    assert set(exports)=={OLD_EXPORT,OLD_EXPORT+MASKED} and len(exports)==2
    context_clears=sum('@llvm.memset.' in line and re.search(r', i(?:32|64)(?: range\([^)]*\))? 120832, i1 false',line) is not None
                       for line in text.splitlines())
    assert context_clears==2,'Original full physical context initialization must be explicit.'
    return receipt,dict(counts=[ordinary.get(k,0) for k in ('fadd','fsub','fcmp','fneg')]+[calls.get('llvm.fabs.v8f32',0)],
        predicates=[predicates.get(k,0) for k in ('oeq','une','uno','oge','ole','ogt','ugt')],
        definitions=len(headers),math_calls=math,math_declarations=math_decls,
        fmax_calls=calls.get('sirius_simd_fmaxf',0),context_clears=context_clears)


def patch_source(source,replacements):
    result=source
    for before,after in replacements:
        assert result.count(before)==1,(before,result.count(before))
        result=result.replace(before,after,1)
    inverse=result
    for before,after in reversed(replacements):
        assert inverse.count(after)==1,(after,inverse.count(after))
        inverse=inverse.replace(after,before,1)
    assert inverse==source
    return result


def emit_scope(contract,inputs):
    owner('derive');owner('emit')
    receipt,shape=inventory();counts=shape['counts'];preds=shape['predicates']
    source=(V/'constrain_scope.cpp').read_text()
    vector=lambda xs:'{'+', '.join(map(str,xs))+'}'
    replacements=[
      ('constexpr StringLiteral OldExport = "sirius_simd_v4_dopri_phase_packet";',f'constexpr StringLiteral OldExport = "{OLD_EXPORT}";'),
      ('constexpr StringLiteral NewExport = "sirius_scope_dp_packet";',f'constexpr StringLiteral NewExport = "{NEW_EXPORT}";'),
      ('constexpr StringLiteral NewMathPrefix = "sirius_scope_v4_";',f'constexpr StringLiteral NewMathPrefix = "{POINTER_PREFIX}";'),
      ('constexpr StringLiteral Source =\n    "/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-inline-private-bank/"\n    "ispc-source/dopri_phase.ispc";',f'constexpr StringLiteral Source = {json.dumps(str(SOURCE))};'),
      ('constexpr StringLiteral PinnedInput =\n    "/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-strict-ir-linking/"\n    "dopri_phase.ll";',f'constexpr StringLiteral PinnedInput = {json.dumps(str(RAW))};'),
      ('    "be998fd307e38ca7944aee9150db28ab88c627690873ee0343be38914df79401";',f'    "{receipt["output"]["sha256"]}";'),
      ('(*Buffer)->getBufferSize() != 53181007',f'(*Buffer)->getBufferSize() != {receipt["output"]["bytes"]}'),
      ('SmallVector<Instruction *, 721> Work;',f'SmallVector<Instruction *, {sum(counts[:3])}> Work;'),
      ('SmallVector<Function *, 56> Definitions;',f'SmallVector<Function *, {shape["definitions"]}> Definitions;'),
      ('Counts != std::array<unsigned, 5>{15, 16, 690, 2, 1276}',f'Counts != std::array<unsigned, 5>{vector(counts)}'),
      ('Predicates != std::array<unsigned, 7>{675, 1, 1, 5, 5, 2, 1}',f'Predicates != std::array<unsigned, 7>{vector(preds)}'),
      ('Definitions.size() != 56 || ExportCount != 2 || Renames.size() != 10',f'Definitions.size() != {shape["definitions"]} || ExportCount != 2 || Renames.size() != {len(shape["math_declarations"])+2}'),
      ('MathCalls != 1581 || FmaxCalls != 1',f'MathCalls != {sum(shape["math_calls"].values())} || FmaxCalls != {shape["fmax_calls"]}'),
      ('if (Name != "sirius-constrain-dp")','if (Name != "sirius-constrain-closed-dp")')]
    derived=patch_source(source,replacements)
    out=V/'constrain_closed_dp.cpp';record=V/'closed_dp_scope_contract.json'
    assert not out.exists() and not record.exists()
    inputs[str(RAW.relative_to(R))]=pin(RAW)
    for p in (W/'inventory.json',W/'derive/owner.json',W/'emit/owner.json'):
        inputs[str(p.relative_to(R))]=pin(p)
    write(out,derived)
    data=dict(phase='scope',source_only=True,pass_name='sirius-constrain-closed-dp',
        source_inputs=inputs,output=dict(path=str(out.relative_to(R)),**pin(out)),
        input=receipt['output'],actual_inventory=shape,source_inverse_exact=True,
        reused_implementation='Exact pinned constrain_scope.cpp; only source/export/input/size/SHA/full inventory/profile constants change. FP mapping/signaling/freezes/types/order/strict/dynamic/exception semantics unchanged.',
        new_exports=[NEW_EXPORT,NEW_EXPORT+MASKED],pointer_math_prefix=POINTER_PREFIX,
        future_strict_output=str(STRICT.relative_to(R)),
        private_scope='Same genuine positive canonical job, full10689-word admission/immutable views, required200/lane0, original shapes/masks/120832-byte initialization/4149-word clear, full MXCSR save/masked RNE/DAZ-FTZ-off/restore; original public/fallback objects remain exact.',
        freeze_proof='Pinned canonical source only unwraps an already-equivalent payload and leaves original full physical initialization/active-mask source unchanged. Reuse the original exact-module FP-only refinement argument; no integer/pointer/mask freeze or new lane domain.',
        qualification='No pass compile/load/conversion/codegen/runtime or transferred numerical/performance/provider acceptance.')
    write(record,json.dumps(data,indent=2)+'\n')
    return dict(output=data['output'],contract=dict(path=str(record.relative_to(R)),**pin(record)))


def emit_value(contract,inputs):
    scope_path=V/'closed_dp_scope_contract.json';scope=doc(scope_path)
    assert scope['phase']=='scope' and scope['source_inverse_exact']
    for name,expected in scope['source_inputs'].items():assert pin(R/name)==expected,name
    scope_cpp=R/scope['output']['path'];assert pin(scope_cpp)=={k:scope['output'][k] for k in ('bytes','sha256')}
    owner('scope')
    converted=doc(W/'scope-result.json')
    assert converted['output']==dict(path=str(STRICT.relative_to(R)),**pin(STRICT))
    assert converted['input']==scope['input'] and converted['scope_contract']==pin(scope_path)
    shape=scope['actual_inventory'];counts={n.removeprefix('sirius_native_v4_'):v for n,v in shape['math_calls'].items()}
    source=(V/'shared_value_abi.cpp').read_text()
    profile_start=source.index('constexpr std::array<Profile,4> Profiles{{')
    profile_end=source.index('\nstruct Recipe ',profile_start)
    profile=('constexpr std::array<Profile,1> Profiles{{\n {"sirius-closed-value-dp","'+STRICT.name+'",'+str(converted['output']['bytes'])+',\n'
             '  "'+converted['output']['sha256']+'","'+NEW_EXPORT+'",0,'+str(sum(counts.values()))+','+
             str(sum(v for k,v in counts.items() if k.startswith('tagged_')))+'}}};')
    recipe_start=source.index('constexpr std::array<Recipe,9> Recipes{{')
    recipe_end=source.index('\n\n[[noreturn]]',recipe_start)
    recipes='constexpr std::array<Recipe,9> Recipes{{\n'+',\n'.join(
        ' {"'+name+'",'+('true' if name.startswith('tagged_') else 'false')+',{'+str(counts.get(name,0))+',0,0}}' for name in RECIPES)+'}};'
    replacements=[
      ('constexpr StringLiteral Directory = "/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-strict-ir-linking/";',f'constexpr StringLiteral Directory = {json.dumps(str(B)+"/")};'),
      ('constexpr StringLiteral Prefix = "sirius_scope_v4_";',f'constexpr StringLiteral Prefix = "{POINTER_PREFIX}";'),
      ('constexpr StringLiteral SharedPrefix = "sirius_shared_value_";',f'constexpr StringLiteral SharedPrefix = "{VALUE_PREFIX}";'),
      (source[profile_start:profile_end],profile),
      (source[recipe_start:recipe_end],recipes),
      ('// remain inaccessible RW; readonly constants justify conservative Other Ref.','// remain inaccessible RW; the closed value boundary has no ordinary pointer root.'),
      ('return MemoryEffects::otherMemOnly(ModRefInfo::Ref) |\n         MemoryEffects::inaccessibleMemOnly(ModRefInfo::ModRef);','return MemoryEffects::inaccessibleMemOnly(ModRefInfo::ModRef);'),
      ('F.setDSOLocal(true);F.setCallingConv(CallingConv::Fast);','F.setDSOLocal(true);F.setCallingConv(CallingConv::X86_RegCall);'),
      ('Call->setCallingConv(CallingConv::Fast);\n  Call->addFnAttr(Attribute::StrictFP);','Call->setCallingConv(CallingConv::X86_RegCall);\n  Call->addFnAttr(Attribute::StrictFP);')]
    derived=patch_source(source,replacements)
    out=V/'closed_dp_value_abi.cpp';record=V/'closed_dp_value_contract.json'
    assert not out.exists() and not record.exists()
    for p in (scope_path,scope_cpp,W/'scope-result.json',W/'scope/owner.json',STRICT):inputs[str(p.relative_to(R))]=pin(p)
    write(out,derived)
    data=dict(phase='value',source_only=True,pass_name='sirius-closed-value-dp',source_inputs=inputs,
        output=dict(path=str(out.relative_to(R)),**pin(out)),input=converted['output'],
        exact_calls_by_recipe=counts,source_inverse_exact=True,
        reused_implementation='Exact pinned shared_value_abi.cpp local bank/capture/bounds/initialization/masked normalization/adapter implementation. Only fixed DP profile and existing cc92/inaccessible-only local math interface are rebound.',
        pointer_math_prefix=POINTER_PREFIX,value_math_prefix=VALUE_PREFIX,
        fixed_math_interface='Nine external hidden noinline strict declarations; dense10/tagged12 integer-vectors+integer mask ->5/6 integer-vector aggregate; X86_RegCall92 and full InaccessibleModRef. No new math body/profile is callable: sole accepted profile Mode0 is raw DP.',
        existing_math={name:inputs[name] for name in PINS if name.startswith('bin/linux-gcc/native-cpu-local-effects/math.')},
        private_scope=scope['private_scope'],qualification='Source only. No lower/cleanup/native, numerical/provider/frame/scientific or benefit qualification transfers.')
    write(record,json.dumps(data,indent=2)+'\n')
    return dict(output=data['output'],contract=dict(path=str(record.relative_to(R)),**pin(record)))


def main():
    assert len(sys.argv)==2 and sys.argv[1] in ('scope','value'), 'Only fixed deferred source derivation phases.'
    contract,inputs=authorities()
    result=(emit_scope if sys.argv[1]=='scope' else emit_value)(contract,inputs)
    print(json.dumps(dict(source_only=True,compiler_run=False,**result)))


if __name__=='__main__':main()
