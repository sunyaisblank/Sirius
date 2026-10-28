"""Deferred, exact-source T/E strict/value pass derivation; never compile/run.

Usage: transport|endpoint scope|value. Every output is new and x-created.
The six changed Transport cast shapes are checked explicitly, never trusted.
"""
from collections import Counter
from pathlib import Path
import hashlib
import json
import re
import sys

R = Path(__file__).resolve().parents[2]
V = Path(__file__).resolve().parent
W = R / 'out/native-cpu-canonical-closed-remaining'
MASKED = '___un_3C_unu_3E_unuun_3C_unu_3E_unuunuunu'
RECIPES = ('dense_add8', 'dense_multiply8', 'dense_divide8', 'dense_sqrt8',
    'tagged_add8', 'tagged_subtract8', 'tagged_multiply8', 'tagged_divide8', 'tagged_sqrt8')
PINS = {
    'constrain_remaining.cpp': (19987, '5e70b3b62c3d5d0badf8070be3d6bcc43a7c10aeb80094b98d74b36d9065d5e7'),
    'remaining_scope_review.json': (11232, 'dbced6f6b516ef0e3e1fcbd672dab66d9c76581a2e624d73f07ba51d7f205055'),
    'constrain_scope.cpp': (14483, 'fa6b225b21edb880054e006504ec1c401155bbd2448c68c250be953ababf8fab'),
    'scope_contract.json': (9197, 'adef7709ec5c94c55442840a30023a5db00972ab4d078ea625950bf8b2391cf9'),
    'shared_value_abi.cpp': (28405, '4cbefca757c267222e0bcde0e6818776ec058ab70fa6913d470785716abab1e6'),
    'canonical_private_remaining.py': (14753, '1cf6545fcc1c11939cabd20b80d4f16023cd020d264ba8cff2c9c37ae1536237'),
    'canonical_private_remaining_contract.json': (58311, 'f2f9b0ec0fb0930485eaf0016d219327044f65602d22a9696b47e9450f98e1cc'),
    'private_finite_closure.json': (5541, '1e1da4cf26685c5f42007a5a4d574562c48ec74c012a973d94057ddf5c4f860c')}
MATH_PINS = {
    'bin/linux-gcc/native-cpu-local-effects/math.bc': (296980, '93e076866edeb28c329a65a81dad1485a2717a8104ec4b7879002c78a9e6a6a1'),
    'bin/linux-gcc/native-cpu-local-effects/math.o': (259216, '69f9d771c5cf2d0de63e32ac70ac534961e9283a02753047784902e3ee54cdbe')}
RAW_PINS = {
    'transport': (82880982, '2f409751b514fa3fd8fe9adba623a6cfd9ae3e3adad056251c9ba09b97529562'),
    'endpoint': (117820770, '2424c39a5178d82d9c6d406af62601c73a5577ae32639fc39029762fbb5c6968')}
# Complete actual source layouts; these are not the DP four-field context.
CONTEXTS = {
    'transport': dict(bytes=92192, output_words=4699,
        IR='{ ptr, [2295 x <8 x i32>], [459 x <8 x i32>], <8 x i32>, [125 x <8 x i32>] }',
        source_fields=['uniform GlobalParams_0* uniform globalParams_0;', 'uint32_t retainedScratch_0[2295];',
            'uint32_t SiriusNormalized_0[459];', 'uint32_t retainedStatus_0;', 'uint32_t retainedCoefficients_0[125];']),
    'endpoint': dict(bytes=117792, output_words=3829,
        IR='{ ptr, [7 x <8 x i32>], <8 x i32>, [3060 x <8 x i32>], [612 x <8 x i32>] }',
        source_fields=['uniform GlobalParams_0* uniform globalParams_0;', 'uint32_t nativeProjection_0[7];',
            'uint32_t retainedStatus_0;', 'uint32_t retainedScratch_0[3060];', 'uint32_t SiriusNormalized_0[612];'])}


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def doc(path):
    return json.loads(path.read_text())


def write(path, text):
    with path.open('x') as stream:
        stream.write(text)


def terminal(name):
    path = W / name / 'owner.json'
    row = doc(path)
    assert row['status'] == 'terminal' and row['passed'] and row['child_exit'] == 0
    assert row['source_unchanged'] and row['whole_inputs_unchanged']
    assert not row['remaining_owned_processes'] and not row['cleanup_errors']
    assert row['observed_births'] and all(x['absent'] for x in row['observed_births'])
    assert path.with_name('inputs-before.json').read_bytes() == path.with_name('inputs-after.json').read_bytes()


def patch_source(original, replacements):
    text = original
    for before, after in replacements:
        assert before != after and text.count(before) == 1, (before, text.count(before))
        text = text.replace(before, after, 1)
    inverse = text
    for before, after in reversed(replacements):
        assert inverse.count(after) == 1, (after, inverse.count(after))
        inverse = inverse.replace(after, before, 1)
    assert inverse == original
    return text


def authorities(stage):
    inputs = {str(Path(__file__).resolve().relative_to(R)): pin(Path(__file__).resolve())}
    for name, (size, sha) in PINS.items():
        path = V / name
        assert pin(path) == dict(bytes=size, sha256=sha), name
        inputs[str(path.relative_to(R))] = pin(path)
    for name, (size, sha) in MATH_PINS.items():
        assert pin(R / name) == dict(bytes=size, sha256=sha), name
        inputs[name] = pin(R / name)
    contract = doc(V / 'canonical_private_remaining_contract.json')
    for name, expected in contract['source_inputs'].items():
        assert pin(R / name) == expected, name
        inputs[name] = expected
    assert sum(p['original_tuple_records'] for p in contract['source_proofs'].values()) == 7292
    proof = contract['source_proofs'][stage]
    for key in ('source_inverse_to_publication_exact', 'source_inverse_to_AS_exact',
        'unchanged_outside_export_binding_and_EvaluateProgram', 'generic_fallback_removed_only_from_new_private_body',
        'canonical_tuple_records_match_exact_tables', 'complete_context_memset_and_full8_physical_packing_unchanged',
        'supported_callers_lane0', 'mixed_families_supported', 'homogeneous_all_active_guard_absent',
        'all_root_validation_before_any_publication', 'original_combined_mask_five_word_store_order'):
        assert proof[key] is True, key
    tuples, layers, roots = (3144, 159, 40) if stage == 'transport' else (4148, 195, 100)
    assert (proof['original_tuple_records'], proof['original_layer_call_sites'], proof['root_positions']) == (tuples, layers, roots)
    assert proof['publication_word_store_sites'] == roots * 5
    assert len(proof['per_layer_family_order']) == layers
    assert proof['per_layer_family_order'] == sorted(proof['per_layer_family_order'], key=lambda row: (row[0], row[1]))
    if stage == 'endpoint':
        for key in ('live_Endpoint11241_guard_unchanged', 'Endpoint_prefix_state_and_whole_state_rollback_unchanged', 'required20_and100_supported_per_lane'):
            assert proof[key] is True, key
        assert proof['required_roots'] == [20, 100]
    else:
        assert proof['required_roots'] == 40
    memory = doc(V / 'private_finite_closure.json')
    assert memory['decision'] == 'GO_REMOVE_ONLY_OTHER_READ'
    assert not memory['finite_FP_obligation']['ordinary_FP_conversion_authorized']
    assert not memory['finite_FP_obligation']['memory_none_authorized']
    for name, expected in memory['source_inputs'].items():
        assert pin(R / name) == expected, name
        inputs[name] = expected
    terminal('derive')
    derived = doc(W / 'derive-result.json')
    source = R / ('bin/linux-gcc/native-cpu-canonical-closed-' + stage + '/ispc-source/' + stage + '.ispc')
    assert derived['source_proofs'] == contract['source_proofs']
    assert pin(source) == contract['derived_sources'][stage]
    assert derived['outputs'][stage] == dict(path=str(source.relative_to(R)), **pin(source))
    body = re.search(r'struct KernelContext_0\s*\{(.*?)\};', source.read_text(), re.S)[1]
    assert [line.strip() for line in body.splitlines() if line.strip()] == CONTEXTS[stage]['source_fields']
    for path in (source, W / 'derive-result.json', W / 'derive/owner.json', W / 'derive/inputs-before.json', W / 'derive/inputs-after.json'):
        inputs[str(path.relative_to(R))] = pin(path)
    return contract, source, inputs


def inventory(stage, source):
    terminal('emit-' + stage)
    raw = source.parent.parent / (stage + '.ll')
    size, sha = RAW_PINS[stage]
    assert pin(raw) == dict(bytes=size, sha256=sha)
    receipt = doc(W / (stage + '-inventory.json'))
    assert receipt['stage'] == stage and receipt['output'] == dict(path=str(raw.relative_to(R)), **pin(raw))
    assert receipt['source'] == dict(path=str(source.relative_to(R)), **pin(source))
    assert receipt['fp_instruction_flags'] == {}
    text = raw.read_text()
    assert re.search(r'^source_filename = ' + re.escape(json.dumps(str(source))) + '$', text, re.M)
    assert '%v8_varying_KernelContext_0 = type ' + CONTEXTS[stage]['IR'] in text
    assert 'strictfp' not in text and 'llvm.experimental.constrained.' not in text
    assert not re.search(r'^\s+.*\b(?:fneg|fadd|fsub|fmul|fdiv|fcmp|call)\s+(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b', text, re.M)
    headers = re.findall(r'^define .*$', text, re.M)
    ops = Counter(re.findall(r'^\s+%[^=]+ = (fadd|fsub|fmul|fdiv|frem|fneg|fcmp|fpext|fptrunc|sitofp|uitofp|fptosi|fptoui)\b', text, re.M))
    predicates = Counter(re.findall(r'= fcmp (\w+) ', text))
    calls = Counter(re.findall(r'\bcall\s+.*?@([^ (]+)\(', text))
    allowed_ops = {'fadd', 'fsub', 'fcmp', 'fneg', 'sitofp'} if stage == 'transport' else {'fadd', 'fsub', 'fcmp', 'fneg', 'fmul'}
    assert set(ops) <= allowed_ops and set(predicates) <= {'oeq', 'une', 'uno', 'oge', 'ole', 'ogt', 'olt'}
    assert dict(ops) == receipt['ordinary_fp_operations'] and dict(predicates) == receipt['comparison_predicates']
    assert dict(calls) == receipt['call_inventory'] and headers == receipt['definition_headers'] and len(headers) == receipt['definitions']
    old = 'sirius_simd_closed_' + stage + '_packet'
    exports = [re.search(r'@([^ (]+)', h)[1] for h in headers if not h.startswith('define internal ')]
    assert len(exports) == 2 and set(exports) == {old, old + MASKED}
    declarations = re.findall(r'^declare .*?@([^ (]+)\(', text, re.M)
    external = {name: calls[name] for name in declarations}
    assert len(external) == len(declarations) and all(external.values())
    original = (V / 'constrain_remaining.cpp').read_text()
    allowed_external = set(re.findall(r'\{"([^"\n]+)", \d+\}', original))
    assert set(external) <= allowed_external, set(external) - allowed_external
    math = {name: count for name, count in calls.items() if name.startswith('sirius_native_v4_')}
    assert math == receipt['math_calls'] and math
    assert {name.removeprefix('sirius_native_v4_') for name in math} <= set(RECIPES)
    assert {name for name in declarations if name.startswith('sirius_native_v4_')} == set(math)
    context_clears = sum('@llvm.memset.' in line and re.search(r', i(?:32|64)(?: range\([^)]*\))? ' + str(CONTEXTS[stage]['bytes']) + r', i1 false', line) is not None for line in text.splitlines())
    assert context_clears == 2
    return raw, receipt, dict(ops=dict(ops), predicates=dict(predicates), definitions=len(headers),
        external_calls=external, total_calls=sum(calls.values()), math_calls=math,
        math_declarations=[n for n in declarations if n.startswith('sirius_native_v4_')],
        context=dict(CONTEXTS[stage]), full_context_zero_memsets=context_clears)


CAST_GUARD = r'''// Six exact closed-Transport casts. Never freeze an integer operand.
// Two local Rational parameters have only14 direct constant-splat callers.
// Four other vectors are existing all-lane0 broadcasts of initialized tables.
LoadInst *coefficientScalar(SIToFPInst &I, const Config &C,
                            StringMap<unsigned> &Tables) {
  if (&C != &Transport || !v8(I.getSrcTy(), false) || !v8(I.getDestTy(), true))
    reject("unexpected signed conversion", &I);
  if (auto *A = dyn_cast<Argument>(I.getOperand(0))) {
    Function *F = A->getParent();
    if (F != I.getFunction() || !F->hasLocalLinkage() ||
        F->getName() != "Rational_0___vyivyi" || F->arg_size() != 3 ||
        A->getArgNo() > 1) reject("unexpected Rational integer source", &I);
    unsigned Callers = 0;
    for (User *U : F->users()) {
      auto *Call = dyn_cast<CallInst>(U);
      if (!Call || Call->getCalledFunction() != F || Call->arg_size() != 3)
        reject("Rational escapes or has an indirect caller", &I);
      ++Callers;
      for (unsigned N = 0; N != 2; ++N) {
        auto *K = dyn_cast<Constant>(Call->getArgOperand(N));
        if (!K || !v8(K->getType(), false)) reject("Rational argument not constant8", &I);
        const ConstantInt *First = nullptr;
        for (unsigned Lane = 0; Lane != 8; ++Lane) {
          auto *Value = dyn_cast_or_null<ConstantInt>(K->getAggregateElement(Lane));
          if (!Value || Value->getSExtValue() < -64448 || Value->getSExtValue() > 64448 ||
              (First && Value->getValue() != First->getValue()))
            reject("Rational argument not defined bounded splat", &I);
          First = Value;
        }
      }
    }
    if (Callers != 14) reject("fourteen exact Rational callers required", &I);
    ++Tables[A->getArgNo() == 0 ? "Rational.arg0" : "Rational.arg1"];
    return nullptr;
  }
  auto *Shuffle = dyn_cast<ShuffleVectorInst>(I.getOperand(0));
  auto *Ins = Shuffle ? dyn_cast<InsertElementInst>(Shuffle->getOperand(0)) : nullptr;
  auto *Index = Ins ? dyn_cast<ConstantInt>(Ins->getOperand(2)) : nullptr;
  auto *Load = Ins ? dyn_cast<LoadInst>(Ins->getOperand(1)) : nullptr;
  if (!I.getFunction()->getName().starts_with("ComputeOwnedRow_0___") ||
      !Shuffle || !v8(Shuffle->getType(), false) ||
      !isa<PoisonValue>(Shuffle->getOperand(1)) || !Ins ||
      !isa<PoisonValue>(Ins->getOperand(0)) || !Index || !Index->isZero() ||
      !Load || !Load->isSimple() || !Load->getType()->isIntegerTy(32) ||
      Shuffle->getShuffleMask().size() != 8)
    reject("coefficient all-lane broadcast shape changed", &I);
  for (int Lane : Shuffle->getShuffleMask()) if (Lane != 0)
    reject("coefficient broadcast demands another lane", &I);
  auto *GEP = dyn_cast<GetElementPtrInst>(Load->getPointerOperand());
  auto *GV = GEP ? dyn_cast<GlobalVariable>(GEP->getPointerOperand()) : nullptr;
  auto *AT = GV ? dyn_cast<ArrayType>(GV->getValueType()) : nullptr;
  if (!GV || !GV->isConstant() || !GV->hasLocalLinkage() ||
      !GV->hasInitializer() || !AT || AT->getNumElements() != 7 ||
      !AT->getElementType()->isIntegerTy(32) || GEP->getNumIndices() != 1 ||
      (GV->getName() != "cn_0" && GV->getName() != "cd_0" &&
       GV->getName() != "en_0" && GV->getName() != "ed_0"))
    reject("coefficient is not an original initialized table", &I);
  for (unsigned K = 0; K != 7; ++K) {
    auto *Value = dyn_cast_or_null<ConstantInt>(GV->getInitializer()->getAggregateElement(K));
    if (!Value || Value->getSExtValue() < -339200 || Value->getSExtValue() > 339200)
      reject("coefficient undefined or outside exact FP32 integer range", &I);
  }
  // Full source/module pins preserve iv0..24: cn/cd iv+1 for iv<6 (1..6),
  // en/ed iv-18 for iv>=20 (2..6). All8 existing broadcast lanes are defined.
  ++Tables[GV->getName()];
  return nullptr;
}
'''
CAST_REWRITE = r'''    for (CastSite Site : Casts) {
      if (Site.Scalar) reject("closed Transport needs the proved full8 vector cast", Site.I);
      IRBuilder<> B(Site.I);
      B.setIsFPConstrained(true);
      B.setDefaultConstrainedExcept(fp::ebStrict);
      B.setDefaultConstrainedRounding(RoundingMode::Dynamic);
      CallInst *Vector = B.CreateConstrainedFPCast(
          Intrinsic::experimental_constrained_sitofp, Site.I->getOperand(0),
          Site.I->getType(), {}, "scope.signed.coefficient", nullptr,
          RoundingMode::Dynamic, fp::ebStrict);
      Vector->takeName(Site.I);
      Site.I->replaceAllUsesWith(Vector);
      Site.I->eraseFromParent();
    }
'''


def scope(stage, source, contract, inputs, raw, receipt, shape):
    original = (V / 'constrain_remaining.cpp').read_text()
    selected = stage.title()
    old_export, new_export = 'sirius_simd_closed_' + stage + '_packet', 'sirius_closed_' + stage + '_packet'
    pointer = 'sirius_closed_' + stage + '_scope_v4_'
    calls_start = original.index('const NamedCount ' + selected + 'Calls[] = {')
    calls_end = original.index('};', calls_start) + 2
    calls = 'const NamedCount ' + selected + 'Calls[] = {\n' + ',\n'.join('    {' + json.dumps(n) + ', ' + str(c) + '}' for n, c in shape['external_calls'].items()) + '};'
    cfg_start = original.index('const Config ' + selected + '{')
    cfg_end = original.index('};', cfg_start) + 2
    op_order, pred_order = ('fadd', 'fsub', 'fmul', 'fcmp', 'fneg', 'sitofp'), ('oeq', 'une', 'uno', 'oge', 'ole', 'ogt', 'olt')
    ops = [shape['ops'].get(n, 0) for n in op_order] + [shape['external_calls'].get('llvm.fabs.v8f32', 0)]
    preds = [shape['predicates'].get(n, 0) for n in pred_order]
    vector = lambda values: '{' + ', '.join(map(str, values)) + '}'
    pass_name = 'sirius-constrain-closed-' + stage
    config = 'const Config ' + selected + '{\n    ' + ',\n    '.join(map(json.dumps, (pass_name, str(raw), receipt['output']['sha256'], str(source)))) + ',\n    ' + json.dumps(old_export) + ', ' + json.dumps(new_export) + ',\n    ' + ', '.join((str(receipt['output']['bytes']), str(shape['definitions']), str(shape['total_calls']), vector(ops), vector(preds), selected + 'Calls')) + '};'
    footer = '''                  if (Name == Transport.Pass)
                    MPM.addPass(ConstrainRemaining(Transport));
                  else if (Name == Endpoint.Pass)
                    MPM.addPass(ConstrainRemaining(Endpoint));
                  else return false;'''
    replacements = [(original[calls_start:calls_end], calls), (original[cfg_start:cfg_end], config),
        ('constexpr StringLiteral NewMath = "sirius_scope_v4_";', 'constexpr StringLiteral NewMath = ' + json.dumps(pointer) + ';'),
        (footer, '                  if (Name != ' + selected + '.Pass) return false;\n                  MPM.addPass(ConstrainRemaining(' + selected + '));')]
    if len(shape['math_declarations']) + 2 != 11:
        replacements.append(('Renames.size() != 11 || ExternalCalls.size() != C.Calls.size()',
            'Renames.size() != ' + str(len(shape['math_declarations']) + 2) + ' || ExternalCalls.size() != C.Calls.size()'))
    cast_proof = None
    if stage == 'transport':
        assert shape['ops'].get('sitofp') == 6 and not shape['ops'].get('fmul', 0)
        start = original.index('// The four actual i32 vectors')
        end = original.index('void requireEndpointMultiply', start)
        loop_start = original.index('    for (CastSite Site : Casts) {')
        loop_end = original.index('    for (CallInst *Call : FmaxSites)', loop_start)
        guard = '''(Tables.size() != 4 || Tables.lookup("cn_0") != 1 ||
         Tables.lookup("cd_0") != 1 || Tables.lookup("en_0") != 1 ||
         Tables.lookup("ed_0") != 1)'''
        new_guard = '''(Tables.size() != 6 || Tables.lookup("cn_0") != 1 ||
         Tables.lookup("cd_0") != 1 || Tables.lookup("en_0") != 1 ||
         Tables.lookup("ed_0") != 1 || Tables.lookup("Rational.arg0") != 1 ||
         Tables.lookup("Rational.arg1") != 1)'''
        replacements += [(original[start:end], CAST_GUARD), (original[loop_start:loop_end], CAST_REWRITE), (guard, new_guard)]
        cast_proof = dict(vector_casts=6, Rational_argument_casts=2, direct_constant_splat_callers=14,
            Rational_coefficient_max_abs=64448, all_lane0_table_broadcasts=4, initialized_tables=['cn_0', 'cd_0', 'en_0', 'ed_0'],
            coefficient_max_abs=339200, exact_FP32_integer_range=True, integer_freeze=False,
            index_source_proof='Original25-iteration loop: cn/cd iv+1 only iv<6 ->1..6; en/ed iv-18 only iv>=20 ->2..6. Full pinned source/IR identity preserves these guards.',
            conversion='Same original six full8 casts -> constrained sitofp.v8f32.v8i32, dynamic/strict; every physical integer lane already defined and exactly representable, including signs. No integer/pointer/mask freeze or noundef strengthening.')
    else:
        assert not shape['ops'].get('sitofp', 0) and shape['ops'].get('fmul') == 1 and shape['predicates'].get('olt') == 3
    derived = patch_source(original, replacements)
    out, record = V / ('constrain_closed_' + stage + '.cpp'), V / ('closed_' + stage + '_scope_contract.json')
    assert not out.exists() and not record.exists()
    inputs[str(raw.relative_to(R))] = pin(raw)
    for path in (W / (stage + '-inventory.json'), W / ('emit-' + stage) / 'owner.json',
        W / ('emit-' + stage) / 'inputs-before.json', W / ('emit-' + stage) / 'inputs-after.json'):
        inputs[str(path.relative_to(R))] = pin(path)
    write(out, derived)
    data = dict(phase='scope', stage=stage, source_only=True, pass_name=pass_name,
        source_inputs=inputs, output=dict(path=str(out.relative_to(R)), **pin(out)), input=receipt['output'], actual_inventory=shape,
        source_inverse_exact=True, new_exports=[new_export, new_export + MASKED], pointer_math_prefix=pointer,
        future_strict_output=str((raw.parent / (stage + '.strict.bc')).relative_to(R)),
        source_proof=contract['source_proofs'][stage], cast_proof=cast_proof,
        reused_implementation='Pinned constrain_remaining.cpp/remaining_scope_review.json (original constrain_scope FP mapping). Only exact source/config/export/inventory bindings and authorized six fully defined Transport vector-cast checks/lowering differ. Original FP-only freezes, quiet/signaling map, types/order/CFG/memory masks, strict calls/effects and Endpoint2^-44 fmul/olt checks remain.',
        private_scope='Genuine positive current canonical Vulkan job, original complete program/admission/extents/grid/alias/immutable-view guards; mixed family/phase masks and exact physical snapshots preserved; full context initialization and active output clear per actual stage. Saved/masked/RNE/DAZ-FTZ-off/full-restored MXCSR and no observer remain mandatory. Old public/fallback/math entries stay exact.',
        future_native_unwind_gate='Before first native emission, final ordinary private entry must carry async uwtable; require exact whole-function FDE/CFA/GP/return/call checks in object and final DSO. This generator performs no context normalization/native emission.',
        qualification='Source-only derived pass; no compiler/load/conversion/numerical/performance/provider acceptance or all32 transfer.')
    write(record, json.dumps(data, indent=2) + '\n')
    return dict(output=data['output'], contract=dict(path=str(record.relative_to(R)), **pin(record)))


def value(stage, contract, inputs):
    scope_path = V / ('closed_' + stage + '_scope_contract.json')
    scoped = doc(scope_path)
    assert scoped['phase'] == 'scope' and scoped['stage'] == stage and scoped['source_inverse_exact']
    for name, expected in scoped['source_inputs'].items():
        assert pin(R / name) == expected, name
    scope_cpp = R / scoped['output']['path']
    assert pin(scope_cpp) == {k: scoped['output'][k] for k in ('bytes', 'sha256')}
    terminal('scope-' + stage)
    converted = doc(W / (stage + '-scope-result.json'))
    B = R / ('bin/linux-gcc/native-cpu-canonical-closed-' + stage)
    strict = B / (stage + '.strict.bc')
    assert converted['output'] == dict(path=str(strict.relative_to(R)), **pin(strict))
    assert converted['input'] == scoped['input'] and converted['scope_contract'] == pin(scope_path)
    counts = {n.removeprefix('sirius_native_v4_'): c for n, c in scoped['actual_inventory']['math_calls'].items()}
    original = (V / 'shared_value_abi.cpp').read_text()
    p0, p1 = original.index('constexpr std::array<Profile,4> Profiles{{'), original.index('\nstruct Recipe ')
    pass_name = 'sirius-closed-value-' + stage
    new_export = 'sirius_closed_' + stage + '_packet'
    profile = ('constexpr std::array<Profile,1> Profiles{{\n {"' + pass_name + '","' + strict.name + '",' + str(converted['output']['bytes']) + ',\n  "' + converted['output']['sha256'] + '","' + new_export + '",0,' + str(sum(counts.values())) + ',' + str(sum(c for n, c in counts.items() if n.startswith('tagged_'))) + '}}};')
    r0, r1 = original.index('constexpr std::array<Recipe,9> Recipes{{'), original.index('\n\n[[noreturn]]')
    recipes = 'constexpr std::array<Recipe,9> Recipes{{\n' + ',\n'.join(' {"' + name + '",' + ('true' if name.startswith('tagged_') else 'false') + ',{' + str(counts.get(name, 0)) + ',0,0}}' for name in RECIPES) + '}};'
    replacements = [
        ('constexpr StringLiteral Directory = "/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-strict-ir-linking/";', 'constexpr StringLiteral Directory = ' + json.dumps(str(B) + '/') + ';'),
        ('constexpr StringLiteral Prefix = "sirius_scope_v4_";', 'constexpr StringLiteral Prefix = ' + json.dumps(scoped['pointer_math_prefix']) + ';'),
        ('constexpr StringLiteral SharedPrefix = "sirius_shared_value_";', 'constexpr StringLiteral SharedPrefix = "sirius_local_value_";'),
        (original[p0:p1], profile), (original[r0:r1], recipes),
        ('// remain inaccessible RW; readonly constants justify conservative Other Ref.', '// remain inaccessible RW; the value boundary has no ordinary pointer root.'),
        ('return MemoryEffects::otherMemOnly(ModRefInfo::Ref) |\n         MemoryEffects::inaccessibleMemOnly(ModRefInfo::ModRef);', 'return MemoryEffects::inaccessibleMemOnly(ModRefInfo::ModRef);'),
        ('F.setDSOLocal(true);F.setCallingConv(CallingConv::Fast);', 'F.setDSOLocal(true);F.setCallingConv(CallingConv::X86_RegCall);'),
        ('Call->setCallingConv(CallingConv::Fast);\n  Call->addFnAttr(Attribute::StrictFP);', 'Call->setCallingConv(CallingConv::X86_RegCall);\n  Call->addFnAttr(Attribute::StrictFP);')]
    derived = patch_source(original, replacements)
    out, record = V / ('closed_' + stage + '_value_abi.cpp'), V / ('closed_' + stage + '_value_contract.json')
    assert not out.exists() and not record.exists()
    for path in (scope_path, scope_cpp, W / (stage + '-scope-result.json'), W / ('scope-' + stage) / 'owner.json',
        W / ('scope-' + stage) / 'inputs-before.json', W / ('scope-' + stage) / 'inputs-after.json', strict):
        inputs[str(path.relative_to(R))] = pin(path)
    write(out, derived)
    data = dict(phase='value', stage=stage, source_only=True, pass_name=pass_name, source_inputs=inputs,
        output=dict(path=str(out.relative_to(R)), **pin(out)), input=converted['output'],
        exact_calls_by_recipe=counts, source_inverse_exact=True, pointer_math_prefix=scoped['pointer_math_prefix'], value_math_prefix='sirius_local_value_',
        reused_implementation='Exact pinned shared_value_abi.cpp bank snapshots, alias/bounds/zero initialization/masked normalization/adapter code. Only actual fixed Mode0 raw-stage profile/recipe counts/names plus already-proved cc92/inaccessible-only local math interface are rebound.',
        fixed_math_interface='Nine external hidden noinline strict declarations:10/12 integer vectors+mask ->5/6 integer-vector aggregate, X86_RegCall92, full InaccessibleModRef and no ordinary effects. Sole accepted Mode0 raw profile; no shared-math body/profile is callable or rebuilt.',
        existing_math={name: inputs[name] for name in MATH_PINS}, private_scope=scoped['private_scope'],
        context=scoped['actual_inventory']['context'], future_native_unwind_gate=scoped['future_native_unwind_gate'],
        qualification='Source only; no broad context pass, lower/cleanup/native/numerical/provider/frame/benefit or all32 acceptance transfer.')
    write(record, json.dumps(data, indent=2) + '\n')
    return dict(output=data['output'], contract=dict(path=str(record.relative_to(R)), **pin(record)))


def main():
    assert len(sys.argv) == 3 and sys.argv[1] in ('transport', 'endpoint') and sys.argv[2] in ('scope', 'value')
    stage, phase = sys.argv[1:]
    contract, source, inputs = authorities(stage)
    raw, receipt, shape = inventory(stage, source)
    result = scope(stage, source, contract, inputs, raw, receipt, shape) if phase == 'scope' else value(stage, contract, inputs)
    print(json.dumps(dict(source_only=True, compiler_run=False, stage=stage, phase=phase, **result)))


if __name__ == '__main__':
    main()
