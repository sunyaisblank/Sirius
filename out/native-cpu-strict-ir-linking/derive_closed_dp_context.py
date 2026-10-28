"""Deferred exact context-pass source preparation; no compiler/pipeline execution."""
from pathlib import Path
import hashlib
import json
import sys

R = Path(__file__).resolve().parents[2]
V = Path(__file__).resolve().parent
W = R / 'out/native-cpu-canonical-closed-dp'
B = R / 'bin/linux-gcc/native-cpu-canonical-closed-dp'
INPUT = B / 'dp.value.bc'
SOURCE_CONTRACT = V / 'canonical_private_dp_contract.json'
SOURCE_CONTRACT_PIN = dict(bytes=11905, sha256='242e1e101927473daea1c0cb8a9238d438e5f569a27d6218e01b56c6b46a674e')
VALUE_PIN = dict(bytes=10646600, sha256='3f69125b95e0a6282a3e7dfa2d7617a0ea933af22dd2aa411314a8b1c76a2cf0')

CPP = r'''// Generated for one exact closed DP value module. No math definition is linked.
#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/ADT/StringExtras.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/Config/llvm-config.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Operator.h"
#include "llvm/IR/PassManager.h"
#include "llvm/IR/Verifier.h"
#include "llvm/IRReader/IRReader.h"
#include "llvm/Passes/PassBuilder.h"
#include "llvm/Plugins/PassPlugin.h"
#include "llvm/Support/ErrorHandling.h"
#include "llvm/Support/MemoryBuffer.h"
#include "llvm/Support/ModRef.h"
#include "llvm/Support/SHA256.h"
#include "llvm/Support/SourceMgr.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Transforms/IPO/AlwaysInliner.h"
#include "llvm/Transforms/IPO/GlobalDCE.h"
#include <map>
#include <string>

#if LLVM_VERSION_MAJOR != 23 || LLVM_VERSION_MINOR != 1 || LLVM_VERSION_PATCH != 2
#error "Requires the pinned LLVM23.1.2 API."
#endif
using namespace llvm;
namespace {
constexpr StringLiteral Input = @INPUT@;
constexpr StringLiteral InputSHA = @SHA@;
constexpr uint64_t InputBytes = @BYTES@, ContextBytes = 120832;
constexpr StringLiteral EntryName = "sirius_closed_dp_packet";
constexpr StringLiteral MaskedName = "sirius_closed_dp_packet___un_3C_unu_3E_unuun_3C_unu_3E_unuunuunu";
constexpr StringLiteral MathPrefix = "sirius_local_value_";

[[noreturn]] void reject(const Twine &Why, const Instruction *I = nullptr) {
  errs() << "sirius-closed-dp-context: " << Why << '\n';
  if (I) { errs() << I->getFunction()->getName() << ": "; I->print(errs()); errs() << '\n'; }
  report_fatal_error("closed DP context normalization refused", false);
}
std::string digest(StringRef S) { SHA256 H; H.update(S); return toHex(H.final(), true); }
std::string fingerprint(const Module &M) {
  std::string S; raw_string_ostream O(S); M.print(O, nullptr); O.flush(); return digest(S);
}
void requirePinned(const Module &M) {
  auto Buffer = MemoryBuffer::getFile(Input);
  if (!Buffer || (*Buffer)->getBufferSize() != InputBytes || digest((*Buffer)->getBuffer()) != InputSHA)
    reject("complete original value input pin differs");
  LLVMContext C; SMDiagnostic D;
  auto Reference = parseIR((*Buffer)->getMemBufferRef(), D, C);
  if (!Reference) reject("cannot independently parse pinned input");
  Reference->setModuleIdentifier(M.getModuleIdentifier());
  auto Expected = fingerprint(*Reference);
  if (fingerprint(M) != Expected) reject("full canonical module fingerprint differs");
  errs() << "sirius-closed-dp-context: input-sha256=" << InputSHA << " module-fingerprint=" << Expected << '\n';
}
using Inventory = std::map<std::string, unsigned>;
Inventory fpInventory(const Module &M) {
  Inventory Out;
  for (const Function &F : M) for (const Instruction &I : instructions(F)) {
    if (auto *Op = dyn_cast<FPMathOperator>(&I)) if (Op->getFastMathFlags().any()) reject("unexpected fast math", &I);
    if (auto *II = dyn_cast<ConstrainedFPIntrinsic>(&I)) ++Out[II->getCalledFunction()->getName().str()];
    else if (I.getOpcode() == Instruction::FNeg) ++Out["representation.fneg"];
    else if (I.getOpcode() == Instruction::FAdd || I.getOpcode() == Instruction::FSub ||
             I.getOpcode() == Instruction::FMul || I.getOpcode() == Instruction::FDiv ||
             I.getOpcode() == Instruction::FRem || I.getOpcode() == Instruction::FCmp ||
             I.getOpcode() == Instruction::FPTrunc || I.getOpcode() == Instruction::FPExt ||
             I.getOpcode() == Instruction::SIToFP || I.getOpcode() == Instruction::UIToFP ||
             I.getOpcode() == Instruction::FPToSI || I.getOpcode() == Instruction::FPToUI)
      reject("ordinary FP operation outside strict copy", &I);
    if (auto *Call = dyn_cast<CallInst>(&I)) if (auto *Target = Call->getCalledFunction()) {
      StringRef Name = Target->getName();
      if (Name.starts_with(MathPrefix) || Name == "sirius_simd_fmaxf" ||
          Name.starts_with("llvm.fabs.") || Name.starts_with("llvm.is.fpclass.")) ++Out[Name.str()];
    }
  }
  return Out;
}
bool words(Type *T) {
  auto *Vec = dyn_cast<FixedVectorType>(T);
  return Vec && Vec->getNumElements() == 8 && Vec->getElementType()->isIntegerTy(32);
}
bool maskedWords(Type *T) {
  auto *Vec = dyn_cast<FixedVectorType>(T);
  return Vec && Vec->getNumElements() == 8 &&
    (Vec->getElementType()->isIntegerTy(32) || Vec->getElementType()->isFloatTy());
}
uint64_t sizeOf(Type *T, const DataLayout &DL) {
  auto Size = DL.getTypeStoreSize(T); if (Size.isScalable()) reject("scalable context access");
  return Size.getFixedValue();
}
bool contextType(Type *T, const DataLayout &DL) {
  auto *S = dyn_cast<StructType>(T);
  if (!S || S->isOpaque() || S->isPacked() || S->getNumElements() != 4 ||
      !S->getElementType(0)->isPointerTy() || !words(S->getElementType(1))) return false;
  auto *Scratch = dyn_cast<ArrayType>(S->getElementType(2));
  auto *Tags = dyn_cast<ArrayType>(S->getElementType(3));
  if (!Scratch || !Tags || Scratch->getNumElements() != 3145 || Tags->getNumElements() != 629 ||
      !words(Scratch->getElementType()) || !words(Tags->getElementType())) return false;
  auto *Layout = DL.getStructLayout(S);
  return DL.getTypeAllocSize(S).getFixedValue() == ContextBytes &&
    Layout->getElementOffset(0) == 0 && Layout->getElementOffset(1) == 32 &&
    Layout->getElementOffset(2) == 64 && Layout->getElementOffset(3) == 100704;
}
void bounds(uint64_t Offset, uint64_t Width, const Instruction *I) {
  if (Offset > ContextBytes || Width > ContextBytes - Offset) reject("context access out of exact bounds", I);
}

struct Context {
  AllocaInst *A;
  MemSetInst *Zero = nullptr;
  SmallVector<Instruction *, 128> Memory;
  SmallVector<IntrinsicInst *, 128> Masked;
};
Context inspect(Function &F, const DataLayout &DL) {
  AllocaInst *A = nullptr;
  for (Instruction &I : instructions(F)) if (auto *Candidate = dyn_cast<AllocaInst>(&I))
    if (contextType(Candidate->getAllocatedType(), DL)) {
      if (A) reject("more than one owned KernelContext after masked-export DCE", Candidate);
      auto *Count = dyn_cast<ConstantInt>(Candidate->getArraySize());
      if (!Count || !Count->isOne() || Candidate->getAddressSpace() != 0) reject("unexpected context allocation", Candidate);
      A = Candidate;
    }
  if (!A) reject("single exact owned KernelContext missing");
  Context Out{A, nullptr, {}, {}};
  SmallVector<Value *, 128> Pending{A}; SmallPtrSet<Value *, 32> Seen;
  while (!Pending.empty()) {
    Value *P = Pending.pop_back_val(); if (!Seen.insert(P).second) continue;
    int64_t SignedOffset = 0;
    if (GetPointerBaseWithConstantOffset(P, SignedOffset, DL) != A || SignedOffset < 0)
      reject("dynamic/nonlocal context address", dyn_cast<Instruction>(P));
    uint64_t Offset = uint64_t(SignedOffset); bounds(Offset, 0, dyn_cast<Instruction>(P));
    for (User *U : P->users()) {
      if (auto *G = dyn_cast<GetElementPtrInst>(U)) {
        if (G->getPointerOperand() != P) reject("context pointer used as index", G);
        Pending.push_back(G); continue;
      }
      if (auto *Load = dyn_cast<LoadInst>(U)) {
        if (Load->getPointerOperand() != P || Load->isAtomic() || Load->isVolatile()) reject("nonordinary context load", Load);
        uint64_t Width = sizeOf(Load->getType(), DL); bounds(Offset, Width, Load);
        if (Load->getType()->isPointerTy() && (Offset != 0 || Width != 8)) reject("unexpected context pointer field", Load);
        Out.Memory.push_back(Load);
        // Offset0 contains the real GlobalParams view. Do NOT walk its loaded
        // external pointer or normalize any view/input/output memory.
        continue;
      }
      if (auto *Store = dyn_cast<StoreInst>(U)) {
        if (Store->getPointerOperand() != P || Store->getValueOperand() == P || Store->isAtomic() || Store->isVolatile())
          reject("context address capture or nonordinary store", Store);
        uint64_t Width = sizeOf(Store->getValueOperand()->getType(), DL); bounds(Offset, Width, Store);
        if (Store->getValueOperand()->getType()->isPointerTy() && (Offset != 0 || Width != 8))
          reject("unexpected stored context pointer field", Store);
        Out.Memory.push_back(Store); continue;
      }
      if (auto *Set = dyn_cast<MemSetInst>(U)) {
        auto *N = dyn_cast<ConstantInt>(Set->getLength());
        if (!N || Set->isVolatile() || Set->getRawDest() != P) reject("unknown context memset", Set);
        bounds(Offset, N->getZExtValue(), Set);
        if (Offset == 0 && N->getZExtValue() == ContextBytes &&
            isa<ConstantInt>(Set->getValue()) && cast<ConstantInt>(Set->getValue())->isZero()) {
          if (Out.Zero && Out.Zero != Set) reject("multiple full context zero initializations", Set);
          Out.Zero = Set;
        }
        Out.Memory.push_back(Set); continue;
      }
      auto *II = dyn_cast<IntrinsicInst>(U);
      if (!II) reject("context pointer escapes to callback/return/unsupported operator", dyn_cast<Instruction>(U));
      if ((II->getIntrinsicID() == Intrinsic::lifetime_start || II->getIntrinsicID() == Intrinsic::lifetime_end) &&
          II->arg_size() == 1 && II->getArgOperand(0) == P && Offset == 0) continue;
      bool Load = II->getIntrinsicID() == Intrinsic::masked_load;
      bool Store = II->getIntrinsicID() == Intrinsic::masked_store;
      if ((!Load && !Store) || II->arg_size() != 3 || II->getArgOperand(Load ? 0 : 1) != P)
        reject("unknown context-derived intrinsic", II);
      Type *T = Load ? II->getType() : II->getArgOperand(0)->getType();
      auto *Mask = dyn_cast<FixedVectorType>(II->getArgOperand(Load ? 1 : 2)->getType());
      if (!maskedWords(T) || !Mask || Mask->getNumElements() != 8 || !Mask->getElementType()->isIntegerTy(1) ||
          sizeOf(T, DL) != 32 || Offset < 32 || Offset % 32)
        reject("masked context access is not an exact32-byte field", II);
      bounds(Offset, 32, II);
      uint64_t Guaranteed = A->getAlign().value();
      while (Offset % Guaranteed) Guaranteed /= 2;
      if (II->getParamAlign(Load ? 0 : 1).valueOrOne().value() > Guaranteed)
        reject("new full access alignment not established by context allocation/offset", II);
      Out.Memory.push_back(II); Out.Masked.push_back(II);
    }
  }
  if (!Out.Zero || Out.Masked.empty()) reject("full zero initialization or target masked fields absent");
  DominatorTree DT(F);
  for (Instruction *I : Out.Memory) if (I != Out.Zero && !DT.dominates(Out.Zero, I))
    reject("full context zero does not dominate every context memory access", I);
  return Out;
}
void normalize(Context &Object) {
  for (IntrinsicInst *II : Object.Masked) {
    bool Load = II->getIntrinsicID() == Intrinsic::masked_load;
    unsigned PI = Load ? 0 : 1; Value *P = II->getArgOperand(PI);
    Type *T = Load ? II->getType() : II->getArgOperand(0)->getType();
    Align Alignment = II->getParamAlign(PI).valueOrOne(); IRBuilder<> At(II);
    Value *Old = At.CreateAlignedLoad(T, P, Alignment, "closed.context.old.initialized");
    Value *Mask = II->getArgOperand(Load ? 1 : 2);
    if (Load) II->replaceAllUsesWith(At.CreateSelect(Mask, Old, II->getArgOperand(2), "closed.context.masked.read"));
    else At.CreateAlignedStore(At.CreateSelect(Mask, II->getArgOperand(0), Old, "closed.context.masked.merge"), P, Alignment);
    II->eraseFromParent();
  }
}
struct ClosedContext : PassInfoMixin<ClosedContext> {
  PreservedAnalyses run(Module &M, ModuleAnalysisManager &AM) {
    requirePinned(M);
    Function *Entry = M.getFunction(EntryName), *Masked = M.getFunction(MaskedName);
    if (!Entry || Entry->isDeclaration() || Entry->getCallingConv() != CallingConv::C ||
        Entry->arg_size() != 6 || !Entry->getReturnType()->isVoidTy() || !Masked || Masked->isDeclaration())
      reject("ordinary/masked private raw ABI differs");
    for (unsigned I = 0; I < 6; ++I) {
      Type *T = Entry->getArg(I)->getType();
      if ((I == 0 || I == 2) ? !T->isPointerTy() : !T->isIntegerTy(32)) reject("ordinary6-argument ABI differs");
    }
    unsigned MathDeclarations = 0;
    for (Function &F : M) if (F.getName().starts_with(MathPrefix)) {
      if (!F.isDeclaration() || F.getCallingConv() != CallingConv::ID(92) ||
          !F.hasFnAttribute(Attribute::StrictFP) ||
          F.getMemoryEffects() != MemoryEffects::inaccessibleMemOnly(ModRefInfo::ModRef))
        reject("math must remain external typed cc92/inaccessibleRW declarations");
      ++MathDeclarations;
    }
    if (MathDeclarations != 9) reject("nine fixed math declarations required before DCE");
    for (User *U : Masked->users()) {
      auto *Call = dyn_cast<CallInst>(U);
      if (!Call || Call->getCalledFunction() != Masked || Call->getFunction() != Entry)
        reject("masked export has an unknown caller/address consumer", dyn_cast<Instruction>(U));
    }
    Masked->setLinkage(GlobalValue::InternalLinkage); Masked->setVisibility(GlobalValue::DefaultVisibility); Masked->setDSOLocal(true);
    // Remove the unused duplicate export/body BEFORE spending work inlining it.
    ModulePassManager Pre; Pre.addPass(GlobalDCEPass()); Pre.run(M, AM);
    unsigned RawDefinitions = 0;
    for (Function &F : M) if (!F.isDeclaration() && &F != Entry) {
      if (!F.hasLocalLinkage() || F.getName().starts_with(MathPrefix) || F.hasFnAttribute(Attribute::OptimizeNone))
        reject("unknown/nonraw/no-opt definition would be forced inline");
      F.removeFnAttr(Attribute::NoInline); F.addFnAttr(Attribute::AlwaysInline); ++RawDefinitions;
    }
    for (Function &F : M) for (Instruction &I : instructions(F)) if (auto *Call = dyn_cast<CallBase>(&I))
      if (Function *Callee = Call->getCalledFunction()) if (!Callee->isDeclaration() && Callee != Entry)
        Call->removeFnAttr(Attribute::NoInline);
    ModulePassManager Inline;
    Inline.addPass(AlwaysInlinerPass()); Inline.addPass(GlobalDCEPass()); Inline.run(M, AM);
    for (Function &F : M) if (!F.isDeclaration() && &F != Entry) reject("raw definition remained after one inlining/DCE pass");
    if (verifyModule(M, &errs())) reject("LLVM rejected raw inlining output");
    // Inlining can expand static FP/math call counts. Preserve its resulting
    // inventory exactly across only the subsequent private-memory rewrite.
    Inventory Before = fpInventory(M);
    Context Object = inspect(*Entry, M.getDataLayout());
    unsigned Count = Object.Masked.size(); normalize(Object);
    if (fpInventory(M) != Before) reject("context normalization changed post-inline FP/math inventory");
    if (verifyModule(M, &errs())) reject("LLVM rejected context normalization output");
    errs() << "sirius-closed-dp-context: raw_definitions_inlined=" << RawDefinitions
           << " owned_contexts=1 context_bytes=120832 dominated_memory_accesses=" << Object.Memory.size()
           << " local_masked_operations_normalized=" << Count << " postinline_FP_inventory_unchanged=1\n";
    return PreservedAnalyses::none();
  }
};
} // namespace
extern "C" LLVM_ATTRIBUTE_WEAK PassPluginLibraryInfo llvmGetPassPluginInfo() {
  return {LLVM_PLUGIN_API_VERSION, "SiriusClosedDPContext", LLVM_VERSION_STRING,
    [](PassBuilder &PB) { PB.registerPipelineParsingCallback(
      [](StringRef Name, ModulePassManager &MPM, ArrayRef<PassBuilder::PipelineElement>) {
        if (Name != "sirius-closed-dp-context") return false; MPM.addPass(ClosedContext()); return true;
      }); }, nullptr};
}
'''


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def doc(path):
    return json.loads(path.read_text())


def terminal(path, success=True):
    row = doc(path)
    assert row['status'] == 'terminal' and row['passed'] is success
    assert row['child_exit'] == (0 if success else 1)
    assert row['source_unchanged'] and row['whole_inputs_unchanged']
    assert not row['remaining_owned_processes'] and not row['cleanup_errors']
    assert row['observed_births'] and all(birth['absent'] for birth in row['observed_births'])
    assert path.with_name('inputs-before.json').read_bytes() == path.with_name('inputs-after.json').read_bytes()


def main():
    assert sys.argv[1:] == ['derive'], 'Only fixed deferred source derivation; no force/compile mode.'
    assert pin(SOURCE_CONTRACT) == SOURCE_CONTRACT_PIN
    source = doc(SOURCE_CONTRACT)
    proof = source['source_proof']
    for key in ('source_inverse_exact', 'canonical_payload_unwrapped_only', 'canonical_tuple_records_match_exact_admitted_table',
        'canonical_root_validation_and_publication_match_table', 'all_validation_before_first_publication',
        'complete_private_context_memset_preserved', 'sole_supported_caller_required200_lane0',
        'full8_physical_payloads_and_tags_unchanged', 'generic_program_fallback_removed_from_new_private_body_only'):
        assert proof[key] is True, key
    assert (proof['layer_calls'], proof['validated_roots'], proof['published_roots']) == (45, 200, 200)
    assert pin(B / 'ispc-source/dopri_phase.ispc') == source['derived_source']
    terminal(W / 'value/owner.json', success=False)
    terminal(W / 'value-audit/owner.json')
    result = doc(W / 'value-result.json')
    assert pin(INPUT) == VALUE_PIN
    expected = dict(path=str(INPUT.relative_to(R)), **VALUE_PIN)
    assert result['output'] == expected
    assert result['original_failed_owner'] == pin(W / 'value/owner.json')
    assert result['accepted_artifact_owner'] == str((W / 'value-audit/owner.json').relative_to(R))
    assert result['lowering_repeated'] is False
    assert sum(result['pointer_calls'].values()) == 1576 and sum(result['value_calls'].values()) == 5
    audited = doc(W / 'value-audit/inputs-before.json')[str(INPUT.relative_to(R))]
    assert {key: audited[key] for key in ('bytes', 'sha256')} == VALUE_PIN
    assert Path(audited['resolved_path']).resolve() == INPUT.resolve()
    value_contract = V / 'closed_dp_value_contract.json'
    value = doc(value_contract)
    assert result['source_value_contract'] == pin(value_contract)
    assert value['phase'] == 'value' and value['source_inverse_exact']
    assert value['value_math_prefix'] == 'sirius_local_value_'
    assert pin(R / value['output']['path']) == {key: value['output'][key] for key in ('bytes', 'sha256')}
    inputs = dict(value['source_inputs'])
    inputs.update(source['source_inputs'])
    for name, identity in inputs.items():
        assert pin(R / name) == identity, name
    for path in (Path(__file__).resolve(), SOURCE_CONTRACT, value_contract, W / 'value-result.json',
                 W / 'value/owner.json', W / 'value/inputs-before.json', W / 'value/inputs-after.json',
                 W / 'value-audit/owner.json', W / 'value-audit/inputs-before.json', W / 'value-audit/inputs-after.json',
                 INPUT, B / 'ispc-source/dopri_phase.ispc'):
        inputs[str(path.relative_to(R))] = pin(path)
    text = CPP
    for token, replacement in (('@INPUT@', json.dumps(str(INPUT))), ('@SHA@', json.dumps(expected['sha256'])),
                               ('@BYTES@', str(expected['bytes']))):
        assert text.count(token) == 1
        text = text.replace(token, replacement, 1)
    assert '@INPUT@' not in text and '@SHA@' not in text and '@BYTES@' not in text
    target, contract = V / 'closed_dp_context.cpp', V / 'closed_dp_context_contract.json'
    assert not target.exists() and not contract.exists()
    with target.open('x') as stream:
        stream.write(text)
    record = dict(phase='context', source_only=True, pass_name='sirius-closed-dp-context',
        input=expected, source_inputs=inputs, output=dict(path=str(target.relative_to(R)), **pin(target)),
        private_scope=value['private_scope'],
        source_proof='Exact private DP source/table/caller/physical initialization proofs preserved; complete120832-byte owned context is zeroed before its GlobalParams store and every arithmetic access. No unproved lane, pointer or data freeze.',
        pipeline='Pre-DCE the internalized unused masked export; mark defined internal RAW helpers alwaysinline; run AlwaysInliner once and post-DCE. Math9 stays external cc92/inaccessibleRW. Then independently prove one exact initialized constant-address context and normalize only its32-byte masked vector fields.',
        initialization='Require the exact source-bound context layout, one allocation, one complete zero memset, and dominance over every context load/store/memset/masked access. Full8 physical payload/tag/status initialization and original defined masks justify load/select/store. Preserve passthroughs and alignment, including inactive lanes.',
        confinement='Offset0 pointer loads/stores retain actual GlobalParams. Stop the walk at loaded external pointers; never normalize global views, output, staged output, publication temporaries, other allocas, callbacks or unknown/dynamic context addresses. No compact IO, slot-copy interface, FP effect rewrite or math definition inlining.',
        inventory='Ordinary raw inlining may expand static call/FP counts. Compare complete constrained/FP-representation/math inventory from post-inline to post-normalization only; never require old1576 to remain after inlining.',
        next_gate='Root focused SROA/local cleanup must eliminate the owned KernelContext and all raw helper/layer escapes before native emission or costly acceptance. Actual memory/address/stack/code/ISA/noFMA evidence is required; source/inlining alone is no speed/SSA or numerical claim.',
        qualification='Source preparation only; no compile/load/pass/cleanup/native/runtime/provider/scientific/release qualification or acceptance transfer.')
    with contract.open('x') as stream:
        json.dump(record, stream, indent=2)
        stream.write('\n')
    print(json.dumps(dict(source_only=True, compiler_run=False, output=record['output'],
        contract=dict(path=str(contract.relative_to(R)), **pin(contract)), pass_name=record['pass_name'])))


if __name__ == '__main__':
    main()
