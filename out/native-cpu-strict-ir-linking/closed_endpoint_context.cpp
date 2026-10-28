// Generated for one exact closed endpoint value module. No math definition is linked.
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
constexpr StringLiteral Input = "/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-canonical-closed-endpoint/endpoint.value.bc";
constexpr StringLiteral InputSHA = "61321c43b452dcc98d0c3ed1a45623012437f2aedf6c47d3e8ebe799b3d02aaa";
constexpr uint64_t InputBytes = 22601688, ContextBytes = 117792;
constexpr StringLiteral EntryName = "sirius_closed_endpoint_packet";
constexpr StringLiteral MaskedName = "sirius_closed_endpoint_packet___un_3C_unu_3E_unuun_3C_unu_3E_unuunuunu";
constexpr StringLiteral MathPrefix = "sirius_local_value_";

[[noreturn]] void reject(const Twine &Why, const Instruction *I = nullptr) {
  errs() << "sirius-closed-endpoint-context: " << Why << '\n';
  if (I) { errs() << I->getFunction()->getName() << ": "; I->print(errs()); errs() << '\n'; }
  report_fatal_error("closed endpoint context normalization refused", false);
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
  errs() << "sirius-closed-endpoint-context: input-sha256=" << InputSHA << " module-fingerprint=" << Expected << '\n';
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
  if (!S || S->isOpaque() || S->isPacked() || S->getNumElements() != 5 ||
      !S->getElementType(0)->isPointerTy() || !words(S->getElementType(2))) return false;
  for (auto Pair : {std::pair<unsigned,unsigned>{1,7}, {3,3060}, {4,612}}) {
    auto *A = dyn_cast<ArrayType>(S->getElementType(Pair.first));
    if (!A || A->getNumElements() != Pair.second || !words(A->getElementType())) return false;
  }
  auto *Layout = DL.getStructLayout(S);
  return DL.getTypeAllocSize(S).getFixedValue() == ContextBytes &&
    Layout->getElementOffset(0) == 0 &&
    Layout->getElementOffset(1) == 32 &&
    Layout->getElementOffset(2) == 256 &&
    Layout->getElementOffset(3) == 288 &&
    Layout->getElementOffset(4) == 98208;
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
    Entry->setUWTableKind(UWTableKind::Async);
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
    unsigned Staging = 0;
    for (Instruction &I : instructions(Entry)) if (auto *A = dyn_cast<AllocaInst>(&I))
      if (A->getName() == "staged_output") {
        auto *T = dyn_cast<ArrayType>(A->getAllocatedType());
        if (!T || T->getNumElements() != 30632 ||
            !T->getElementType()->isIntegerTy(32) || A->getAlign() != Align(64) ||
            !isa<ConstantInt>(A->getArraySize()) || !cast<ConstantInt>(A->getArraySize())->isOne())
          reject("original transactional staging allocation differs", A);
        ++Staging;
      }
    if (Staging != 1) reject("one exact original transactional staging allocation required");
    Inventory Before = fpInventory(M);
    Context Object = inspect(*Entry, M.getDataLayout());
    unsigned Count = Object.Masked.size(); normalize(Object);
    if (fpInventory(M) != Before) reject("context normalization changed post-inline FP/math inventory");
    if (verifyModule(M, &errs())) reject("LLVM rejected context normalization output");
    errs() << "sirius-closed-endpoint-context: raw_definitions_inlined=" << RawDefinitions
           << " owned_contexts=1 context_bytes=117792 dominated_memory_accesses=" << Object.Memory.size()
           << " local_masked_operations_normalized=" << Count << " postinline_FP_inventory_unchanged=1\n";
    return PreservedAnalyses::none();
  }
};
} // namespace
extern "C" LLVM_ATTRIBUTE_WEAK PassPluginLibraryInfo llvmGetPassPluginInfo() {
  return {LLVM_PLUGIN_API_VERSION, "SiriusClosedEndpointContext", LLVM_VERSION_STRING,
    [](PassBuilder &PB) { PB.registerPipelineParsingCallback(
      [](StringRef Name, ModulePassManager &MPM, ArrayRef<PassBuilder::PipelineElement>) {
        if (Name != "sirius-closed-endpoint-context") return false; MPM.addPass(ClosedContext()); return true;
      }); }, nullptr};
}
