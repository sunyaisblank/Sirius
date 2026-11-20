// Normalize only the owned SIMD16 kernel context. The transactional output
// allocation, external views, and strict value math have separate lifetimes.
#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/MapVector.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/ADT/StringExtras.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/AsmParser/Parser.h"
#include "llvm/Config/llvm-config.h"
#include "llvm/IR/CFG.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Verifier.h"
#include "llvm/Passes/PassBuilder.h"
#include "llvm/Plugins/PassPlugin.h"
#include "llvm/Support/ErrorHandling.h"
#include "llvm/Support/MemoryBuffer.h"
#include "llvm/Support/ModRef.h"
#include "llvm/Support/SHA256.h"
#include "llvm/Support/SourceMgr.h"
#include "llvm/Transforms/IPO/AlwaysInliner.h"
#include "llvm/Transforms/IPO/GlobalDCE.h"
#include <optional>

#if LLVM_VERSION_MAJOR != 23 || LLVM_VERSION_MINOR != 1 || LLVM_VERSION_PATCH != 2
#error "Context lowering requires the pinned LLVM23.1.2 API."
#endif

using namespace llvm;

namespace {
constexpr unsigned Lanes = 16;
constexpr unsigned ColumnBytes = Lanes * 4;
constexpr StringLiteral Prefix = "sirius_wide_strict_";
constexpr StringLiteral MaskedSuffix = "___un_3C_unu_3E_unuun_3C_unu_3E_unuunuunu";
constexpr StringLiteral MathPrefix = "sirius_wide_value_";

struct Profile {
  StringLiteral Entry, SHA;
  uint64_t InputBytes, ContextBytes;
  unsigned OutputWords, MathDeclarations;
  ArrayRef<unsigned> Columns;
};
constexpr unsigned TransportColumns[]{2295, 459, 1, 125};
constexpr unsigned EndpointColumns[]{7, 1, 3060, 612};
constexpr unsigned DopriColumns[]{1, 3145, 629};
constexpr Profile Transport{"transport", "7f07c1a9a34111086b8ba61c0c9a0a2a2886a382353401dc0c9b189e345ec24d",
    131087821, 184384, 4699, 8, TransportColumns};
constexpr Profile Endpoint{"endpoint", "98dcdd9b28c9e669e9c2ebd3d991b36228412ca361793eb4aebcb79b1e7b63af",
    191430599, 235584, 3829, 9, EndpointColumns};
constexpr Profile Dopri{"dp", "e38d743e768a66429d01623e544961068ad086ea9837a83aaf92460585a47fa0",
    84578608, 241664, 4149, 5, DopriColumns};

[[noreturn]] void reject(const Twine &Why, const Instruction *I = nullptr) {
  errs() << "SIMD16 context: " << Why << '\n';
  if (I) {
    I->print(errs());
    errs() << '\n';
  }
  report_fatal_error("unadmitted SIMD16 context", false);
}

std::string digest(StringRef Text) {
  return toHex(SHA256::hash(arrayRefFromStringRef(Text)), true);
}

std::string fingerprint(const Module &M) {
  std::string Text;
  raw_string_ostream Stream(Text);
  M.print(Stream, nullptr);
  return digest(Stream.str());
}

bool valueIdentity(StringRef Text, StringRef InputPath, const Profile &P) {
  auto [Header, Body] = Text.split('\n');
  StringRef Stage = &P == &Dopri ? StringRef("dopri_phase") : StringRef(P.Entry);
  StringRef Directory = InputPath.take_front(InputPath.rfind('/') + 1);
  const auto Current = (Twine("; ModuleID = '") + Directory + Stage + ".strict.ll'").str();
  const auto Sealed = (Twine("; ModuleID = '/home/astra/.project/Sirius/out/simd16-value-azgcxqtz/") +
                       Stage + ".strict.ll'").str();
  if ((Header != Current && Header != Sealed) || Sealed.size() + 1 + Body.size() != P.InputBytes)
    return false;
  // opt writes its input filename in this comment. Check every other byte
  // against the sealed value module when replay uses a new task directory.
  SHA256 Hash;
  Hash.update(Sealed);
  Hash.update("\n");
  Hash.update(Body);
  return toHex(Hash.final(), true) == P.SHA;
}

const Profile &admit(const Module &M) {
  auto Input = MemoryBuffer::getFile(M.getModuleIdentifier());
  if (!Input) reject("sealed value input is required");
  StringRef Text = (*Input)->getBuffer();
  for (const Profile *P : {&Transport, &Endpoint, &Dopri}) {
    if (!valueIdentity(Text, M.getModuleIdentifier(), *P)) continue;
    LLVMContext Context;
    SMDiagnostic Diagnostic;
    auto Reference = parseAssembly((*Input)->getMemBufferRef(), Diagnostic, Context);
    if (!Reference) reject("sealed value input could not be parsed");
    Reference->setModuleIdentifier(M.getModuleIdentifier());
    if (fingerprint(M) != fingerprint(*Reference)) reject("value input transformed before context binding");
    return *P;
  }
  reject("sealed value input SHA-256 or length differs");
}

bool words(Type *T) {
  auto *V = dyn_cast<FixedVectorType>(T);
  return V && V->getNumElements() == Lanes && V->getElementType()->isIntegerTy(32);
}

bool contextType(Type *T, const DataLayout &DL, const Profile &P) {
  auto *S = dyn_cast<StructType>(T);
  if (!S || S->isOpaque() || S->isPacked() || S->getNumElements() != P.Columns.size() + 1 ||
      !S->getElementType(0)->isPointerTy()) return false;
  const auto *Layout = DL.getStructLayout(S);
  uint64_t Offset = ColumnBytes;
  for (unsigned K = 0; K != P.Columns.size(); ++K) {
    unsigned Count = P.Columns[K];
    Type *Field = S->getElementType(K + 1);
    if (Count != 1) {
      auto *Array = dyn_cast<ArrayType>(Field);
      if (!Array || Array->getNumElements() != Count) return false;
      Field = Array->getElementType();
    }
    if (!words(Field) || Layout->getElementOffset(K + 1) != Offset) return false;
    Offset += uint64_t(Count) * ColumnBytes;
  }
  return Offset == P.ContextBytes && DL.getTypeAllocSize(S).getFixedValue() == P.ContextBytes;
}

AllocaInst *contextAllocation(Function &F, const DataLayout &DL, const Profile &P) {
  AllocaInst *Found = nullptr;
  for (Instruction &I : instructions(F)) {
    auto *A = dyn_cast<AllocaInst>(&I);
    if (!A || !contextType(A->getAllocatedType(), DL, P)) continue;
    auto *Count = dyn_cast<ConstantInt>(A->getArraySize());
    if (Found || !Count || !Count->isOne() || A->getAddressSpace() != 0 || A->getAlign() < Align(ColumnBytes))
      reject("unexpected, duplicate, or underaligned context allocation", A);
    Found = A;
  }
  if (!Found) reject("exact owned context allocation absent");
  return Found;
}

// Derived from the admitted SIMD8 Transport reset and coefficient proofs.
#include "transport_context.inc"

uint64_t sizeOf(Type *T, const DataLayout &DL, const Instruction *I) {
  TypeSize S = DL.getTypeStoreSize(T);
  if (S.isScalable()) reject("scalable context access", I);
  return S.getFixedValue();
}

void bounds(uint64_t Offset, uint64_t Bytes, const Profile &P, const Instruction *I) {
  if (Offset > P.ContextBytes || Bytes > P.ContextBytes - Offset)
    reject("context access outside its allocation", I);
}

struct Context {
  MemSetInst *Zero = nullptr;
  SmallVector<Instruction *, 128> Memory;
  SmallVector<IntrinsicInst *, 128> Masked;
};

Context inspect(Function &F, AllocaInst *A, const Profile &P, const DataLayout &DL) {
  Context Out;
  SmallVector<Value *, 128> Pending{A};
  SmallPtrSet<Value *, 32> Seen;
  while (!Pending.empty()) {
    Value *Pointer = Pending.pop_back_val();
    if (!Seen.insert(Pointer).second) continue;
    int64_t SignedOffset = 0;
    if (GetPointerBaseWithConstantOffset(Pointer, SignedOffset, DL) != A || SignedOffset < 0)
      reject("dynamic or nonlocal context address", dyn_cast<Instruction>(Pointer));
    uint64_t Offset = SignedOffset;
    bounds(Offset, 0, P, dyn_cast<Instruction>(Pointer));
    for (User *U : Pointer->users()) {
      if (auto *G = dyn_cast<GetElementPtrInst>(U)) {
        if (G->getPointerOperand() != Pointer) reject("context address used as a GEP index", G);
        Pending.push_back(G);
        continue;
      }
      if (auto *Load = dyn_cast<LoadInst>(U)) {
        if (Load->getPointerOperand() != Pointer || !Load->isSimple()) reject("nonordinary context load", Load);
        uint64_t Width = sizeOf(Load->getType(), DL, Load);
        bounds(Offset, Width, P, Load);
        if (Load->getType()->isPointerTy() && (Offset || Width != 8)) reject("unexpected context view field", Load);
        // The loaded view itself leads to external storage. Do not walk it.
        Out.Memory.push_back(Load);
        continue;
      }
      if (auto *Store = dyn_cast<StoreInst>(U)) {
        if (Store->getPointerOperand() != Pointer || Store->getValueOperand() == Pointer || !Store->isSimple())
          reject("context address escape or nonordinary store", Store);
        Type *T = Store->getValueOperand()->getType();
        uint64_t Width = sizeOf(T, DL, Store);
        bounds(Offset, Width, P, Store);
        if (T->isPointerTy() && (Offset || Width != 8)) reject("unexpected stored view field", Store);
        Out.Memory.push_back(Store);
        continue;
      }
      if (auto *Zero = dyn_cast<MemSetInst>(U)) {
        auto *Length = dyn_cast<ConstantInt>(Zero->getLength());
        auto *Value = dyn_cast<ConstantInt>(Zero->getValue());
        if (!Length || !Value || !Value->isZero() || Zero->isVolatile() || Zero->getRawDest() != Pointer)
          reject("unknown context initialization", Zero);
        bounds(Offset, Length->getZExtValue(), P, Zero);
        if (!Offset && Length->getZExtValue() == P.ContextBytes) {
          if (Out.Zero && Out.Zero != Zero) reject("more than one full context initialization", Zero);
          Out.Zero = Zero;
        }
        Out.Memory.push_back(Zero);
        continue;
      }
      auto *I = dyn_cast<IntrinsicInst>(U);
      if (!I) reject("context address escapes", dyn_cast<Instruction>(U));
      if ((I->getIntrinsicID() == Intrinsic::lifetime_start || I->getIntrinsicID() == Intrinsic::lifetime_end) &&
          I->arg_size() == 1 && I->getArgOperand(0) == Pointer && !Offset) continue;
      const bool Load = I->getIntrinsicID() == Intrinsic::masked_load;
      const unsigned Arg = Load ? 0 : 1;
      if ((!Load && I->getIntrinsicID() != Intrinsic::masked_store) || I->arg_size() != 3 || I->getArgOperand(Arg) != Pointer)
        reject("unexpected context intrinsic", I);
      auto *Vector = dyn_cast<FixedVectorType>(Load ? I->getType() : I->getArgOperand(0)->getType());
      auto *Mask = dyn_cast<FixedVectorType>(I->getArgOperand(Load ? 1 : 2)->getType());
      if (!Vector || Vector->getNumElements() != Lanes ||
          (!Vector->getElementType()->isIntegerTy(32) && !Vector->getElementType()->isFloatTy()) ||
          !Mask || Mask->getNumElements() != Lanes || !Mask->getElementType()->isIntegerTy(1) ||
          Offset < ColumnBytes || Offset % ColumnBytes || sizeOf(Vector, DL, I) != ColumnBytes ||
          I->getParamAlign(Arg).valueOrOne() > commonAlignment(A->getAlign(), Offset))
        reject("masked context access cannot be widened to an aligned exact column", I);
      bounds(Offset, ColumnBytes, P, I);
      Out.Memory.push_back(I);
      Out.Masked.push_back(I);
    }
  }
  if (!Out.Zero || Out.Masked.empty()) reject("expected initialized context and masked accesses are absent");
  DominatorTree Dominators(F);
  for (Instruction *I : Out.Memory)
    if (I != Out.Zero && !Dominators.dominates(Out.Zero, I)) reject("full zero fails to dominate a context access", I);
  return Out;
}

void normalize(const Context &Object) {
  for (IntrinsicInst *I : Object.Masked) {
    const bool Load = I->getIntrinsicID() == Intrinsic::masked_load;
    const unsigned Arg = Load ? 0 : 1;
    Value *Pointer = I->getArgOperand(Arg);
    Type *T = Load ? I->getType() : I->getArgOperand(0)->getType();
    Align Alignment = I->getParamAlign(Arg).valueOrOne();
    IRBuilder<> Builder(I);
    Value *Old = Builder.CreateAlignedLoad(T, Pointer, Alignment, "wide.context.old");
    Value *Mask = I->getArgOperand(Load ? 1 : 2);
    if (Load) I->replaceAllUsesWith(Builder.CreateSelect(Mask, Old, I->getArgOperand(2), "wide.context.read"));
    else Builder.CreateAlignedStore(Builder.CreateSelect(Mask, I->getArgOperand(0), Old, "wide.context.merge"), Pointer, Alignment);
    I->eraseFromParent();
  }
}

SmallVector<const Instruction *> strictOperations(const Function &F) {
  SmallVector<const Instruction *> Inventory;
  for (const Instruction &I : instructions(F)) {
    if (isa<ConstrainedFPIntrinsic>(I) || I.getOpcode() == Instruction::FNeg) Inventory.push_back(&I);
    if (auto *Call = dyn_cast<CallInst>(&I)) {
      auto *Target = Call->getCalledFunction();
      if (Target && (Target->getName().starts_with(MathPrefix) || Target->getIntrinsicID() == Intrinsic::fabs ||
          Target->getName() == "sirius_simd_fmaxf")) Inventory.push_back(&I);
    }
  }
  return Inventory;
}

struct LowerContext : PassInfoMixin<LowerContext> {
  PreservedAnalyses run(Module &M, ModuleAnalysisManager &AM) {
    const Profile &P = admit(M);
    const auto Name = (Prefix + P.Entry + "_packet").str();
    Function *Entry = M.getFunction(Name), *Masked = M.getFunction(Name + MaskedSuffix.str());
    if (!Entry || Entry->isDeclaration() || Entry->getCallingConv() != CallingConv::C || Entry->arg_size() != 6 ||
        !Entry->getReturnType()->isVoidTy() || !Masked || Masked->isDeclaration()) reject("private packet ABI differs");
    for (unsigned I = 0; I != 6; ++I) {
      Type *T = Entry->getArg(I)->getType();
      if ((I == 0 || I == 2) ? !T->isPointerTy() : !T->isIntegerTy(32)) reject("six-argument packet ABI differs");
    }
    unsigned Helpers = 0;
    for (Function &F : M) {
      if (!F.getName().starts_with(MathPrefix)) continue;
      if (!F.isDeclaration() || F.getCallingConv() != CallingConv::X86_RegCall ||
          !F.hasFnAttribute(Attribute::StrictFP) || F.getMemoryEffects() != MemoryEffects::inaccessibleMemOnly())
        reject("value helper ABI changed before context lowering");
      ++Helpers;
    }
    if (Helpers != P.MathDeclarations) reject("number of stage value helper declarations differs");
    if (&P == &Transport) markCoefficientDomains(M);
    for (User *U : Masked->users()) {
      auto *Call = dyn_cast<CallInst>(U);
      if (!Call || Call->getCalledFunction() != Masked || Call->getFunction() != Entry)
        reject("masked export has an external or indirect consumer", dyn_cast<Instruction>(U));
    }
    Masked->setLinkage(GlobalValue::InternalLinkage);
    Masked->setVisibility(GlobalValue::DefaultVisibility);
    Masked->setDSOLocal(true);
    Entry->setUWTableKind(UWTableKind::Async);
    ModulePassManager Dead;
    Dead.addPass(GlobalDCEPass());
    Dead.run(M, AM);
    for (Function &F : M) {
      if (F.isDeclaration() || &F == Entry) continue;
      if (!F.hasLocalLinkage() || F.hasFnAttribute(Attribute::OptimizeNone)) reject("non-private or no-opt inlining target");
      F.removeFnAttr(Attribute::NoInline);
      F.addFnAttr(Attribute::AlwaysInline);
    }
    for (Function &F : M) for (Instruction &I : instructions(F)) {
      if (auto *Call = dyn_cast<CallBase>(&I))
        if (auto *Callee = Call->getCalledFunction(); Callee && !Callee->isDeclaration() && Callee != Entry)
          Call->removeFnAttr(Attribute::NoInline);
    }
    ModulePassManager Inline;
    Inline.addPass(AlwaysInlinerPass());
    Inline.addPass(GlobalDCEPass());
    Inline.run(M, AM);
    for (Function &F : M)
      if (!F.isDeclaration() && &F != Entry) reject("private definition remained after inlining and DCE");
    if (verifyModule(M, &errs())) reject("LLVM rejected inlined value module");

    unsigned Staging = 0;
    for (Instruction &I : instructions(*Entry)) {
      auto *A = dyn_cast<AllocaInst>(&I);
      if (!A || A->getName() != "staged_output") continue;
      auto *T = dyn_cast<ArrayType>(A->getAllocatedType());
      auto *Count = dyn_cast<ConstantInt>(A->getArraySize());
      if (!T || T->getNumElements() != P.OutputWords * Lanes || !T->getElementType()->isIntegerTy(32) ||
          A->getAlign() != Align(ColumnBytes) || !Count || !Count->isOne()) reject("transactional staging allocation differs", A);
      ++Staging;
    }
    if (Staging != 1) reject("one exact transactional staging allocation required");
    const auto Before = strictOperations(*Entry);
    const DataLayout &DL = M.getDataLayout();
    AllocaInst *A = contextAllocation(*Entry, DL, P);
    if (&P == &Transport) {
      if (expandReset(*Entry, A, DL) != 918) reject("Transport reset expansion differs");
      isolateCoefficient(*Entry, A, DL);
    }
    Context Object = inspect(*Entry, A, P, DL);
    normalize(Object);
    if (strictOperations(*Entry) != Before || verifyModule(M, &errs()))
      reject("context lowering changed strict FP/math or invalidated LLVM IR");
    errs() << "SIMD16 context: " << P.Entry << " context_bytes=" << P.ContextBytes
           << " staging_words=" << P.OutputWords * Lanes << " bounded_accesses=" << Object.Memory.size()
           << " normalized=" << Object.Masked.size() << " strict_operations=" << Before.size() << '\n';
    return PreservedAnalyses::none();
  }
};
} // namespace

extern "C" LLVM_ATTRIBUTE_WEAK PassPluginLibraryInfo llvmGetPassPluginInfo() {
  return {LLVM_PLUGIN_API_VERSION, "SiriusWideContext", LLVM_VERSION_STRING, [](PassBuilder &Builder) {
    Builder.registerPipelineParsingCallback(
      [](StringRef Name, ModulePassManager &Passes, ArrayRef<PassBuilder::PipelineElement> Inner) {
        if (Name != "sirius-wide-context" || !Inner.empty()) return false;
        Passes.addPass(LowerContext());
        return true;
      });
  }, nullptr};
}
