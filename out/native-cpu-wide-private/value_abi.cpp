// Replace the private pointer calls in the sealed strict SIMD16 modules with
// the integer value ABI. Context and external row staging remain in place.
#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/STLExtras.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/ADT/StringExtras.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/AsmParser/Parser.h"
#include "llvm/Config/llvm-config.h"
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
#include <array>
#include <memory>

#if LLVM_VERSION_MAJOR != 23 || LLVM_VERSION_MINOR != 1 || LLVM_VERSION_PATCH != 2
#error "The sealed strict IR and value ABI require LLVM 23.1.2."
#endif

using namespace llvm;

namespace {
constexpr unsigned Lanes = 16;
constexpr unsigned ColumnBytes = Lanes * 4;
constexpr StringLiteral PointerPrefix = "sirius_wide_strict_pointer_";
constexpr StringLiteral ValuePrefix = "sirius_wide_value_";

struct Profile {
  StringLiteral Stage, SHA;
  uint64_t Bytes;
  unsigned Calls, Tagged;
};
constexpr Profile Profiles[]{
    {"transport", "7e4f828c8183fcbab9e38d5dce790e217fc3e1b47c37b924d8630cc3d02be4d9", 122298881, 2698, 2682},
    {"endpoint", "9d53adccd9180e28b379e0de463f5afe81310702e5b966d9745c7677ee7f881e", 178541326, 3954, 3924},
    {"dopri_phase", "14540760dba2075f7febed7834d7a7978009a176f5c0ca0ecfb54c25f9f5b935", 79407257, 1576, 1575},
};
constexpr StringLiteral Recipes[]{"dense_add16", "dense_multiply16", "dense_divide16", "dense_sqrt16",
    "tagged_add16", "tagged_subtract16", "tagged_multiply16", "tagged_divide16", "tagged_sqrt16"};

[[noreturn]] void refuse(const Twine &Reason, const Instruction *I = nullptr) {
  errs() << "SIMD16 value ABI: " << Reason << '\n';
  if (I) {
    I->print(errs());
    errs() << '\n';
  }
  report_fatal_error("unadmitted SIMD16 packed bank", false);
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

const Profile &admit(const Module &M) {
  auto Input = MemoryBuffer::getFile(M.getModuleIdentifier());
  if (!Input) refuse("strict input file is required");
  StringRef Text = (*Input)->getBuffer();
  const std::string SHA = digest(Text);
  for (const Profile &P : Profiles) {
    if (Text.size() != P.Bytes || SHA != P.SHA) continue;
    LLVMContext Context;
    SMDiagnostic Diagnostic;
    auto Reference = parseAssembly((*Input)->getMemBufferRef(), Diagnostic, Context);
    if (!Reference) refuse("sealed strict input could not be parsed");
    Reference->setModuleIdentifier(M.getModuleIdentifier());
    if (fingerprint(M) != fingerprint(*Reference)) refuse("strict input was transformed before the value ABI pass");
    return P;
  }
  refuse("strict input SHA-256 or length changed");
}

bool vectorOf(Type *T, unsigned Bits) {
  auto *V = dyn_cast<FixedVectorType>(T);
  return V && V->getNumElements() == Lanes && V->getElementType()->isIntegerTy(Bits);
}

struct Bank {
  AllocaInst *Root;
  CallInst *Math;
  unsigned Columns;
  SmallVector<StoreInst *, 18> Stores;
  SmallVector<MemSetInst *, 2> Zeros;
  SmallVector<IntrinsicInst *, 12> Masked;
  SmallVector<Instruction *, 2> Starts, Ends;
};

uint64_t offset(const Bank &B, Value *Pointer, uint64_t Bytes, const DataLayout &Layout) {
  int64_t Position = 0;
  uint64_t Extent = uint64_t(B.Columns) * ColumnBytes;
  if (GetPointerBaseWithConstantOffset(Pointer, Position, Layout) != B.Root ||
      Position < 0 || uint64_t(Position) > Extent || Bytes > Extent - uint64_t(Position))
    refuse("private bank access outside its constant allocation", B.Math);
  return Position;
}

uint64_t sizeOf(Type *T, const DataLayout &Layout, const Instruction *I) {
  TypeSize Size = Layout.getTypeStoreSize(T);
  if (Size.isScalable()) refuse("scalable private access", I);
  return Size.getFixedValue();
}

Bank inspectBank(AllocaInst *Root, CallInst *Math, unsigned Columns, const DataLayout &Layout) {
  const auto *Array = dyn_cast<ArrayType>(Root->getAllocatedType());
  const auto *Count = dyn_cast<ConstantInt>(Root->getArraySize());
  if (!Array || !Count || !Count->isOne() || Root->getAddressSpace() != 0 ||
      !Array->getElementType()->isIntegerTy(32) || Array->getNumElements() != Columns * Lanes)
    refuse("wrong private bank extent", Root);
  Bank B{Root, Math, Columns, {}, {}, {}, {}, {}};
  SmallVector<Value *, 32> Pending{Root};
  SmallPtrSet<Value *, 32> Seen;
  while (!Pending.empty()) {
    Value *Pointer = Pending.pop_back_val();
    if (!Seen.insert(Pointer).second) continue;
    for (User *User : Pointer->users()) {
      if (auto *GEP = dyn_cast<GetElementPtrInst>(User)) {
        if (GEP->getPointerOperand() != Pointer) refuse("bank address used as GEP index", GEP);
        offset(B, GEP, 0, Layout);
        Pending.push_back(GEP);
        continue;
      }
      if (auto *Load = dyn_cast<LoadInst>(User)) {
        if (Load->getPointerOperand() != Pointer || !Load->isSimple()) refuse("nonordinary bank load", Load);
        offset(B, Pointer, sizeOf(Load->getType(), Layout, Load), Layout);
        continue;
      }
      if (auto *Store = dyn_cast<StoreInst>(User)) {
        if (Store->getPointerOperand() != Pointer || !Store->isSimple()) refuse("bank address stored or used atomically", Store);
        offset(B, Pointer, sizeOf(Store->getValueOperand()->getType(), Layout, Store), Layout);
        B.Stores.push_back(Store);
        continue;
      }
      if (User == Math) {
        if (Pointer != Root) refuse("pointer helper receives an interior bank address", Math);
        continue;
      }
      if (auto *Zero = dyn_cast<MemSetInst>(User)) {
        const auto *Length = dyn_cast<ConstantInt>(Zero->getLength());
        const auto *Byte = dyn_cast<ConstantInt>(Zero->getValue());
        if (Zero->getRawDest() != Pointer || Zero->isVolatile() || !Length || !Byte || !Byte->isZero())
          refuse("nonzero or unknown bank initialization", Zero);
        offset(B, Pointer, Length->getZExtValue(), Layout);
        B.Zeros.push_back(Zero);
        continue;
      }
      auto *Intrinsic = dyn_cast<IntrinsicInst>(User);
      if (!Intrinsic) refuse("private bank escapes", dyn_cast<Instruction>(User));
      switch (Intrinsic->getIntrinsicID()) {
      case Intrinsic::lifetime_start: B.Starts.push_back(Intrinsic); break;
      case Intrinsic::lifetime_end: B.Ends.push_back(Intrinsic); break;
      case Intrinsic::masked_load:
      case Intrinsic::masked_store: {
        const bool Load = Intrinsic->getIntrinsicID() == Intrinsic::masked_load;
        if (Intrinsic->arg_size() != 3 || Intrinsic->getArgOperand(Load ? 0 : 1) != Pointer ||
            !vectorOf(Load ? Intrinsic->getType() : Intrinsic->getArgOperand(0)->getType(), 32) ||
            !vectorOf(Intrinsic->getArgOperand(Load ? 1 : 2)->getType(), 1))
          refuse("unexpected masked bank operation", Intrinsic);
        offset(B, Pointer, ColumnBytes, Layout);
        B.Masked.push_back(Intrinsic);
        break;
      }
      default: refuse("unsupported intrinsic uses a private bank", Intrinsic);
      }
    }
  }
  return B;
}

bool initialized(const Bank &B, unsigned Column, Instruction *At, bool ZeroOnly,
                 const DominatorTree &DT, const DataLayout &Layout) {
  const uint64_t Start = uint64_t(Column) * ColumnBytes;
  for (MemSetInst *Zero : B.Zeros) {
    const uint64_t Length = cast<ConstantInt>(Zero->getLength())->getZExtValue();
    const uint64_t Position = offset(B, Zero->getRawDest(), Length, Layout);
    if (Position <= Start && Length >= ColumnBytes && Start - Position <= Length - ColumnBytes && DT.dominates(Zero, At))
      return true;
  }
  for (StoreInst *Store : B.Stores) {
    Value *V = Store->getValueOperand();
    if (!vectorOf(V->getType(), 32) || offset(B, Store->getPointerOperand(), ColumnBytes, Layout) != Start ||
        !DT.dominates(Store, At)) continue;
    if (!ZeroOnly || (isa<Constant>(V) && cast<Constant>(V)->isNullValue())) return true;
  }
  return false;
}

bool overlap(uint64_t A, uint64_t Bytes, uint64_t Start, uint64_t End) {
  return Bytes && A < End && Start < A + Bytes;
}

void requireLifetime(const Bank &B, Instruction *At, const DominatorTree &DT) {
  if (llvm::any_of(B.Ends, [&](Instruction *End) { return DT.dominates(End, At); }))
    refuse("bank lifetime ended before access", At);
  if (!B.Starts.empty() && !llvm::any_of(B.Starts, [&](Instruction *Start) { return DT.dominates(Start, At); }))
    refuse("bank lifetime start does not dominate access", At);
}

void proveInitialization(const Bank &B, unsigned Destination, const DominatorTree &DT, const DataLayout &Layout) {
  requireLifetime(B, B.Math, DT);
  for (unsigned Column = 0; Column != B.Columns; ++Column)
    if (!initialized(B, Column, B.Math, false, DT, Layout)) refuse("bank is not initialized before math", B.Math);
  const unsigned Columns = B.Columns == 15 ? 5 : 1;
  const uint64_t Start = uint64_t(Destination) * Columns * ColumnBytes;
  const uint64_t End = Start + Columns * ColumnBytes;
  for (unsigned Column = 0; Column != Columns; ++Column)
    if (!initialized(B, Destination * Columns + Column, B.Math, true, DT, Layout))
      refuse("destination is not initially zero", B.Math);

  for (StoreInst *Store : B.Stores) {
    Value *V = Store->getValueOperand();
    const uint64_t Width = sizeOf(V->getType(), Layout, Store);
    if (overlap(offset(B, Store->getPointerOperand(), Width, Layout), Width, Start, End) &&
        !(isa<Constant>(V) && cast<Constant>(V)->isNullValue()))
      refuse("ordinary nonzero write overlaps the zero destination", Store);
  }
  for (IntrinsicInst *I : B.Masked) {
    const bool Load = I->getIntrinsicID() == Intrinsic::masked_load;
    const unsigned PointerArg = Load ? 0 : 1;
    const uint64_t Position = offset(B, I->getArgOperand(PointerArg), ColumnBytes, Layout);
    requireLifetime(B, I, DT);
    if (Position % ColumnBytes || !initialized(B, Position / ColumnBytes, I, false, DT, Layout))
      refuse("masked access would read an uninitialized full column", I);
    if (I->getParamAlign(PointerArg).valueOrOne() > commonAlignment(B.Root->getAlign(), Position))
      refuse("full-column access would overstate its alignment", I);
    if (!Load && overlap(Position, ColumnBytes, Start, End)) refuse("masked store overlaps the zero destination", I);
  }
}

void normalize(const Bank &B) {
  for (IntrinsicInst *I : B.Masked) {
    const bool Load = I->getIntrinsicID() == Intrinsic::masked_load;
    const unsigned PointerArg = Load ? 0 : 1;
    IRBuilder<> Builder(I);
    Value *Pointer = I->getArgOperand(PointerArg);
    Type *Type = Load ? I->getType() : I->getArgOperand(0)->getType();
    Align Alignment = I->getParamAlign(PointerArg).valueOrOne();
    Value *Old = Builder.CreateAlignedLoad(Type, Pointer, Alignment, "wide.bank.old");
    Value *Mask = I->getArgOperand(Load ? 1 : 2);
    if (Load) {
      I->replaceAllUsesWith(Builder.CreateSelect(Mask, Old, I->getArgOperand(2), "wide.bank.loaded"));
    } else {
      Value *Selected = Builder.CreateSelect(Mask, I->getArgOperand(0), Old, "wide.bank.merged");
      Builder.CreateAlignedStore(Selected, Pointer, Alignment);
    }
    I->eraseFromParent();
  }
}

const MemoryEffects MathEffects = MemoryEffects::inaccessibleMemOnly();

Function *declareCore(Module &M, StringRef Recipe, bool Tagged) {
  LLVMContext &C = M.getContext();
  Type *Vector = FixedVectorType::get(Type::getInt32Ty(C), Lanes);
  SmallVector<Type *, 13> Args(Tagged ? 12 : 10, Vector);
  Args.push_back(Type::getInt32Ty(C));
  SmallVector<Type *, 6> Results(Tagged ? 6 : 5, Vector);
  std::string Name = (ValuePrefix + Recipe).str();
  if (M.getNamedValue(Name)) refuse("value core name collision");
  auto *Core = Function::Create(FunctionType::get(StructType::get(C, Results), Args, false),
                                GlobalValue::ExternalLinkage, Name, M);
  Core->setCallingConv(CallingConv::X86_RegCall);
  Core->setVisibility(GlobalValue::HiddenVisibility);
  Core->setDSOLocal(true);
  Core->setMemoryEffects(MathEffects);
  Core->addFnAttr(Attribute::StrictFP);
  Core->addFnAttr(Attribute::NoUnwind);
  Core->addFnAttr(Attribute::NoInline);
  return Core;
}

Value *slot(IRBuilder<> &B, Value *Pointer, unsigned Column) {
  return B.CreateConstGEP1_64(B.getInt8Ty(), Pointer, uint64_t(Column) * ColumnBytes, "wide.bank.slot");
}

void bindMath(CallInst &Old, Function &Core, bool Tagged) {
  IRBuilder<> B(&Old);
  Type *V = FixedVectorType::get(B.getInt32Ty(), Lanes);
  const unsigned A = Tagged ? 1 : 0;
  const unsigned Destination = Tagged ? 0 : 2;
  Value *Payload = Old.getArgOperand(0);
  SmallVector<Value *, 13> Args;
  for (unsigned Operand : {A, A + 1})
    for (unsigned Column = 0; Column != 5; ++Column)
      Args.push_back(B.CreateAlignedLoad(V, slot(B, Payload, Operand * 5 + Column), Align(4)));
  if (Tagged) {
    Value *Tags = Old.getArgOperand(1);
    for (unsigned Operand : {1U, 2U})
      Args.push_back(B.CreateAlignedLoad(V, slot(B, Tags, Operand), Align(4)));
  }
  Args.push_back(Old.getArgOperand(Old.arg_size() - 1));
  CallInst *Call = B.CreateCall(&Core, Args, "wide.value");
  Call->setCallingConv(CallingConv::X86_RegCall);
  Call->setMemoryEffects(MathEffects);
  Call->addFnAttr(Attribute::StrictFP);
  Call->setDoesNotThrow();
  for (unsigned Column = 0; Column != 5; ++Column)
    B.CreateAlignedStore(B.CreateExtractValue(Call, Column), slot(B, Payload, Destination * 5 + Column), Align(4));
  if (Tagged) B.CreateAlignedStore(B.CreateExtractValue(Call, 5), Old.getArgOperand(1), Align(4));
  Old.eraseFromParent();
}

SmallVector<const Instruction *> strictOperations(const Module &M) {
  SmallVector<const Instruction *> Result;
  for (const Function &F : M) {
    for (const Instruction &I : instructions(F)) {
      if (isa<ConstrainedFPIntrinsic>(I) || I.getOpcode() == Instruction::FNeg) Result.push_back(&I);
      if (const auto *Call = dyn_cast<CallInst>(&I)) {
        const Function *Target = Call->getCalledFunction();
        if (Target && (Target->getIntrinsicID() == Intrinsic::fabs || Target->getName() == "sirius_simd_fmaxf"))
          Result.push_back(&I);
      }
    }
  }
  return Result;
}

struct ValueABI : PassInfoMixin<ValueABI> {
  PreservedAnalyses run(Module &M, ModuleAnalysisManager &) {
    const Profile &P = admit(M);
    const DataLayout &Layout = M.getDataLayout();
    const auto FP = strictOperations(M);
    SmallVector<Function *, 9> Helpers;
    SmallVector<Bank, 8> Banks;
    SmallPtrSet<AllocaInst *, 32> Roots;
    DenseMap<Function *, std::unique_ptr<DominatorTree>> Dominators;
    unsigned Calls = 0, TaggedCalls = 0, Masked = 0;
    for (Function &F : M) {
      if (!F.getName().starts_with(PointerPrefix)) continue;
      const StringRef Recipe = F.getName().drop_front(PointerPrefix.size());
      const bool Tagged = Recipe.starts_with("tagged_");
      const unsigned Pointers = Tagged ? 2 : 1;
      if (!llvm::is_contained(Recipes, Recipe) || !F.isDeclaration() || !F.doesNotThrow() ||
          F.getVisibility() != GlobalValue::HiddenVisibility || F.arg_size() != Pointers + 4 || !F.getReturnType()->isVoidTy())
        refuse("unexpected sealed pointer declaration");
      for (User *U : F.users()) {
        auto *Call = dyn_cast<CallInst>(U);
        if (!Call || Call->getCalledFunction() != &F || !Call->hasFnAttr(Attribute::StrictFP))
          refuse("pointer helper has an indirect or unscoped use");
        Function *Worker = Call->getFunction();
        auto &DT = Dominators[Worker];
        if (!DT) DT = std::make_unique<DominatorTree>(*Worker);
        for (unsigned Arg = 0; Arg != Pointers; ++Arg) {
          auto *Root = dyn_cast<AllocaInst>(Call->getArgOperand(Arg));
          if (!Root || !Roots.insert(Root).second) refuse("bank reused across helpers", Call);
          Banks.push_back(inspectBank(Root, Call, Arg == 0 ? 15 : 3, Layout));
          proveInitialization(Banks.back(), Tagged ? 0 : 2, *DT, Layout);
          Masked += Banks.back().Masked.size();
        }
        ++Calls;
        TaggedCalls += Tagged;
      }
      Helpers.push_back(&F);
    }
    if (Calls != P.Calls || TaggedCalls != P.Tagged) refuse("sealed helper or tag inventory changed");
    for (const Bank &B : Banks) normalize(B);
    for (Function *Helper : Helpers) {
      const StringRef Recipe = Helper->getName().drop_front(PointerPrefix.size());
      const bool Tagged = Recipe.starts_with("tagged_");
      Function *Core = declareCore(M, Recipe, Tagged);
      while (!Helper->use_empty()) bindMath(*cast<CallInst>(*Helper->user_begin()), *Core, Tagged);
      Helper->eraseFromParent();
    }
    if (strictOperations(M) != FP || verifyModule(M, &errs())) refuse("value binding altered the strict arithmetic or invalidated IR");
    errs() << "SIMD16 value ABI: " << P.Stage << " calls=" << Calls << " tags=" << TaggedCalls
           << " masked=" << Masked << '\n';
    return PreservedAnalyses::none();
  }
};
} // namespace

extern "C" LLVM_ATTRIBUTE_WEAK PassPluginLibraryInfo llvmGetPassPluginInfo() {
  return {LLVM_PLUGIN_API_VERSION, "SiriusWideValueABI", LLVM_VERSION_STRING,
          [](PassBuilder &Builder) {
            Builder.registerPipelineParsingCallback(
                [](StringRef Name, ModulePassManager &Passes, ArrayRef<PassBuilder::PipelineElement> Inner) {
                  if (Name != "sirius-wide-value" || !Inner.empty()) return false;
                  Passes.addPass(ValueABI());
                  return true;
                });
          }, nullptr};
}
