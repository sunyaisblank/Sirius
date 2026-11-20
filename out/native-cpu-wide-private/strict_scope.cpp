// Convert only the three sealed, ordinary ISPC SIMD16 modules. The private
// worker still needs the value ABI, context and same-thread FP scope gates.
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/ADT/StringExtras.h"
#include "llvm/ADT/StringMap.h"
#include "llvm/AsmParser/Parser.h"
#include "llvm/Config/llvm-config.h"
#include "llvm/IR/FPEnv.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Operator.h"
#include "llvm/IR/Verifier.h"
#include "llvm/Passes/PassBuilder.h"
#include "llvm/Plugins/PassPlugin.h"
#include "llvm/Support/ErrorHandling.h"
#include "llvm/Support/MemoryBuffer.h"
#include "llvm/Support/ModRef.h"
#include "llvm/Support/SHA256.h"
#include "llvm/Support/SourceMgr.h"
#include <array>

#if LLVM_VERSION_MAJOR != 23 || LLVM_VERSION_MINOR != 1 || LLVM_VERSION_PATCH != 2
#error "The reviewed partition and native lowering require LLVM 23.1.2."
#endif

using namespace llvm;

namespace {
constexpr unsigned Lanes = 16;
constexpr StringLiteral Prefix = "sirius_wide_";
constexpr StringLiteral MathPrefix = "sirius_wide_pointer_";
constexpr StringLiteral MaskedSuffix = "___un_3C_unu_3E_unuun_3C_unu_3E_unuunuunu";

struct Stage {
  StringLiteral Filename, Entry, SHA;
  uint64_t Bytes;
  unsigned Definitions;
  // fadd, fsub, fmul, fcmp, fneg, sitofp
  std::array<unsigned, 6> Operations;
  // Quiet: oeq, une, uno. Signaling: oge, ole, ogt, olt, ugt.
  std::array<unsigned, 8> Predicates;
  unsigned PointerCalls, TaggedCalls;
};

constexpr Stage Transport{"transport", "transport", "b1cd845390f15b4dd8acdcfa930fc5f647bd1f8c7c97b149dc011693f8c06773",
                          121732773, 157, {79, 13, 0, 362, 7, 6}, {304, 12, 10, 22, 2, 12, 0, 0}, 2698, 2682};
constexpr Stage Endpoint{"endpoint", "endpoint", "1031d1ec4c86444c78d00419c9d89bad070cd6a6c1afe2c8f0b7f1db53588f37",
                         177735224, 190, {96, 13, 1, 471, 16, 0}, {418, 14, 11, 8, 3, 14, 3, 0}, 3954, 3924};
constexpr Stage Dopri{"dopri_phase", "dp", "65414cdf7d2343cc88e818457921c26d087ec1325cbfd384ac9f3474f0f4f4dd",
                      78952924, 49, {9, 16, 0, 666, 2, 0}, {656, 0, 0, 4, 5, 0, 0, 1}, 1576, 1575};

[[noreturn]] void reject(const Twine &Why, const Instruction *I = nullptr) {
  errs() << "SIMD16 strict scope: " << Why << '\n';
  if (I) {
    I->print(errs());
    errs() << '\n';
  }
  report_fatal_error("unadmitted SIMD16 IR", false);
}

std::string digest(StringRef Text) { return toHex(SHA256::hash(arrayRefFromStringRef(Text)), true); }

std::string fingerprint(const Module &M) {
  std::string Text;
  raw_string_ostream Stream(Text);
  M.print(Stream, nullptr);
  return digest(Stream.str());
}

const Stage &admit(const Module &M) {
  auto Raw = MemoryBuffer::getFile(M.getModuleIdentifier());
  if (!Raw) reject("original module file is required");
  StringRef Data = (*Raw)->getBuffer();
  std::string SHA = digest(Data);
  for (const Stage *S : {&Transport, &Endpoint, &Dopri}) {
    if (Data.size() != S->Bytes || SHA != S->SHA) continue;
    StringRef Path = M.getModuleIdentifier();
    auto Directory = Path.take_front(Path.rfind('/') + 1);
    if (Path != (Directory + S->Filename + ".ll").str() ||
        M.getSourceFileName() != (Directory + "ispc-source/" + S->Filename + ".ispc").str())
      reject("original raw/source path changed");
    LLVMContext Context;
    SMDiagnostic Diagnostic;
    auto Reference = parseAssembly((*Raw)->getMemBufferRef(), Diagnostic, Context);
    if (!Reference) reject("pinned reference could not be parsed");
    Reference->setModuleIdentifier(M.getModuleIdentifier());
    if (fingerprint(M) != fingerprint(*Reference)) reject("parsed input differs from the pinned raw module");
    return *S;
  }
  reject("ordinary IR SHA-256 or length differs from the sealed products");
}

bool wide(Type *T, bool Float) {
  auto *V = dyn_cast<FixedVectorType>(T);
  return V && V->getNumElements() == Lanes &&
         (Float ? V->getElementType()->isFloatTy() : V->getElementType()->isIntegerTy(32));
}

unsigned predicate(CmpInst::Predicate P) {
  switch (P) {
  case CmpInst::FCMP_OEQ: return 0;
  case CmpInst::FCMP_UNE: return 1;
  case CmpInst::FCMP_UNO: return 2;
  case CmpInst::FCMP_OGE: return 3;
  case CmpInst::FCMP_OLE: return 4;
  case CmpInst::FCMP_OGT: return 5;
  case CmpInst::FCMP_OLT: return 6;
  case CmpInst::FCMP_UGT: return 7;
  default: reject("unexpected floating-point predicate");
  }
}

bool allowedExternal(const Function &F) {
  StringRef N = F.getName();
  if (!F.isIntrinsic()) {
    if (N == "sirius_simd_fmaxf") return true;
    if (!N.consume_front(MathPrefix) || !N.consume_back("16")) return false;
    for (StringRef Recipe : {"dense_add", "dense_multiply", "dense_divide", "dense_sqrt",
                             "tagged_add", "tagged_multiply", "tagged_subtract", "tagged_divide", "tagged_sqrt"})
      if (N == Recipe) return true;
    return false;
  }
  switch (F.getIntrinsicID()) {
  case Intrinsic::fabs:
  case Intrinsic::memset:
  case Intrinsic::cttz:
  case Intrinsic::ctlz:
  case Intrinsic::masked_load:
  case Intrinsic::masked_store:
  case Intrinsic::lifetime_start:
  case Intrinsic::lifetime_end:
  case Intrinsic::experimental_noalias_scope_decl:
  case Intrinsic::usub_sat:
    return true;
  default:
    return N == "llvm.x86.avx512.mask.gather.qpi.512" || N == "llvm.x86.avx512.mask.scatter.qpi.512";
  }
}

bool extendedMask(const Value *Mask) {
  const auto *ZExt = dyn_cast<ZExtInst>(Mask);
  return ZExt && ZExt->getSrcTy()->isIntegerTy(Lanes) && ZExt->getDestTy()->isIntegerTy(32);
}

bool admitPackedCall(CallBase &Call, StringRef Name, SmallPtrSetImpl<const AllocaInst *> &Banks) {
  const bool Tagged = Name.drop_front(MathPrefix.size()).starts_with("tagged_");
  const unsigned First = Tagged ? 2 : 1;
  if (!Call.getType()->isVoidTy() || Call.arg_size() != First + 4)
    reject("unexpected raw pointer ABI", &Call);
  for (unsigned Arg = 0; Arg != First; ++Arg) {
    const auto *Bank = dyn_cast<AllocaInst>(Call.getArgOperand(Arg));
    const auto *Array = Bank ? dyn_cast<ArrayType>(Bank->getAllocatedType()) : nullptr;
    const auto *Count = Bank ? dyn_cast<ConstantInt>(Bank->getArraySize()) : nullptr;
    const unsigned Words = (Arg == 0 ? 15 : 3) * Lanes;
    if (!Array || !Count || !Count->isOne() || Bank->getAddressSpace() != 0 ||
        !Array->getElementType()->isIntegerTy(32) || Array->getNumElements() != Words ||
        !Banks.insert(Bank).second)
      reject("raw helper does not own exact, separate payload and tag banks", &Call);
  }
  const std::array<unsigned, 3> Slots = Tagged ? std::array<unsigned, 3>{0, 1, 2}
                                               : std::array<unsigned, 3>{2, 0, 1};
  for (unsigned Arg = 0; Arg != Slots.size(); ++Arg) {
    const auto *Slot = dyn_cast<ConstantInt>(Call.getArgOperand(First + Arg));
    if (!Slot || !Slot->getType()->isIntegerTy(32) || Slot->getZExtValue() != Slots[Arg])
      reject("raw math slot tuple changed", &Call);
  }
  const Value *Mask = Call.getArgOperand(First + 3);
  if (!extendedMask(Mask)) {
    // The Transport dense path merges two originally zero-extended masks.
    const auto *Phi = dyn_cast<PHINode>(Mask);
    if (!Phi || Phi->getNumIncomingValues() != 2 ||
        !extendedMask(Phi->getIncomingValue(0)) || !extendedMask(Phi->getIncomingValue(1)))
      reject("raw math mask can carry bits outside sixteen lanes", &Call);
  }
  return Tagged;
}

// Keep the original integer operand. Constants and the four table broadcasts
// must prove every element defined and exactly representable before conversion.
void admitCoefficient(SIToFPInst &I, const Stage &S, StringMap<unsigned> &Tables) {
  if (&S != &Transport || !wide(I.getSrcTy(), false) || !wide(I.getDestTy(), true))
    reject("unexpected integer conversion", &I);
  if (auto *Argument = dyn_cast<llvm::Argument>(I.getOperand(0))) {
    const Function *F = Argument->getParent();
    if (F != I.getFunction() || !F->hasLocalLinkage() || F->getName() != "Rational_0___vyivyi" ||
        F->arg_size() != 3 || Argument->getArgNo() > 1 || F->getNumUses() != 14)
      reject("unexpected Rational coefficient source", &I);
    for (const User *User : F->users()) {
      auto *Call = dyn_cast<CallInst>(User);
      if (!Call || Call->getCalledFunction() != F || Call->arg_size() != 3)
        reject("Rational has an indirect or incompatible caller", &I);
      for (unsigned Arg = 0; Arg != 2; ++Arg) {
        auto *K = dyn_cast<Constant>(Call->getArgOperand(Arg));
        auto *First = K && wide(K->getType(), false) ? dyn_cast_or_null<ConstantInt>(K->getAggregateElement(0U)) : nullptr;
        if (!First || First->getSExtValue() < -64448 || First->getSExtValue() > 64448)
          reject("Rational coefficient is outside the proved range", &I);
        for (unsigned Lane = 1; Lane != Lanes; ++Lane)
          if (K->getAggregateElement(Lane) != First) reject("Rational coefficient is not a defined splat", &I);
      }
    }
    ++Tables[Argument->getArgNo() ? "Rational.arg1" : "Rational.arg0"];
    return;
  }

  auto *Shuffle = dyn_cast<ShuffleVectorInst>(I.getOperand(0));
  auto *Insert = Shuffle ? dyn_cast<InsertElementInst>(Shuffle->getOperand(0)) : nullptr;
  auto *Index = Insert ? dyn_cast<ConstantInt>(Insert->getOperand(2)) : nullptr;
  auto *Load = Insert ? dyn_cast<LoadInst>(Insert->getOperand(1)) : nullptr;
  if (!I.getFunction()->getName().starts_with("ComputeOwnedRow_0___") ||
      !Shuffle || !isa<PoisonValue>(Shuffle->getOperand(1)) || !Insert ||
      !isa<PoisonValue>(Insert->getOperand(0)) || !Index || !Index->isZero() ||
      !Load || !Load->isSimple() || !Load->getType()->isIntegerTy(32) ||
      Shuffle->getShuffleMask().size() != Lanes)
    reject("coefficient is not a broadcast of an initialized scalar", &I);
  for (int Lane : Shuffle->getShuffleMask())
    if (Lane != 0) reject("coefficient broadcast reaches an undefined element", &I);
  auto *GEP = dyn_cast<GetElementPtrInst>(Load->getPointerOperand());
  auto *Table = GEP ? dyn_cast<GlobalVariable>(GEP->getPointerOperand()) : nullptr;
  auto *Type = Table ? dyn_cast<ArrayType>(Table->getValueType()) : nullptr;
  if (!Table || !Table->isConstant() || !Table->hasLocalLinkage() || !Table->hasInitializer() ||
      !Type || Type->getNumElements() != 7 || !Type->getElementType()->isIntegerTy(32) ||
      GEP->getNumIndices() != 1 ||
      (Table->getName() != "cn_0" && Table->getName() != "cd_0" &&
       Table->getName() != "en_0" && Table->getName() != "ed_0"))
    reject("unknown Transport coefficient table", &I);
  for (unsigned K = 0; K != 7; ++K) {
    auto *Value = dyn_cast_or_null<ConstantInt>(Table->getInitializer()->getAggregateElement(K));
    if (!Value || Value->getSExtValue() < -339200 || Value->getSExtValue() > 339200)
      reject("undefined or inexact Transport coefficient", &I);
  }
  // The sealed loop bounds the indices to 1..6 (cn/cd) or 2..6 (en/ed).
  ++Tables[Table->getName()];
}

void admitMultiply(Instruction &I, const Stage &S) {
  auto *K = dyn_cast<Constant>(I.getOperand(1));
  if (&S != &Endpoint || !K) reject("unexpected Endpoint scale operand", &I);
  for (unsigned Lane = 0; Lane != Lanes; ++Lane) {
    auto *Value = dyn_cast_or_null<ConstantFP>(K->getAggregateElement(Lane));
    if (!Value || Value->getValueAPF().bitcastToAPInt().getZExtValue() != 0x29800000U)
      reject("Endpoint scale is outside the original 2^-44 operation", &I);
  }
}

void constrain(Instruction &Old) {
  IRBuilder<> Builder(&Old);
  Builder.setIsFPConstrained(true);
  Builder.setDefaultConstrainedExcept(fp::ebStrict);
  Builder.setDefaultConstrainedRounding(RoundingMode::Dynamic);
  CallInst *New;
  if (isa<SIToFPInst>(Old)) {
    New = Builder.CreateConstrainedFPCast(Intrinsic::experimental_constrained_sitofp,
        Old.getOperand(0), Old.getType(), {}, "", nullptr, RoundingMode::Dynamic, fp::ebStrict);
  } else {
    // As in the eight-lane scope, only FP operands are frozen. The live-lane
    // context proof and untouched MXCSR restore belong to the downstream gates.
    Value *A = Builder.CreateFreeze(Old.getOperand(0), "wide.fp.freeze");
    Value *B = Builder.CreateFreeze(Old.getOperand(1), "wide.fp.freeze");
    if (auto *Cmp = dyn_cast<FCmpInst>(&Old)) {
      New = Builder.CreateConstrainedFPCmp(predicate(Cmp->getPredicate()) >= 3
          ? Intrinsic::experimental_constrained_fcmps : Intrinsic::experimental_constrained_fcmp,
          Cmp->getPredicate(), A, B, "", fp::ebStrict);
    } else {
      Intrinsic::ID Id = Old.getOpcode() == Instruction::FAdd ? Intrinsic::experimental_constrained_fadd
          : Old.getOpcode() == Instruction::FSub ? Intrinsic::experimental_constrained_fsub
          : Intrinsic::experimental_constrained_fmul;
      New = Builder.CreateConstrainedFPBinOp(Id, A, B, {}, "", nullptr, RoundingMode::Dynamic, fp::ebStrict);
    }
  }
  New->takeName(&Old);
  Old.replaceAllUsesWith(New);
  Old.eraseFromParent();
}

struct StrictWide : PassInfoMixin<StrictWide> {
  PreservedAnalyses run(Module &M, ModuleAnalysisManager &) {
    const Stage &S = admit(M);
    if (M.getTargetTriple().str() != "x86_64-unknown-linux-gnu" || M.getDataLayoutStr() !=
        "e-m:e-p270:32:32-p271:32:32-p272:64:64-i64:64-i128:128-f80:128-n8:16:32:64-S128" ||
        !M.alias_empty() || !M.ifunc_empty() || !M.getModuleInlineAsm().empty())
      reject("unexpected raw target, layout or linkage");

    SmallVector<Function *> Definitions, Private;
    SmallVector<Instruction *> Work;
    SmallVector<CallInst *> Fmax;
    SmallPtrSet<const AllocaInst *, 32> Banks;
    std::array<unsigned, 6> Operations{};
    std::array<unsigned, 8> Predicates{};
    StringMap<unsigned> Tables;
    unsigned Exports = 0, PointerCalls = 0, TaggedCalls = 0;
    std::string Entry = (Prefix + S.Entry + "_packet").str();
    for (Function &F : M) {
      if (F.hasFnAttribute(Attribute::StrictFP)) reject("input already has strictfp");
      if (F.isDeclaration()) {
        if (!allowedExternal(F)) reject(Twine("unexpected declaration ") + F.getName());
        if (F.getName().starts_with(MathPrefix)) Private.push_back(&F);
        continue;
      }
      Definitions.push_back(&F);
      if (!F.hasLocalLinkage()) {
        if (F.getName() != Entry && F.getName() != Entry + MaskedSuffix.str()) reject("unexpected export");
        ++Exports;
        Private.push_back(&F);
      }
      for (Instruction &I : instructions(F)) {
        if (auto *Op = dyn_cast<FPMathOperator>(&I))
          if (Op->getFastMathFlags().any() || I.getMetadata(LLVMContext::MD_fpmath))
            reject("fast or approximate floating-point operation", &I);
        const unsigned Op = I.getOpcode();
        if (Op == Instruction::FAdd || Op == Instruction::FSub || Op == Instruction::FMul || Op == Instruction::FCmp) {
          if (!wide(I.getOperand(0)->getType(), true) || I.getOperand(1)->getType() != I.getOperand(0)->getType())
            reject("unexpected SIMD16 FP operands", &I);
          ++Operations[Op == Instruction::FAdd ? 0 : Op == Instruction::FSub ? 1 : Op == Instruction::FMul ? 2 : 3];
          if (Op == Instruction::FMul) admitMultiply(I, S);
          if (auto *Cmp = dyn_cast<FCmpInst>(&I)) ++Predicates[predicate(Cmp->getPredicate())];
          Work.push_back(&I);
        } else if (Op == Instruction::FNeg) {
          if (!wide(I.getType(), true)) reject("unexpected fneg type", &I);
          ++Operations[4];
        } else if (auto *Cast = dyn_cast<SIToFPInst>(&I)) {
          admitCoefficient(*Cast, S, Tables);
          ++Operations[5];
          Work.push_back(&I);
        } else if ((I.isBinaryOp() && I.getType()->isFPOrFPVectorTy()) ||
                   (isa<CastInst>(I) && Op != Instruction::BitCast &&
                    (I.getType()->isFPOrFPVectorTy() || I.getOperand(0)->getType()->isFPOrFPVectorTy()))) {
          reject("unsupported floating-point arithmetic or cast", &I);
        }
        if (auto *Call = dyn_cast<CallBase>(&I)) {
          auto *Target = Call->getCalledFunction();
          if (!isa<CallInst>(Call) || !Target || Call->hasOperandBundles())
            reject("unsupported indirect, assembly, invoke or bundled call", &I);
          if (Target->getName().starts_with(MathPrefix)) {
            ++PointerCalls;
            TaggedCalls += admitPackedCall(*Call, Target->getName(), Banks);
          }
          if (Target->getIntrinsicID() == Intrinsic::fabs && !wide(Call->getType(), true))
            reject("unexpected fabs vector type", &I);
          if (Target->getName() == "sirius_simd_fmaxf") {
            if (Call->arg_size() != 2 || !Call->getType()->isFloatTy() ||
                !Call->getArgOperand(0)->getType()->isFloatTy() || !Call->getArgOperand(1)->getType()->isFloatTy())
              reject("opaque fmaxf ABI changed", &I);
            Fmax.push_back(cast<CallInst>(Call));
          }
        }
      }
    }

    if (Operations != S.Operations || Predicates != S.Predicates || Definitions.size() != S.Definitions ||
        Exports != 2 || PointerCalls != S.PointerCalls || TaggedCalls != S.TaggedCalls)
      reject("SIMD16 operation or call inventory changed");
    if (&S == &Transport) {
      for (StringRef Name : {"Rational.arg0", "Rational.arg1", "cn_0", "cd_0", "en_0", "ed_0"})
        if (Tables.lookup(Name) != 1) reject("Transport coefficient site missing");
    }

    for (Function *F : Definitions) {
      F->addFnAttr(Attribute::StrictFP);
      F->setMemoryEffects(F->getMemoryEffects() | MemoryEffects::inaccessibleMemOnly());
    }
    for (Instruction *I : Work) constrain(*I);
    for (CallInst *Call : Fmax) {
      IRBuilder<> Builder(Call);
      for (unsigned Arg = 0; Arg != 2; ++Arg)
        Call->setArgOperand(Arg, Builder.CreateFreeze(Call->getArgOperand(Arg), "wide.fp.freeze"));
    }
    for (Function *F : Definitions) {
      for (Instruction &I : instructions(F)) {
        if (auto *Call = dyn_cast<CallBase>(&I)) {
          Call->addFnAttr(Attribute::StrictFP);
          if (!Call->getCalledFunction()->isIntrinsic())
            Call->setMemoryEffects(Call->getMemoryEffects() | MemoryEffects::inaccessibleMemOnly());
        }
      }
    }
    for (Function *F : Private) {
      std::string Name = ("sirius_wide_strict_" + F->getName().drop_front(Prefix.size())).str();
      if (M.getNamedValue(Name)) reject("private symbol collision");
      F->setName(Name);
      F->setVisibility(GlobalValue::HiddenVisibility);
      F->setDSOLocal(true);
    }
    if (verifyModule(M, &errs())) reject("LLVM rejected the constrained module");
    errs() << "SIMD16 strict scope: " << S.Filename << " converted=" << Work.size()
           << " pointer_calls=" << PointerCalls << " tagged_calls=" << TaggedCalls
           << " fmaxf=" << Fmax.size() << '\n';
    return PreservedAnalyses::none();
  }
};
} // namespace

extern "C" LLVM_ATTRIBUTE_WEAK PassPluginLibraryInfo llvmGetPassPluginInfo() {
  return {LLVM_PLUGIN_API_VERSION, "SiriusStrictWide", LLVM_VERSION_STRING,
          [](PassBuilder &Builder) {
            Builder.registerPipelineParsingCallback(
                [](StringRef Name, ModulePassManager &Passes, ArrayRef<PassBuilder::PipelineElement> Inner) {
                  if (Name != "sirius-wide-strict" || !Inner.empty()) return false;
                  Passes.addPass(StrictWide());
                  return true;
                });
          }, nullptr};
}
