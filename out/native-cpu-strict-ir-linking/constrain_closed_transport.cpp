// Fixed LLVM23.1.2 private Transport/Endpoint copies, source preparation only.
// Source-semantic proof and exact inventories: remaining_scope_review.json.
// No public object, pointer, integer input, control mask or memory op changes.
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
#include "llvm/IR/PassManager.h"
#include "llvm/IR/Verifier.h"
#include "llvm/Passes/PassBuilder.h"
#include "llvm/Plugins/PassPlugin.h"
#include "llvm/Support/ErrorHandling.h"
#include "llvm/Support/MemoryBuffer.h"
#include "llvm/Support/ModRef.h"
#include "llvm/Support/SHA256.h"
#include "llvm/Support/SourceMgr.h"
#include "llvm/Support/raw_ostream.h"
#include <array>

#if LLVM_VERSION_MAJOR != 23 || LLVM_VERSION_MINOR != 1 || LLVM_VERSION_PATCH != 2
#error "The reviewed APIs, parser and lowering require LLVM 23.1.2."
#endif

using namespace llvm;

namespace {
struct NamedCount { StringLiteral Name; unsigned Count; };
const NamedCount TransportCalls[] = {
    {"sirius_native_v4_tagged_add8", 947},
    {"sirius_native_v4_tagged_multiply8", 1482},
    {"sirius_native_v4_tagged_subtract8", 10},
    {"sirius_native_v4_tagged_divide8", 240},
    {"sirius_native_v4_tagged_sqrt8", 3},
    {"sirius_native_v4_dense_divide8", 3},
    {"sirius_native_v4_dense_add8", 7},
    {"sirius_native_v4_dense_multiply8", 6},
    {"sirius_simd_fmaxf", 10},
    {"llvm.memset.p0.i32", 18},
    {"llvm.cttz.i64", 11},
    {"llvm.ctlz.v8i32", 2},
    {"llvm.masked.load.v8i32.p0", 29579},
    {"llvm.masked.store.v8i32.p0", 36979},
    {"llvm.x86.avx512.mask.gather3div8.si", 10},
    {"llvm.x86.avx512.mask.scatterdiv8.si", 16},
    {"llvm.lifetime.start.p0", 5380},
    {"llvm.lifetime.end.p0", 5380},
    {"llvm.experimental.noalias.scope.decl", 26514},
    {"llvm.fabs.v8f32", 308},
    {"llvm.usub.sat.v8i32", 892}};
const NamedCount EndpointCalls[] = {
    {"sirius_native_v4_tagged_add8", 1607},
    {"sirius_native_v4_tagged_multiply8", 1963},
    {"sirius_native_v4_tagged_subtract8", 197},
    {"sirius_native_v4_tagged_divide8", 154},
    {"sirius_native_v4_tagged_sqrt8", 3},
    {"sirius_native_v4_dense_divide8", 7},
    {"sirius_native_v4_dense_sqrt8", 2},
    {"sirius_native_v4_dense_add8", 13},
    {"sirius_native_v4_dense_multiply8", 13},
    {"sirius_simd_fmaxf", 12},
    {"llvm.memset.p0.i32", 37},
    {"llvm.cttz.i64", 13},
    {"llvm.ctlz.v8i32", 5},
    {"llvm.masked.load.v8i32.p0", 45730},
    {"llvm.masked.store.v8i32.p0", 52028},
    {"llvm.masked.store.v8f32.p0", 8},
    {"llvm.x86.avx512.mask.gather3div8.si", 230},
    {"llvm.x86.avx512.mask.scatterdiv8.si", 48},
    {"llvm.lifetime.start.p0", 7894},
    {"llvm.lifetime.end.p0", 7894},
    {"llvm.experimental.noalias.scope.decl", 39154},
    {"llvm.fabs.v8f32", 311},
    {"llvm.usub.sat.v8i32", 404},
    {"llvm.memset.p0.i64", 1}};

struct Config {
  StringLiteral Pass, Input, SHA, Source, OldExport, NewExport;
  uint64_t Bytes;
  unsigned Definitions, TotalCalls;
  // add, sub, mul, compare, negate, signed conversion, fabs
  std::array<unsigned, 7> Ops;
  // oeq, une, uno, oge, ole, ogt, olt
  std::array<unsigned, 7> Predicates;
  ArrayRef<NamedCount> Calls;
};
const Config Transport{
    "sirius-constrain-closed-transport",
    "/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-canonical-closed-transport/transport.ll",
    "2f409751b514fa3fd8fe9adba623a6cfd9ae3e3adad056251c9ba09b97529562",
    "/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-canonical-closed-transport/ispc-source/transport.ispc",
    "sirius_simd_closed_transport_packet", "sirius_closed_transport_packet",
    82880982, 156, 107992, {79, 13, 0, 362, 7, 6, 308}, {304, 12, 10, 22, 2, 12, 0}, TransportCalls};
const Config Endpoint{
    "sirius-constrain-endpoint",
    "/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-strict-ir-linking/endpoint.ll",
    "f1af24efea9ee6632063486343addb746610f53cee9d66ce4761dad2520bd014",
    "/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-inline-private-bank/ispc-source/endpoint.ispc",
    "sirius_simd_v4_endpoint_packet", "sirius_scope_endpoint_packet",
    117127171, 204, 157944, {99, 13, 1, 484, 16, 0, 311},
    {431, 10, 12, 9, 3, 16, 3}, EndpointCalls};
constexpr StringLiteral OldMath = "sirius_native_v4_";
constexpr StringLiteral NewMath = "sirius_closed_transport_scope_v4_";
constexpr StringLiteral MaskedSuffix =
    "___un_3C_unu_3E_unuun_3C_unu_3E_unuunuunu";

[[noreturn]] void reject(const Twine &Why, const Instruction *I = nullptr) {
  errs() << "sirius-constrain-remaining: " << Why << '\n';
  if (I) {
    errs() << "in " << I->getFunction()->getName() << ": ";
    I->print(errs());
    errs() << '\n';
  }
  report_fatal_error("closed remaining-stage conversion refused", false);
}

std::string digest(StringRef Text) {
  SHA256 Hash;
  Hash.update(Text);
  return toHex(Hash.final(), true);
}
std::string fingerprint(const Module &M) {
  std::string Text;
  raw_string_ostream OS(Text);
  M.print(OS, nullptr);
  OS.flush();
  return digest(Text);
}
void requirePinnedModule(const Module &M, const Config &C) {
  // Producer also checks the raw SHA before any rewrite.  Full parsed-module
  // equality binds every initializer, use, mask and CFG edge, not just counts.
  auto Buffer = MemoryBuffer::getFile(C.Input);
  if (!Buffer || (*Buffer)->getBufferSize() != C.Bytes ||
      digest((*Buffer)->getBuffer()) != C.SHA)
    reject("original module source pin absent or changed");
  LLVMContext ReferenceContext;
  SMDiagnostic Diagnostic;
  auto Reference = parseAssembly((*Buffer)->getMemBufferRef(), Diagnostic,
                                 ReferenceContext);
  if (!Reference)
    reject("pinned original cannot be parsed with LLVM23.1.2");
  Reference->setModuleIdentifier(M.getModuleIdentifier());
  std::string Expected = fingerprint(*Reference);
  if (fingerprint(M) != Expected)
    reject("full module fingerprint differs from pinned original");
  errs() << C.Pass << ": source-sha256=" << C.SHA
         << " module-fingerprint=" << Expected << '\n';
}
bool v8(Type *T, bool FP) {
  auto *V = dyn_cast<FixedVectorType>(T);
  return V && V->getNumElements() == 8 &&
         (FP ? V->getElementType()->isFloatTy()
             : V->getElementType()->isIntegerTy(32));
}
unsigned predicateIndex(CmpInst::Predicate P) {
  switch (P) {
  case CmpInst::FCMP_OEQ: return 0;
  case CmpInst::FCMP_UNE: return 1;
  case CmpInst::FCMP_UNO: return 2;
  case CmpInst::FCMP_OGE: return 3;
  case CmpInst::FCMP_OLE: return 4;
  case CmpInst::FCMP_OGT: return 5;
  case CmpInst::FCMP_OLT: return 6;
  default: reject("predicate outside the reviewed source vocabulary");
  }
}
const NamedCount *external(StringRef Name, const Config &C) {
  for (const NamedCount &Row : C.Calls)
    if (Row.Name == Name)
      return &Row;
  return nullptr;
}
std::string renamed(StringRef Name, const Config &C) {
  return Name.starts_with(OldMath)
             ? (NewMath + Name.drop_front(OldMath.size())).str()
             : (C.NewExport + Name.drop_front(C.OldExport.size())).str();
}

// Six exact closed-Transport casts. Never freeze an integer operand.
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
void requireEndpointMultiply(Instruction &I, const Config &C) {
  auto *K = dyn_cast<Constant>(I.getOperand(1));
  if (&C != &Endpoint || !K)
    reject("unexpected multiplication site", &I);
  for (unsigned Lane = 0; Lane != 8; ++Lane) {
    auto *FP = dyn_cast_or_null<ConstantFP>(K->getAggregateElement(Lane));
    if (!FP || FP->getValueAPF().bitcastToAPInt().getZExtValue() != 0x29800000U)
      reject("Endpoint multiply is not the original FP32 scale by2^-44", &I);
  }
}

struct CastSite { SIToFPInst *I; LoadInst *Scalar; };
struct ConstrainRemaining : PassInfoMixin<ConstrainRemaining> {
  const Config &C;
  explicit ConstrainRemaining(const Config &C) : C(C) {}
  PreservedAnalyses run(Module &M, ModuleAnalysisManager &) {
    requirePinnedModule(M, C);
    if (M.getSourceFileName() != C.Source ||
        M.getTargetTriple().str() != "x86_64-unknown-linux-gnu" ||
        M.getDataLayoutStr() !=
            "e-m:e-p270:32:32-p271:32:32-p272:64:64-i64:64-i128:128-"
            "f80:128-n8:16:32:64-S128" ||
        !M.alias_empty() || !M.ifunc_empty() || M.hasModuleInlineAsm())
      reject("source/target/layout is outside the fixed scope");

    SmallVector<Instruction *, 600> Work;
    SmallVector<CastSite, 4> Casts;
    SmallVector<CallInst *, 12> FmaxSites;
    SmallVector<Function *, 204> Definitions;
    SmallVector<Function *, 11> Renames;
    std::array<unsigned, 7> Ops{}, Predicates{};
    StringMap<unsigned> ExternalCalls, Tables;
    unsigned TotalCalls = 0, DeclarationCount = 0, ExportCount = 0;

    // All checks precede every IR mutation.  Only exact reviewed intrinsics
    // and opaque numerical targets are admitted; fenv/inline asm is rejected.
    for (Function &F : M) {
      if (F.hasFnAttribute(Attribute::StrictFP))
        reject("input already strict; conversion must run first");
      if (F.isDeclaration()) {
        ++DeclarationCount;
        if (!external(F.getName(), C))
          reject(Twine("unexpected external/observer ") + F.getName());
        if (F.getName().starts_with(OldMath))
          Renames.push_back(&F);
        continue;
      }
      Definitions.push_back(&F);
      if (!F.hasLocalLinkage()) {
        if (F.getName() != C.OldExport &&
            F.getName() != (C.OldExport + MaskedSuffix).str())
          reject("unexpected defined export");
        ++ExportCount;
        Renames.push_back(&F);
      }
      for (Instruction &I : instructions(F)) {
        if (auto *FM = dyn_cast<FPMathOperator>(&I))
          if (FM->getFastMathFlags().any() ||
              I.getMetadata(LLVMContext::MD_fpmath))
            reject("fast/approximate FP is unsupported", &I);
        unsigned Op = I.getOpcode();
        if (Op == Instruction::FAdd || Op == Instruction::FSub ||
            Op == Instruction::FMul || Op == Instruction::FCmp) {
          if (!v8(I.getOperand(0)->getType(), true) ||
              I.getOperand(1)->getType() != I.getOperand(0)->getType())
            reject("FP operand type/width changed", &I);
          unsigned Slot = Op == Instruction::FAdd ? 0 :
                          Op == Instruction::FSub ? 1 :
                          Op == Instruction::FMul ? 2 : 3;
          ++Ops[Slot];
          if (Op == Instruction::FMul)
            requireEndpointMultiply(I, C);
          if (auto *Cmp = dyn_cast<FCmpInst>(&I))
            ++Predicates[predicateIndex(Cmp->getPredicate())];
          Work.push_back(&I);
        } else if (Op == Instruction::FNeg) {
          if (!v8(I.getType(), true)) reject("unexpected fneg type", &I);
          ++Ops[4]; // Exception-free representation operation remains exact.
        } else if (auto *SI = dyn_cast<SIToFPInst>(&I)) {
          ++Ops[5];
          Casts.push_back({SI, coefficientScalar(*SI, C, Tables)});
        } else if (I.isBinaryOp() && I.getType()->isFPOrFPVectorTy()) {
          reject("unknown FP arithmetic opcode", &I);
        } else if (isa<CastInst>(I) && Op != Instruction::BitCast &&
                   (I.getType()->isFPOrFPVectorTy() ||
                    I.getOperand(0)->getType()->isFPOrFPVectorTy())) {
          reject("unknown FP conversion", &I);
        }
        if (auto *CB = dyn_cast<CallBase>(&I)) {
          Function *Target = CB->getCalledFunction();
          if (!isa<CallInst>(CB) || !Target || CB->hasOperandBundles())
            reject("indirect/assembly/invoke/bundled call", &I);
          ++TotalCalls;
          if (Target->isDeclaration()) {
            if (!external(Target->getName(), C))
              reject("unexpected external call or FP observer", &I);
            ++ExternalCalls[Target->getName()];
          }
          if (Target->getName() == "sirius_simd_fmaxf") {
            if (CB->arg_size() != 2 || !CB->getType()->isFloatTy())
              reject("opaque fmaxf ABI changed", &I);
            for (Value *A : CB->args())
              if (!A->getType()->isFloatTy())
                reject("opaque fmaxf argument type changed", &I);
            FmaxSites.push_back(cast<CallInst>(CB));
          }
          if (Target->getName() == "llvm.fabs.v8f32") {
            if (!v8(CB->getType(), true)) reject("fabs type changed", &I);
            ++Ops[6];
          }
        }
      }
    }
    if (Ops != C.Ops || Predicates != C.Predicates ||
        Definitions.size() != C.Definitions || TotalCalls != C.TotalCalls ||
        DeclarationCount != C.Calls.size() || ExportCount != 2 ||
        Renames.size() != 10 || ExternalCalls.size() != C.Calls.size())
      reject("exact stage opcode/call/definition inventory changed");
    for (const NamedCount &Row : C.Calls)
      if (ExternalCalls.lookup(Row.Name) != Row.Count)
        reject(Twine("external call count changed: ") + Row.Name);
    if (&C == &Transport &&
        (Tables.size() != 6 || Tables.lookup("cn_0") != 1 ||
         Tables.lookup("cd_0") != 1 || Tables.lookup("en_0") != 1 ||
         Tables.lookup("ed_0") != 1 || Tables.lookup("Rational.arg0") != 1 ||
         Tables.lookup("Rational.arg1") != 1))
      reject("four coefficient-table cast joins changed");
    for (Function *F : Renames)
      if (M.getNamedValue(renamed(F->getName(), C)))
        reject("private-name collision");

    for (Function *F : Definitions) {
      F->addFnAttr(Attribute::StrictFP);
      F->setMemoryEffects(F->getMemoryEffects() |
                          MemoryEffects::inaccessibleMemOnly());
    }
    for (Instruction *Old : Work) {
      IRBuilder<> B(Old);
      B.setIsFPConstrained(true);
      B.setDefaultConstrainedExcept(fp::ebStrict);
      B.setDefaultConstrainedRounding(RoundingMode::Dynamic);
      Value *A = B.CreateFreeze(Old->getOperand(0), "scope.fp.freeze");
      Value *D = B.CreateFreeze(Old->getOperand(1), "scope.fp.freeze");
      CallInst *New;
      if (auto *Cmp = dyn_cast<FCmpInst>(Old)) {
        New = B.CreateConstrainedFPCmp(
            predicateIndex(Cmp->getPredicate()) >= 3
                ? Intrinsic::experimental_constrained_fcmps
                : Intrinsic::experimental_constrained_fcmp,
            Cmp->getPredicate(), A, D, "", fp::ebStrict);
      } else {
        Intrinsic::ID ID = Old->getOpcode() == Instruction::FAdd
                               ? Intrinsic::experimental_constrained_fadd
                           : Old->getOpcode() == Instruction::FSub
                               ? Intrinsic::experimental_constrained_fsub
                               : Intrinsic::experimental_constrained_fmul;
        New = B.CreateConstrainedFPBinOp(ID, A, D, {}, "", nullptr,
                                         RoundingMode::Dynamic, fp::ebStrict);
      }
      New->takeName(Old);
      Old->replaceAllUsesWith(New);
      Old->eraseFromParent();
    }
    for (CastSite Site : Casts) {
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
    for (CallInst *Call : FmaxSites) {
      IRBuilder<> B(Call);
      for (unsigned K = 0; K != 2; ++K)
        Call->setArgOperand(K, B.CreateFreeze(Call->getArgOperand(K),
                                             "scope.fp.freeze"));
    }
    for (Function *F : Definitions)
      for (Instruction &I : instructions(F))
        if (auto *CB = dyn_cast<CallBase>(&I)) {
          CB->addFnAttr(Attribute::StrictFP);
          if (!CB->getCalledFunction()->isIntrinsic())
            CB->setMemoryEffects(CB->getMemoryEffects() |
                                 MemoryEffects::inaccessibleMemOnly());
        }
    for (Function *F : Renames) {
      F->setName(renamed(F->getName(), C));
      F->setVisibility(GlobalValue::HiddenVisibility);
      F->setDSOLocal(true);
    }
    if (verifyModule(M, &errs())) reject("LLVM rejected converted private copy");
    return PreservedAnalyses::none();
  }
};
} // namespace

extern "C" LLVM_ATTRIBUTE_WEAK PassPluginLibraryInfo llvmGetPassPluginInfo() {
  return {LLVM_PLUGIN_API_VERSION, "SiriusScopedRemaining", LLVM_VERSION_STRING,
          [](PassBuilder &PB) {
            PB.registerPipelineParsingCallback(
                [](StringRef Name, ModulePassManager &MPM,
                   ArrayRef<PassBuilder::PipelineElement> Inner) {
                  if (!Inner.empty()) return false;
                  if (Name != Transport.Pass) return false;
                  MPM.addPass(ConstrainRemaining(Transport));
                  return true;
                });
          }, nullptr};
}
