// Source-only LLVM 23.1.2 new-PM pass for the separately named DP copy.
// This is NOT an ordinary/raw entrypoint replacement.  See scope_contract.json.
// The exact-module source proof (not a machine-checked lane analysis) is:
// 120832 context bytes are zeroed, admission guarantees each active row's
// 4149-word clear, and original masks/initialization define every live value.
// FP-only freeze refines excluded undef/poison choices, leaving defined bits
// unchanged.  Extra private flags are contained by the full-MXCSR scope.
// There is no force switch and no permission for a different input module.
#include "llvm/ADT/StringExtras.h"
#include "llvm/ADT/SmallVector.h"
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
#error "This source requires the reviewed LLVM 23.1.2 API and lowering."
#endif

using namespace llvm;

namespace {
constexpr StringLiteral OldExport = "sirius_simd_v4_dopri_phase_packet";
constexpr StringLiteral NewExport = "sirius_scope_dp_packet";
constexpr StringLiteral OldMathPrefix = "sirius_native_v4_";
constexpr StringLiteral NewMathPrefix = "sirius_scope_v4_";
constexpr StringLiteral Source =
    "/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-inline-private-bank/"
    "ispc-source/dopri_phase.ispc";
constexpr StringLiteral PinnedInput =
    "/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-strict-ir-linking/"
    "dopri_phase.ll";
constexpr StringLiteral PinnedSHA =
    "be998fd307e38ca7944aee9150db28ab88c627690873ee0343be38914df79401";

[[noreturn]] void reject(const Twine &Why, const Instruction *I = nullptr) {
  errs() << "sirius-constrain-dp: " << Why << '\n';
  if (I) {
    errs() << "in " << I->getFunction()->getName() << ": ";
    I->print(errs());
    errs() << '\n';
  }
  report_fatal_error("closed DP scope conversion refused", false);
}

bool v8f32(Type *T) {
  auto *V = dyn_cast<FixedVectorType>(T);
  return V && V->getNumElements() == 8 && V->getElementType()->isFloatTy();
}

std::string digest(StringRef Text) {
  SHA256 Hash;
  Hash.update(Text);
  return toHex(Hash.final(), true);
}

std::string moduleFingerprint(const Module &M) {
  std::string Text;
  raw_string_ostream OS(Text);
  M.print(OS, nullptr);
  OS.flush();
  return digest(Text);
}

void requirePinnedModule(const Module &M) {
  // The producer also checks this SHA BEFORE any rewrite.  Comparing complete
  // canonical LLVM prints prevents an earlier pass or same-inventory edit
  // from inheriting the source proof.  A separate context avoids type renames.
  auto Buffer = MemoryBuffer::getFile(PinnedInput);
  if (!Buffer || (*Buffer)->getBufferSize() != 53181007 ||
      digest((*Buffer)->getBuffer()) != PinnedSHA)
    reject("original DP source pin is absent or changed");
  LLVMContext ReferenceContext;
  SMDiagnostic Diagnostic;
  auto Reference = parseAssembly((*Buffer)->getMemBufferRef(), Diagnostic,
                                 ReferenceContext);
  if (!Reference)
    reject("pinned original DP cannot be parsed by LLVM 23.1.2");
  Reference->setModuleIdentifier(M.getModuleIdentifier());
  std::string Expected = moduleFingerprint(*Reference);
  if (moduleFingerprint(M) != Expected)
    reject("full module fingerprint differs from the pinned original");
  errs() << "sirius-constrain-dp: source-sha256=" << PinnedSHA
         << " module-fingerprint=" << Expected << '\n';
}

// This is the exact quiet/signaling classification used by LLVM 23.1.2
// translateX86FSETCC, not a claim that ordinary FCmp specifies exceptions.
bool signaling(CmpInst::Predicate P) {
  switch (P) {
  case CmpInst::FCMP_OEQ:
  case CmpInst::FCMP_UNE:
  case CmpInst::FCMP_UNO:
    return false;
  case CmpInst::FCMP_OGE:
  case CmpInst::FCMP_OGT:
  case CmpInst::FCMP_OLE:
  case CmpInst::FCMP_UGT:
    return true;
  default:
    reject("comparison predicate outside the sealed DP vocabulary");
  }
}

bool mathName(StringRef N) {
  if (!N.consume_front(OldMathPrefix))
    return false;
  return N == "tagged_add8" || N == "tagged_subtract8" ||
         N == "tagged_multiply8" || N == "tagged_divide8" ||
         N == "dense_add8" || N == "dense_multiply8" ||
         N == "dense_divide8" || N == "dense_sqrt8";
}

bool allowedIntrinsic(Intrinsic::ID ID) {
  switch (ID) {
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
    return false;
  }
}

bool allowedX86MemoryIntrinsic(StringRef N) {
  return N == "llvm.x86.avx512.mask.gather3div8.si" ||
         N == "llvm.x86.avx512.mask.scatterdiv8.si";
}

struct ConstrainDPScope : PassInfoMixin<ConstrainDPScope> {
  PreservedAnalyses run(Module &M, ModuleAnalysisManager &) {
    requirePinnedModule(M);
    if (M.getSourceFileName() != Source ||
        M.getTargetTriple().str() != "x86_64-unknown-linux-gnu" ||
        M.getDataLayoutStr() !=
            "e-m:e-p270:32:32-p271:32:32-p272:64:64-i64:64-i128:128-"
            "f80:128-n8:16:32:64-S128" ||
        !M.alias_empty() || !M.ifunc_empty() || M.hasModuleInlineAsm())
      reject("module is outside the pinned standalone DP copy");

    SmallVector<Instruction *, 721> Work;
    SmallVector<CallInst *, 1> FmaxSites;
    SmallVector<Function *, 56> Definitions;
    SmallVector<Function *, 10> Renames;
    std::array<unsigned, 5> Counts{}; // add, sub, cmp, neg, fabs
    std::array<unsigned, 7> Predicates{}; // oeq, une, uno, oge, ole, ogt, ugt
    unsigned MathCalls = 0, FmaxCalls = 0, ExportCount = 0;

    // Validate EVERYTHING before modifying any IR.  No partial output is
    // accepted if an operation, call, type or exact source identity differs.
    for (Function &F : M) {
      if (F.hasFnAttribute(Attribute::StrictFP))
        reject("input already has strictfp; conversion must run first");
      if (F.isDeclaration()) {
        if (F.isIntrinsic()) {
          if (!allowedIntrinsic(F.getIntrinsicID()) &&
              !allowedX86MemoryIntrinsic(F.getName()))
            reject(Twine("unknown intrinsic ") + F.getName());
        } else if (!mathName(F.getName()) &&
                   F.getName() != "sirius_simd_fmaxf") {
          reject(Twine("unknown external function ") + F.getName());
        }
        if (mathName(F.getName()))
          Renames.push_back(&F);
        continue;
      }
      Definitions.push_back(&F);
      if (!F.hasLocalLinkage()) {
        if (F.getName() != OldExport &&
            !F.getName().starts_with((OldExport + "___").str()))
          reject(Twine("unexpected defined export ") + F.getName());
        ++ExportCount;
        Renames.push_back(&F);
      }
      for (Instruction &I : instructions(F)) {
        if (auto *FM = dyn_cast<FPMathOperator>(&I))
          if (FM->getFastMathFlags().any() ||
              I.getMetadata(LLVMContext::MD_fpmath))
            reject("fast/approximate FP metadata is outside this conversion", &I);

        unsigned Op = I.getOpcode();
        if (Op == Instruction::FAdd || Op == Instruction::FSub ||
            Op == Instruction::FCmp) {
          if (!v8f32(I.getOperand(0)->getType()) ||
              I.getOperand(1)->getType() != I.getOperand(0)->getType())
            reject("FP operands are not unchanged eight-lane float32", &I);
          if (auto *Cmp = dyn_cast<FCmpInst>(&I)) {
            (void)signaling(Cmp->getPredicate());
            switch (Cmp->getPredicate()) {
            case CmpInst::FCMP_OEQ: ++Predicates[0]; break;
            case CmpInst::FCMP_UNE: ++Predicates[1]; break;
            case CmpInst::FCMP_UNO: ++Predicates[2]; break;
            case CmpInst::FCMP_OGE: ++Predicates[3]; break;
            case CmpInst::FCMP_OLE: ++Predicates[4]; break;
            case CmpInst::FCMP_OGT: ++Predicates[5]; break;
            case CmpInst::FCMP_UGT: ++Predicates[6]; break;
            default: reject("unexpected comparison predicate", &I);
            }
          }
          ++Counts[Op == Instruction::FAdd ? 0 :
                   Op == Instruction::FSub ? 1 : 2];
          Work.push_back(&I);
        } else if (Op == Instruction::FNeg) {
          if (!v8f32(I.getType()))
            reject("unexpected fneg type", &I);
          ++Counts[3]; // Representation operation: retained, no freeze.
        } else if (I.isBinaryOp() && I.getType()->isFPOrFPVectorTy()) {
          reject("unsupported floating-point arithmetic opcode", &I);
        } else if (isa<CastInst>(I) && Op != Instruction::BitCast &&
                   (I.getType()->isFPOrFPVectorTy() ||
                    I.getOperand(0)->getType()->isFPOrFPVectorTy())) {
          reject("unexpected floating-point conversion", &I);
        }
        if (auto *CB = dyn_cast<CallBase>(&I)) {
          Function *Callee = CB->getCalledFunction();
          if (!isa<CallInst>(CB) || !Callee || CB->hasOperandBundles())
            reject("indirect/assembly/invoke/bundled call is unsupported", &I);
          if (mathName(Callee->getName()))
            ++MathCalls;
          if (Callee->getName() == "sirius_simd_fmaxf") {
            ++FmaxCalls;
            if (CB->arg_size() != 2 || !CB->getType()->isFloatTy())
              reject("opaque fmaxf ABI changed", &I);
            for (Value *A : CB->args())
              if (!A->getType()->isFloatTy())
                reject("opaque fmaxf argument type changed", &I);
            FmaxSites.push_back(cast<CallInst>(CB));
          }
          if (Callee->getIntrinsicID() == Intrinsic::fabs) {
            if (!v8f32(CB->getType()))
              reject("unexpected fabs type", &I);
            ++Counts[4]; // Bit representation operation, retained exactly.
          }
        }
      }
    }
    if (Counts != std::array<unsigned, 5>{15, 16, 690, 2, 1276} ||
        Predicates != std::array<unsigned, 7>{675, 1, 1, 5, 5, 2, 1} ||
        Definitions.size() != 56 || ExportCount != 2 || Renames.size() != 10 ||
        MathCalls != 1581 || FmaxCalls != 1)
      reject("sealed DP inventory changed");
    for (Function *F : Renames) {
      StringRef N = F->getName();
      std::string Name = N.starts_with(OldExport)
                             ? (NewExport + N.drop_front(OldExport.size())).str()
                             : (NewMathPrefix + N.drop_front(OldMathPrefix.size())).str();
      if (M.getNamedValue(Name))
        reject("private-name collision");
    }
    for (Function *F : Definitions) {
      F->addFnAttr(Attribute::StrictFP);
      // Constrained operations read/write the inaccessible FP environment.
      // Old ordinary-IR memory(read) annotations do not describe that effect.
      // Widen only this effect; preserve all original accessible-memory facts.
      F->setMemoryEffects(F->getMemoryEffects() |
                          MemoryEffects::inaccessibleMemOnly());
    }
    for (Instruction *I : Work) {
      Instruction &Old = *I;
      IRBuilder<> B(&Old);
      B.setIsFPConstrained(true);
      B.setDefaultConstrainedExcept(fp::ebStrict);
      B.setDefaultConstrainedRounding(RoundingMode::Dynamic);
      Value *A[2] = {Old.getOperand(0), Old.getOperand(1)};
      for (unsigned K = 0; K != 2; ++K)
        A[K] = B.CreateFreeze(A[K], "scope.fp.freeze");
      CallInst *New;
      if (auto *Cmp = dyn_cast<FCmpInst>(&Old)) {
        New = B.CreateConstrainedFPCmp(
            signaling(Cmp->getPredicate())
                ? Intrinsic::experimental_constrained_fcmps
                : Intrinsic::experimental_constrained_fcmp,
            Cmp->getPredicate(), A[0], A[1], "", fp::ebStrict);
      } else {
        New = B.CreateConstrainedFPBinOp(
            Old.getOpcode() == Instruction::FAdd
                ? Intrinsic::experimental_constrained_fadd
                : Intrinsic::experimental_constrained_fsub,
            A[0], A[1], {}, "", nullptr, RoundingMode::Dynamic, fp::ebStrict);
      }
      New->takeName(&Old);
      Old.replaceAllUsesWith(New);
      Old.eraseFromParent();
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
      StringRef N = F->getName();
      std::string Name = N.starts_with(OldExport)
                             ? (NewExport + N.drop_front(OldExport.size())).str()
                             : (NewMathPrefix + N.drop_front(OldMathPrefix.size())).str();
      F->setName(Name);
      F->setVisibility(GlobalValue::HiddenVisibility);
      F->setDSOLocal(true);
    }
    if (verifyModule(M, &errs()))
      reject("LLVM verification rejected the converted copy");
    return PreservedAnalyses::none();
  }
};
} // namespace

extern "C" LLVM_ATTRIBUTE_WEAK PassPluginLibraryInfo llvmGetPassPluginInfo() {
  return {LLVM_PLUGIN_API_VERSION, "SiriusScopedDP", LLVM_VERSION_STRING,
          [](PassBuilder &PB) {
            PB.registerPipelineParsingCallback(
                [](StringRef Name, ModulePassManager &MPM,
                   ArrayRef<PassBuilder::PipelineElement>) {
                  if (Name != "sirius-constrain-dp")
                    return false;
                  MPM.addPass(ConstrainDPScope());
                  return true;
                });
          }, nullptr};
}
