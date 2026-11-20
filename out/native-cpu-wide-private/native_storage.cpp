// Read-only storage audit of the optimized SIMD16 entries before codegen.
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/Config/llvm-config.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Verifier.h"
#include "llvm/Passes/PassBuilder.h"
#include "llvm/Plugins/PassPlugin.h"
#include "llvm/Support/ErrorHandling.h"

#if LLVM_VERSION_MAJOR != 23 || LLVM_VERSION_MINOR != 1 || LLVM_VERSION_PATCH != 2
#error "SIMD16 native storage audit requires LLVM 23.1.2."
#endif

using namespace llvm;

namespace {
constexpr unsigned Lanes = 16;
constexpr unsigned ColumnBytes = Lanes * 4;

[[noreturn]] void reject(const Twine &Why, const Instruction *I = nullptr) {
  errs() << "SIMD16 native storage: " << Why << '\n';
  if (I) {
    I->print(errs());
    errs() << '\n';
  }
  report_fatal_error("unadmitted SIMD16 native storage", false);
}

bool words(Type *T) {
  auto *Vector = dyn_cast<FixedVectorType>(T);
  return Vector && Vector->getNumElements() == Lanes && Vector->getElementType()->isIntegerTy(32);
}

bool vector32(Type *T) {
  auto *Vector = dyn_cast<FixedVectorType>(T);
  return Vector && Vector->getNumElements() == Lanes &&
      (Vector->getElementType()->isIntegerTy(32) || Vector->getElementType()->isFloatTy());
}

#include "endpoint_local.inc"

uint64_t extent(AllocaInst *A, const DataLayout &DL) {
  const auto *Count = dyn_cast<ConstantInt>(A->getArraySize());
  TypeSize Size = DL.getTypeAllocSize(A->getAllocatedType());
  if (!Count || !Count->isOne() || A->getAddressSpace() != 0 || Size.isScalable())
    reject("dynamic or nonlocal allocation", A);
  return Size.getFixedValue();
}

unsigned endpointLocal(AllocaInst *A, Function &F, const DataLayout &DL) {
  StringRef Name = A->getName();
  const bool Root = Name == "root_1.i";
  if (!Root && Name != "programBase_3.i" && Name != "count_2.i" && Name != "layers_2.i" &&
      Name != "prefix_5.i" && Name != "component_4.i") reject("unknown Endpoint residual", A);
  unsigned Size = Root ? 5 * ColumnBytes : ColumnBytes;
  if (extent(A, DL) != Size || A->getAlign() != Align(ColumnBytes)) reject("Endpoint local size/alignment differs", A);
  if (!Root) {
    if (!words(A->getAllocatedType())) reject("Endpoint selector vector differs", A);
  } else {
    auto *S = dyn_cast<StructType>(A->getAllocatedType());
    if (!S || S->isOpaque() || S->isPacked() || S->getNumElements() != 5 || !words(S->getElementType(4)) ||
        DL.getStructLayout(S)->getElementOffset(4) != 4 * ColumnBytes)
      reject("Endpoint projection root differs", A);
    for (unsigned K = 0; K != 4; ++K)
      if (!vector32(S->getElementType(K)) || !cast<FixedVectorType>(S->getElementType(K))->getElementType()->isFloatTy() ||
          DL.getStructLayout(S)->getElementOffset(K) != K * ColumnBytes) reject("Endpoint root field differs", A);
  }
  local(A, Size, F, DL);
  return Size;
}

struct StorageAudit : PassInfoMixin<StorageAudit> {
  PreservedAnalyses run(Module &M, ModuleAnalysisManager &) {
    Function *Entry = nullptr;
    for (Function &F : M) {
      if (F.isDeclaration()) continue;
      if (Entry) reject("more than one native entry");
      Entry = &F;
    }
    if (!Entry || Entry->getCallingConv() != CallingConv::C || Entry->arg_size() != 6 || Entry->isVarArg() ||
        Entry->getVisibility() != GlobalValue::HiddenVisibility || !Entry->getReturnType()->isVoidTy() ||
        !Entry->hasFnAttribute(Attribute::StrictFP) || Entry->getUWTableKind() != UWTableKind::Async)
      reject("sole hidden strict C entry with async unwind is required");
    const auto Name = Entry->getName();
    const bool Transport = Name == "sirius_wide_strict_transport_packet";
    const bool Endpoint = Name == "sirius_wide_strict_endpoint_packet";
    const bool Dopri = Name == "sirius_wide_strict_dp_packet";
    if (!Transport && !Endpoint && !Dopri) reject("unknown packet entry");
    const unsigned Words = Transport ? 4699 : Endpoint ? 3829 : 4149;
    unsigned Staging = 0, Residual = 0, ResidualBytes = 0;
    const DataLayout &DL = M.getDataLayout();
    for (Instruction &I : instructions(*Entry)) {
      auto *A = dyn_cast<AllocaInst>(&I);
      if (!A) continue;
      if (A->getName() == "staged_output") {
        auto *T = dyn_cast<ArrayType>(A->getAllocatedType());
        if (!T || T->getNumElements() != Words * Lanes || !T->getElementType()->isIntegerTy(32) ||
            A->getAlign() != Align(ColumnBytes) || extent(A, DL) != Words * ColumnBytes)
          reject("transactional staging ABI differs", A);
        ++Staging;
        continue;
      }
      if (Transport && A->getName() == "wide.coefficient.private8000") {
        auto *T = dyn_cast<ArrayType>(A->getAllocatedType());
        if (!T || T->getNumElements() != 125 || !words(T->getElementType()) || extent(A, DL) != 8000 ||
            A->getAlign() != Align(ColumnBytes)) reject("Transport coefficient cache ABI differs", A);
        ResidualBytes += 8000;
      } else if (Endpoint) {
        ResidualBytes += endpointLocal(A, *Entry, DL);
      } else {
        reject("unexpected private allocation after storage promotion", A);
      }
      ++Residual;
    }
    if (Staging != 1 || verifyModule(M, &errs())) reject("native storage audit did not confirm the entry");
    errs() << "SIMD16 native storage: " << Name << " staging_bytes=" << Words * ColumnBytes
           << " residual_allocas=" << Residual << " residual_bytes=" << ResidualBytes << '\n';
    return PreservedAnalyses::all();
  }
};
} // namespace

extern "C" LLVM_ATTRIBUTE_WEAK PassPluginLibraryInfo llvmGetPassPluginInfo() {
  return {LLVM_PLUGIN_API_VERSION, "SiriusWideNativeStorage", LLVM_VERSION_STRING, [](PassBuilder &Builder) {
    Builder.registerPipelineParsingCallback(
      [](StringRef Name, ModulePassManager &Passes, ArrayRef<PassBuilder::PipelineElement> Inner) {
        if (Name != "sirius-wide-native-storage" || !Inner.empty()) return false;
        Passes.addPass(StorageAudit());
        return true;
      });
  }, nullptr};
}
