// Exact DP-only layer projection. Source contract: layer_projection_contract.json.
#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/ADT/StringExtras.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/Config/llvm-config.h"
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
#include "llvm/Support/SHA256.h"
#include "llvm/Support/SourceMgr.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Transforms/Utils/Cloning.h"
#include <map>
#include <set>
#include <string>

#if LLVM_VERSION_MAJOR != 23 || LLVM_VERSION_MINOR != 1 || LLVM_VERSION_PATCH != 2
#error "Requires the pinned LLVM23.1.2 API."
#endif
using namespace llvm;
namespace {
constexpr StringLiteral Input = "/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-local-effects/dp.cleanup.bc";
constexpr StringLiteral InputSHA = "a2f52cebde91547a09f0693776f7e017ec1f447893895b1e360e4cd460b5bf6e";
constexpr uint64_t InputBytes = 8644532, ContextBytes = 120832;
constexpr unsigned ExpectedLayers = 41, MaximumSlots = 1153;

[[noreturn]] void reject(const Twine &Why, const Instruction *I = nullptr) {
  errs() << "sirius-project-dp: " << Why << '\n';
  if (I) { errs() << I->getFunction()->getName() << ": "; I->print(errs()); errs() << '\n'; }
  report_fatal_error("exact DP layer projection refused", false);
}
std::string digest(StringRef S) { SHA256 H; H.update(S); return toHex(H.final(), true); }
std::string fingerprint(const Module &M) {
  std::string S; raw_string_ostream O(S); M.print(O, nullptr); O.flush(); return digest(S);
}
void requirePinned(const Module &M) {
  auto Buffer = MemoryBuffer::getFile(Input);
  if (!Buffer || (*Buffer)->getBufferSize() != InputBytes || digest((*Buffer)->getBuffer()) != InputSHA)
    reject("original complete input pin differs");
  LLVMContext C; SMDiagnostic D;
  auto Reference = parseIR((*Buffer)->getMemBufferRef(), D, C);
  if (!Reference) reject("cannot independently parse original input");
  Reference->setModuleIdentifier(M.getModuleIdentifier());
  auto Expected = fingerprint(*Reference);
  if (fingerprint(M) != Expected) reject("full canonical module fingerprint differs");
  errs() << "sirius-project-dp: input-sha256=" << InputSHA << " module-fingerprint=" << Expected << '\n';
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
             I.getOpcode() == Instruction::FCmp || I.getOpcode() == Instruction::FPTrunc ||
             I.getOpcode() == Instruction::FPExt) reject("ordinary FP operation", &I);
  }
  return Out;
}
unsigned valueCalls(const Module &M) {
  unsigned N = 0;
  for (const Function &F : M) for (const Instruction &I : instructions(F))
    if (auto *C = dyn_cast<CallInst>(&I)) if (auto *Target = C->getCalledFunction())
      N += Target->getName().starts_with("sirius_local_value_");
  return N;
}
bool vectorWord(Type *T) {
  auto *V = dyn_cast<FixedVectorType>(T);
  return V && V->getNumElements() == 8 &&
    (V->getElementType()->isIntegerTy(32) || V->getElementType()->isFloatTy());
}

// This is a readonly closure of an existing global-view pointer, not a new
// observer. Only the exact ReadOriginal helper and existing read gathers may
// receive it. Stores, captures, integerization and unknown calls fail closed.
void readonlyView(Value *Start, SmallPtrSetImpl<Value *> &Seen) {
  SmallVector<Value *, 16> Pending{Start};
  while (!Pending.empty()) {
    Value *P = Pending.pop_back_val(); if (!Seen.insert(P).second) continue;
    for (User *U : P->users()) {
      if (auto *G = dyn_cast<GetElementPtrInst>(U)) {
        if (G->getPointerOperand() != P) reject("view pointer used as index", G);
        Pending.push_back(G); continue;
      }
      if (auto *L = dyn_cast<LoadInst>(U)) {
        if (L->getPointerOperand() != P || L->isVolatile() || L->isAtomic()) reject("nonordinary view read", L);
        if (L->getType()->isPointerTy()) Pending.push_back(L);
        continue;
      }
      if (auto *C = dyn_cast<CallInst>(U)) {
        Function *Callee = C->getCalledFunction();
        if (!Callee) reject("indirect view escape", C);
        StringRef Name = Callee->getName();
        if (Name.starts_with("llvm.x86.avx512.mask.gather") && C->arg_size() == 5 && C->getArgOperand(1) == P) continue;
        if (Name.starts_with("ReadOriginal_0___") && !Callee->isDeclaration()) {
          unsigned Matches = 0;
          for (unsigned N = 0; N < C->arg_size(); ++N) if (C->getArgOperand(N) == P) {
            ++Matches; readonlyView(Callee->getArg(N), Seen);
          }
          if (Matches == 1) continue;
        }
        reject("unknown call receives readonly global view", C);
      }
      reject("global view write/capture/unsupported use", dyn_cast<Instruction>(U));
    }
  }
}

struct Slot { bool Read = false, Write = false; };
struct Access { Instruction *I; uint64_t Offset; };
struct Layer {
  Function *F; Argument *Context; unsigned Number;
  std::map<uint64_t, Slot> Slots;
  SmallVector<Access, 64> Accesses;
  SmallVector<LoadInst *, 4> ViewLoads;
  SmallVector<Instruction *, 64> PointerNodes;
};
Layer inspect(Function &F, unsigned Number, const DataLayout &DL) {
  if (F.isDeclaration() || F.isVarArg() || !F.getReturnType()->isVoidTy() || F.hasPersonalityFn() ||
      F.getCallingConv() != CallingConv::Fast || !F.hasFnAttribute(Attribute::StrictFP)) reject("layer ABI/strict scope differs");
  Argument *Context = nullptr;
  for (Argument &A : F.args()) if (A.getType()->isPointerTy()) {
    if (Context || !A.getName().starts_with("kernelContext_")) reject("unknown layer pointer formal");
    Context = &A;
  }
  if (!Context) reject("layer has no exact context formal");
  for (User *U : F.users()) {
    auto *C = dyn_cast<CallInst>(U);
    if (!C || C->getCalledFunction() != &F) reject("layer address escapes", dyn_cast<Instruction>(U));
  }
  Layer L{&F, Context, Number, {}, {}, {}, {}};
  SmallVector<Value *, 64> Pending{Context}; SmallPtrSet<Value *, 32> Seen, ViewSeen;
  while (!Pending.empty()) {
    Value *P = Pending.pop_back_val(); if (!Seen.insert(P).second) continue;
    int64_t Offset = 0;
    if (GetPointerBaseWithConstantOffset(P, Offset, DL) != Context || Offset < 0 || uint64_t(Offset) >= ContextBytes)
      reject("dynamic/out-of-range context pointer", dyn_cast<Instruction>(P));
    for (User *U : P->users()) {
      if (auto *G = dyn_cast<GetElementPtrInst>(U)) {
        if (G->getPointerOperand() != P) reject("context used as GEP index", G);
        L.PointerNodes.push_back(G); Pending.push_back(G); continue;
      }
      if (auto *Load = dyn_cast<LoadInst>(U)) {
        if (Load->getPointerOperand() != P || Load->isAtomic() || Load->isVolatile()) reject("nonordinary context load", Load);
        if (Offset == 0 && Load->getType()->isPointerTy()) {
          readonlyView(Load, ViewSeen); L.ViewLoads.push_back(Load); continue;
        }
      }
      bool Read = false, Write = false; Type *T = nullptr; Instruction *I = dyn_cast<Instruction>(U);
      if (auto *Load = dyn_cast<LoadInst>(U)) { Read = true; T = Load->getType(); }
      else if (auto *Store = dyn_cast<StoreInst>(U)) {
        if (Store->getPointerOperand() != P || Store->isAtomic() || Store->isVolatile()) reject("context pointer stored/escaped", Store);
        Write = true; T = Store->getValueOperand()->getType();
      } else if (auto *II = dyn_cast<IntrinsicInst>(U)) {
        if (II->getIntrinsicID() == Intrinsic::masked_load && II->arg_size() == 3 && II->getArgOperand(0) == P) {
          Read = true; T = II->getType();
        } else if (II->getIntrinsicID() == Intrinsic::masked_store && II->arg_size() == 3 && II->getArgOperand(1) == P) {
          Write = true; T = II->getArgOperand(0)->getType();
        } else reject("unknown context-derived intrinsic/call", II);
      } else reject("context pointer escape or unsupported use", I);
      if (!T || !vectorWord(T) || Offset < 32 || Offset % 32 || uint64_t(Offset) + 32 > ContextBytes)
        reject("context field is not an exact complete initialized eight-word slot", I);
      Slot &S = L.Slots[Offset]; S.Read |= Read; S.Write |= Write;
      L.Accesses.push_back({I, uint64_t(Offset)});
    }
  }
  if (L.Slots.empty() || L.Slots.size() > MaximumSlots) reject("projection exceeds exact64-op layer bound");
  return L;
}
Value *ioSlot(IRBuilder<> &B, Value *Base, unsigned N) {
  return B.CreateConstGEP1_64(B.getInt8Ty(), Base, uint64_t(N) * 32, "project.io.slot");
}

void project(Layer &L) {
  Function &F = *L.F; Module &M = *F.getParent(); LLVMContext &C = M.getContext();
  Type *V = FixedVectorType::get(Type::getInt32Ty(C), 8);
  SmallVector<Type *, 5> Types{PointerType::get(C, 0)};
  bool View = !L.ViewLoads.empty(); if (View) Types.push_back(PointerType::get(C, 0));
  for (Argument &A : F.args()) if (&A != L.Context) Types.push_back(A.getType());
  std::string Name = "sirius_projected_layer_" + std::to_string(L.Number);
  if (M.getNamedValue(Name)) reject("projected-layer name collision");
  Function *Core = Function::Create(FunctionType::get(Type::getVoidTy(C), Types, false), GlobalValue::InternalLinkage, Name, M);
  Core->setCallingConv(CallingConv::Fast); Core->setDSOLocal(true);
  BasicBlock *Entry = BasicBlock::Create(C, "project.initialize", Core); IRBuilder<> B(Entry);
  std::map<uint64_t, AllocaInst *> Locals; unsigned N = 0;
  for (auto &[Offset, S] : L.Slots) {
    auto *A = B.CreateAlloca(V, nullptr, "project.local." + std::to_string(Offset)); A->setAlignment(Align(64));
    Locals[Offset] = A;
    B.CreateAlignedStore(B.CreateAlignedLoad(V, ioSlot(B, Core->getArg(0), N++), Align(32)), A, Align(64));
  }
  ValueToValueMapTy Map; unsigned Arg = 1 + View;
  Map[L.Context] = PoisonValue::get(L.Context->getType()); // No such pointer may survive below.
  for (Argument &A : F.args()) if (&A != L.Context) Map[&A] = Core->getArg(Arg++);
  SmallVector<ReturnInst *, 4> Returns;
  CloneFunctionInto(Core, &F, Map, CloneFunctionChangeType::LocalChangesOnly, Returns);
  B.CreateBr(cast<BasicBlock>(Map[&F.getEntryBlock()]));
  Core->setAttributes(AttributeList::get(C, F.getAttributes().getFnAttrs(), AttributeSet{}, {}));
  Core->removeFnAttr(Attribute::AlwaysInline); Core->addFnAttr(Attribute::NoInline); Core->addFnAttr(Attribute::StrictFP);
  for (LoadInst *Old : L.ViewLoads) {
    auto *Load = cast<LoadInst>(Map[Old]); Load->replaceAllUsesWith(Core->getArg(1)); Load->eraseFromParent();
  }
  unsigned Normalized = 0;
  for (Access A : L.Accesses) {
    Instruction *I = cast<Instruction>(Map[A.I]); Value *P = Locals.at(A.Offset);
    if (auto *Load = dyn_cast<LoadInst>(I)) { if (Load->getAlign() > Align(64)) reject("unreviewed load alignment", Load); Load->setOperand(0, P); }
    else if (auto *Store = dyn_cast<StoreInst>(I)) { if (Store->getAlign() > Align(64)) reject("unreviewed store alignment", Store); Store->setOperand(1, P); }
    else {
      auto *II = cast<IntrinsicInst>(I); bool Load = II->getIntrinsicID() == Intrinsic::masked_load;
      unsigned PI = Load ? 0 : 1; Align Alignment = II->getParamAlign(PI).valueOrOne();
      if (Alignment > Align(64)) reject("unreviewed masked alignment", II);
      IRBuilder<> At(II); Type *T = Load ? II->getType() : II->getArgOperand(0)->getType();
      Value *Old = At.CreateAlignedLoad(T, P, Alignment, "project.old.initialized");
      Value *Mask = II->getArgOperand(Load ? 1 : 2);
      if (Load) II->replaceAllUsesWith(At.CreateSelect(Mask, Old, II->getArgOperand(2), "project.masked.read"));
      else At.CreateAlignedStore(At.CreateSelect(Mask, II->getArgOperand(0), Old, "project.masked.merge"), P, Alignment);
      II->eraseFromParent(); ++Normalized;
    }
  }
  // Every old context-address computation is now dead, including its poison
  // placeholder. Unknown uses were rejected before cloning, never waived.
  bool Changed;
  do {
    Changed = false;
    for (Instruction *Old : L.PointerNodes) if (auto *I = dyn_cast_or_null<Instruction>(Map[Old]))
      if (I->getParent() && I->use_empty()) { I->eraseFromParent(); Changed = true; }
  } while (Changed);
  for (Instruction *Old : L.PointerNodes) if (Map[Old]) reject("old context-derived pointer remains in core");
  for (auto &[Offset, A] : Locals) for (User *U : A->users()) {
    if (auto *Load = dyn_cast<LoadInst>(U)) { if (Load->getPointerOperand() == A) continue; }
    if (auto *Store = dyn_cast<StoreInst>(U)) { if (Store->getPointerOperand() == A) continue; }
    reject("projected local escaped its loads/stores", dyn_cast<Instruction>(U));
  }
  unsigned Written = 0;
  for (auto &[Offset, S] : L.Slots) Written += S.Write;
  for (ReturnInst *Ret : Returns) {
    IRBuilder<> Out(Ret); N = 0;
    for (auto &[Offset, S] : L.Slots) {
      if (S.Write) Out.CreateAlignedStore(Out.CreateAlignedLoad(V, Locals.at(Offset), Align(64)), ioSlot(Out, Core->getArg(0), N), Align(32));
      ++N;
    }
  }
  // The old-name adapter becomes tiny and alwaysinline. No aggregate return
  // or hundreds-vector native ABI is introduced: one compact IO pointer.
  F.deleteBody(); F.setLinkage(GlobalValue::InternalLinkage); F.setVisibility(GlobalValue::DefaultVisibility); F.setDSOLocal(true);
  F.removeFnAttr(Attribute::NoInline); F.addFnAttr(Attribute::AlwaysInline); F.addFnAttr(Attribute::StrictFP);
  IRBuilder<> Adapter(BasicBlock::Create(C, "project.adapter", &F));
  auto *IO = Adapter.CreateAlloca(ArrayType::get(V, L.Slots.size()), nullptr, "project.io"); IO->setAlignment(Align(64));
  Adapter.CreateLifetimeStart(IO); N = 0;
  for (auto &[Offset, S] : L.Slots) {
    Value *P = Adapter.CreateConstGEP1_64(Adapter.getInt8Ty(), L.Context, Offset);
    Adapter.CreateAlignedStore(Adapter.CreateAlignedLoad(V, P, Align(1)), ioSlot(Adapter, IO, N++), Align(32));
  }
  SmallVector<Value *, 5> Args{IO};
  if (View) Args.push_back(Adapter.CreateAlignedLoad(PointerType::get(C, 0), L.Context, Align(1), "project.readonly.view"));
  for (Argument &A : F.args()) if (&A != L.Context) Args.push_back(&A);
  auto *Call = Adapter.CreateCall(Core, Args); Call->setCallingConv(CallingConv::Fast); Call->addFnAttr(Attribute::StrictFP);
  N = 0;
  for (auto &[Offset, S] : L.Slots) {
    if (S.Write) Adapter.CreateAlignedStore(Adapter.CreateAlignedLoad(V, ioSlot(Adapter, IO, N), Align(32)),
        Adapter.CreateConstGEP1_64(Adapter.getInt8Ty(), L.Context, Offset), Align(1));
    ++N;
  }
  Adapter.CreateLifetimeEnd(IO); Adapter.CreateRetVoid();
  errs() << "sirius-project-dp: layer=" << L.Number << " slots=" << L.Slots.size()
         << " compact_bytes=" << L.Slots.size() * 32 << " written_slots=" << Written
         << " local_masked_operations_normalized=" << Normalized << " whole_context_in_core=0\n";
}

struct ProjectDP : PassInfoMixin<ProjectDP> {
  PreservedAnalyses run(Module &M, ModuleAnalysisManager &) {
    requirePinned(M); auto Before = fpInventory(M); SmallVector<Layer, 41> Layers;
    if (valueCalls(M) != 1581) reject("original shared math call inventory differs");
    for (Function &F : M) if (F.getName().starts_with("EvaluateStaticLayer_0_")) {
      StringRef Tail = F.getName().drop_front(StringRef("EvaluateStaticLayer_0_").size());
      auto Pair = Tail.split('_'); unsigned Number = 0;
      if (Pair.first.getAsInteger(10, Number) || Number < 4 || Number > 44) reject("unexpected surviving DP layer identity");
      Layers.push_back(inspect(F, Number, M.getDataLayout()));
    }
    if (Layers.size() != ExpectedLayers) reject("surviving layer count differs");
    std::set<unsigned> Numbers; for (auto &L : Layers) Numbers.insert(L.Number);
    if (Numbers.size() != ExpectedLayers) reject("duplicate layer identity");
    for (Layer &L : Layers) project(L);
    for (StringRef Suffix : {StringRef(""), StringRef("___un_3C_unu_3E_unuun_3C_unu_3E_unuunuunu")}) {
      auto *F = M.getFunction((Twine("sirius_local_dp_packet") + Suffix).str());
      std::string Name = (Twine("sirius_projected_dp_packet") + Suffix).str();
      if (!F || F->isDeclaration() || M.getNamedValue(Name)) reject("raw export missing/collision");
      F->setName(Name);
    }
    if (fpInventory(M) != Before) reject("projection changed complete strict FP inventory");
    if (valueCalls(M) != 1581) reject("projection changed shared math calls");
    if (verifyModule(M, &errs())) reject("LLVM rejected projected module");
    errs() << "sirius-project-dp: projected_layers=41 old_public_objects_changed=0 FP_inventory_unchanged=1\n";
    return PreservedAnalyses::none();
  }
};
} // namespace
extern "C" LLVM_ATTRIBUTE_WEAK PassPluginLibraryInfo llvmGetPassPluginInfo() {
  return {LLVM_PLUGIN_API_VERSION, "SiriusDPLayerProjection", LLVM_VERSION_STRING,
    [](PassBuilder &PB) { PB.registerPipelineParsingCallback(
      [](StringRef Name, ModulePassManager &MPM, ArrayRef<PassBuilder::PipelineElement>) {
        if (Name != "sirius-project-dp") return false; MPM.addPass(ProjectDP()); return true;
      }); }, nullptr};
}
