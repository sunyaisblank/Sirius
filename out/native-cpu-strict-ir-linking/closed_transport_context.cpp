// Generated for one exact closed transport value module. No math definition is linked.
#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/MapVector.h"
#include "llvm/IR/CFG.h"
#include <optional>
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
constexpr StringLiteral Input = "/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-canonical-closed-transport/transport.value.bc";
constexpr StringLiteral InputSHA = "7d810098cf634829e40b4bfc3d4758ac276378a85114369b1cd95ca4ae0630a4";
constexpr uint64_t InputBytes = 16012504, ContextBytes = 92192;
constexpr StringLiteral EntryName = "sirius_closed_transport_packet";
constexpr StringLiteral MaskedName = "sirius_closed_transport_packet___un_3C_unu_3E_unuun_3C_unu_3E_unuunuunu";
constexpr StringLiteral MathPrefix = "sirius_local_value_";

[[noreturn]] void reject(const Twine &Why, const Instruction *I = nullptr) {
  errs() << "sirius-closed-transport-context: " << Why << '\n';
  if (I) { errs() << I->getFunction()->getName() << ": "; I->print(errs()); errs() << '\n'; }
  report_fatal_error("closed transport context normalization refused", false);
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
  errs() << "sirius-closed-transport-context: input-sha256=" << InputSHA << " module-fingerprint=" << Expected << '\n';
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
      !S->getElementType(0)->isPointerTy() || !words(S->getElementType(3))) return false;
  for (auto Pair : {std::pair<unsigned,unsigned>{1,2295}, {2,459}, {4,125}}) {
    auto *A = dyn_cast<ArrayType>(S->getElementType(Pair.first));
    if (!A || A->getNumElements() != Pair.second || !words(A->getElementType())) return false;
  }
  auto *Layout = DL.getStructLayout(S);
  return DL.getTypeAllocSize(S).getFixedValue() == ContextBytes &&
    Layout->getElementOffset(0) == 0 &&
    Layout->getElementOffset(1) == 32 &&
    Layout->getElementOffset(2) == 73472 &&
    Layout->getElementOffset(3) == 88160 &&
    Layout->getElementOffset(4) == 88192;
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

// Four exact source counters; this is not a general dynamic-GEP permission.
constexpr StringLiteral DomainKey = "sirius.closed.coefficient.source.counter";
using Assignment = DenseMap<const Value *, APInt>;
std::optional<APInt> integerValue(Value *V, const Assignment &A, unsigned Depth = 0) {
  if (Depth > 40) return std::nullopt;
  if (auto *C = dyn_cast<ConstantInt>(V)) return C->getValue();
  auto It = A.find(V); if (It != A.end()) return It->second;
  auto *I = dyn_cast<Instruction>(V); if (!I || !I->getType()->isIntegerTy()) return std::nullopt;
  auto Get = [&](Value *X) { return integerValue(X, A, Depth + 1); };
  if (auto *C = dyn_cast<CastInst>(I)) {
    auto X = Get(C->getOperand(0)); if (!X) return std::nullopt;
    unsigned Width = C->getType()->getIntegerBitWidth();
    switch (C->getOpcode()) {
      case Instruction::Trunc: return X->trunc(Width);
      case Instruction::ZExt: return X->zext(Width);
      case Instruction::SExt: return X->sext(Width);
      default: return std::nullopt;
    }
  }
  if (auto *S = dyn_cast<SelectInst>(I)) {
    auto C = Get(S->getCondition()); if (!C) return std::nullopt;
    return Get(S->getOperand(C->isZero() ? 2 : 1));
  }
  auto *B = dyn_cast<BinaryOperator>(I);
  auto *C = dyn_cast<ICmpInst>(I);
  if (!B && !C) return std::nullopt;
  auto L = Get(I->getOperand(0)), R = Get(I->getOperand(1));
  if (!L || !R || L->getBitWidth() != R->getBitWidth()) return std::nullopt;
  if (C) return APInt(1, ICmpInst::compare(*L, *R, C->getPredicate()));
  switch (B->getOpcode()) {
    case Instruction::Add: return *L + *R;
    case Instruction::Sub: return *L - *R;
    case Instruction::Mul: return *L * *R;
    case Instruction::And: return *L & *R;
    case Instruction::Or: return *L | *R;
    case Instruction::Xor: return *L ^ *R;
    case Instruction::Shl:
      if (R->uge(L->getBitWidth())) return std::nullopt; return L->shl(R->getZExtValue());
    case Instruction::LShr:
      if (R->uge(L->getBitWidth())) return std::nullopt; return L->lshr(R->getZExtValue());
    default: return std::nullopt;
  }
}
std::optional<APInt> pointerOffset(Value *P, Value *Base, const Assignment &A,
                                  const DataLayout &DL, unsigned Depth = 0) {
  if (P == Base) return APInt(64, 0);
  if (Depth > 20) return std::nullopt;
  auto *G = dyn_cast<GEPOperator>(P); if (!G || G->getPointerAddressSpace() != 0) return std::nullopt;
  auto Offset = pointerOffset(G->getPointerOperand(), Base, A, DL, Depth + 1);
  if (!Offset) return std::nullopt;
  SmallMapVector<Value *, APInt, 4> Variables; APInt Constant(64, 0);
  if (!G->collectOffset(DL, 64, Variables, Constant)) return std::nullopt;
  *Offset += Constant;
  for (auto &Term : Variables) {
    auto Index = integerValue(Term.first, A); if (!Index) return std::nullopt;
    *Offset += Index->sextOrTrunc(64) * Term.second;
  }
  return Offset;
}
PHINode *namedCounter(Function &F, StringRef Name) {
  for (Instruction &I : instructions(F)) if (I.getName() == Name) {
    auto *P = dyn_cast<PHINode>(&I); if (!P) reject("source counter is not a PHI", &I); return P;
  }
  reject("exact source counter absent: " + Name);
}
void proveCounter(PHINode *P, uint64_t Start, Value *Limit, unsigned Identity) {
  if (!P->getType()->isIntegerTy(64) || P->getNumIncomingValues() != 2) reject("source counter ABI differs", P);
  unsigned Init = 2;
  for (unsigned J = 0; J < 2; ++J) if (auto *C = dyn_cast<ConstantInt>(P->getIncomingValue(J)))
    if (C->getZExtValue() == Start) Init = J;
  if (Init == 2) reject("source counter initial value differs", P);
  auto *Step = dyn_cast<BinaryOperator>(P->getIncomingValue(1 - Init));
  auto *One = Step ? dyn_cast<ConstantInt>(Step->getOperand(1)) : nullptr;
  if (!Step || Step->getOpcode() != Instruction::Add || Step->getOperand(0) != P || !One || !One->isOne() ||
      !Step->hasNoUnsignedWrap() || !Step->hasNoSignedWrap() ||
      Step->getParent() != P->getIncomingBlock(1 - Init)) reject("source counter unit recurrence differs", P);
  auto *Branch = dyn_cast<BranchInst>(Step->getParent()->getTerminator());
  auto *Test = Branch && Branch->isConditional() ? dyn_cast<ICmpInst>(Branch->getCondition()) : nullptr;
  if (!Test || Test->getPredicate() != ICmpInst::ICMP_EQ || Test->getOperand(0) != Step ||
      Test->getOperand(1) != Limit || Branch->getSuccessor(1) != P->getParent() ||
      Branch->getSuccessor(0) == P->getParent()) reject("source counter exact terminating branch differs", P);
  P->setMetadata(DomainKey, MDNode::get(P->getContext(),
      {ConstantAsMetadata::get(ConstantInt::get(Type::getInt32Ty(P->getContext()), Identity))}));
}
void markCoefficientDomains(Module &M) {
  Function *F = M.getFunction("ComputeOwnedRow_0___s_5B_vyRowID_5D_unuun_3C_s_5B_vyKernelContext_0_5D__3E_");
  if (!F || F->isDeclaration()) reject("pinned coefficient source function absent");
  PHINode *Cache = namedCounter(*F, "indvars.iv2007"), *Stage = namedCounter(*F, "indvars.iv2034"),
          *Inner = namedCounter(*F, "indvars.iv2014"), *Last = namedCounter(*F, "indvars.iv2042");
  Type *I64 = Type::getInt64Ty(M.getContext());
  proveCounter(Cache, 0, ConstantInt::get(I64, 25), 0);
  proveCounter(Stage, 0, ConstantInt::get(I64, 7), 1);
  proveCounter(Inner, 1, Stage, 2);
  proveCounter(Last, 2, ConstantInt::get(I64, 7), 3);
  DominatorTree DT(*F); bool InnerGuard = false;
  for (BasicBlock &BB : *F) if (auto *B = dyn_cast<BranchInst>(BB.getTerminator())) if (B->isConditional())
    if (auto *C = dyn_cast<ICmpInst>(B->getCondition())) if (C->getPredicate() == ICmpInst::ICMP_UGT &&
        C->getOperand(0) == Stage && C->getOperand(1) == ConstantInt::get(I64, 1) &&
        DT.dominates(BasicBlockEdge(&BB, B->getSuccessor(0)), Inner->getParent())) InnerGuard = true;
  if (!InnerGuard) reject("inner coefficient loop lacks dominating stage>1 edge");
}
AllocaInst *ownedAllocation(Function &F, const DataLayout &DL) {
  AllocaInst *A = nullptr;
  for (Instruction &I : instructions(F)) if (auto *X = dyn_cast<AllocaInst>(&I))
    if (contextType(X->getAllocatedType(), DL)) { if (A) reject("multiple exact contexts", X); A = X; }
  if (!A || !isa<ConstantInt>(A->getArraySize()) || !cast<ConstantInt>(A->getArraySize())->isOne() ||
      A->getAddressSpace() != 0 || A->getAlign() != Align(32)) reject("single exact32-aligned owned context missing");
  return A;
}
unsigned expandReset(Function &F, AllocaInst *A, const DataLayout &DL) {
  BasicBlock *Loop = nullptr; PHINode *IV = nullptr; BasicBlock *Pre = nullptr, *Exit = nullptr;
  SmallVector<IntrinsicInst *, 6> Stores; Value *Mask = nullptr;
  // Exactly one source reset block has a self-edge and six zero masked stores.
  for (BasicBlock &BB : F) {
    auto *B = dyn_cast<BranchInst>(BB.getTerminator());
    if (!B || !B->isConditional() || B->getSuccessor(1) != &BB || B->getSuccessor(0) == &BB) continue;
    SmallVector<IntrinsicInst *, 6> Candidate;
    for (Instruction &I : BB) if (auto *II = dyn_cast<IntrinsicInst>(&I))
      if (II->getIntrinsicID() == Intrinsic::masked_store) Candidate.push_back(II);
    if (Candidate.size() != 6) continue;
    if (Loop) reject("more than one possible exact reset self-loop");
    Loop = &BB; Exit = B->getSuccessor(0); Stores = Candidate;
  }
  if (!Loop || !Exit->phis().empty()) reject("exact reset loop/PHI-free exit absent");
  auto *B = cast<BranchInst>(Loop->getTerminator());
  auto *Test = dyn_cast<ICmpInst>(B->getCondition());
  auto *Step = Test ? dyn_cast<BinaryOperator>(Test->getOperand(0)) : nullptr;
  auto *Bound = Test ? dyn_cast<ConstantInt>(Test->getOperand(1)) : nullptr;
  auto *Three = Step ? dyn_cast<ConstantInt>(Step->getOperand(1)) : nullptr;
  IV = Step ? dyn_cast<PHINode>(Step->getOperand(0)) : nullptr;
  if (!Test || Test->getPredicate() != ICmpInst::ICMP_EQ || !Bound || Bound->getZExtValue() != 459 ||
      !Step || Step->getOpcode() != Instruction::Add || !Three || Three->getZExtValue() != 3 ||
      !Step->hasNoUnsignedWrap() || !Step->hasNoSignedWrap() || !IV || IV->getParent() != Loop ||
      !IV->getType()->isIntegerTy(64) || IV->getNumIncomingValues() != 2) reject("reset IV0/+3/exit459 shape differs");
  for (unsigned J = 0; J < 2; ++J) {
    if (IV->getIncomingBlock(J) == Loop) { if (IV->getIncomingValue(J) != Step) reject("reset backedge differs"); }
    else { Pre = IV->getIncomingBlock(J); auto *Z = dyn_cast<ConstantInt>(IV->getIncomingValue(J));
      if (!Z || !Z->isZero()) reject("reset initial IV differs"); }
  }
  auto *PB = Pre ? dyn_cast<BranchInst>(Pre->getTerminator()) : nullptr;
  if (!PB || !PB->isUnconditional() || PB->getSuccessor(0) != Loop || pred_size(Loop) != 2)
    reject("reset preheader/predecessors differ");
  Mask = Stores.front()->getArgOperand(2);
  if (auto *MI = dyn_cast<Instruction>(Mask)) if (MI->getParent() == Loop) reject("reset mask is not invariant", MI);
  for (Instruction &I : *Loop) {
    if (&I == IV || &I == B || &I == Test) {}
    else if (auto *G = dyn_cast<GetElementPtrInst>(&I)) { if (G->getPointerAddressSpace() != 0) reject("reset GEP address space differs", G); }
    else if (auto *Op = dyn_cast<BinaryOperator>(&I)) {
      if (!Op->getType()->isIntegerTy(64) || (Op->getOpcode() != Instruction::Add && Op->getOpcode() != Instruction::Mul))
        reject("reset includes unsupported integer operation", Op);
    } else if (auto *II = dyn_cast<IntrinsicInst>(&I)) {
      auto *Zero = II->arg_size() ? dyn_cast<Constant>(II->getArgOperand(0)) : nullptr;
      if (II->getIntrinsicID() != Intrinsic::masked_store || II->arg_size() != 3 || !Zero || !Zero->isNullValue() ||
          !words(Zero->getType()) || II->getArgOperand(2) != Mask) reject("reset includes nonzero/nonlocal/other effect", II);
    } else reject("reset includes an unknown instruction", &I);
    if (!I.getType()->isVoidTy()) for (User *U : I.users()) {
      auto *UI = dyn_cast<Instruction>(U); if (!UI || UI->getParent() != Loop) reject("reset value has live-out consumer", &I);
    }
  }
  // Check every original dynamic address for all153 iterations before replacing
  // only this memory-initialization bookkeeping with its identical918 stores.
  for (unsigned K = 0; K < 153; ++K) {
    Assignment Values; Values.try_emplace(IV, APInt(64, 3 * K));
    for (unsigned J = 0; J < 6; ++J) {
      auto Offset = pointerOffset(Stores[J]->getArgOperand(1), A, Values, DL);
      uint64_t Expected = J % 2 ? 73472 + 32 * (3 * K + J / 2) : 160 + 160 * (3 * K + J / 2);
      if (!Offset || *Offset != APInt(64, Expected)) reject("reset affine field address/order differs", Stores[J]);
    }
  }
  IRBuilder<> At(PB);
  for (unsigned K = 0; K < 459; ++K) for (unsigned J = 0; J < 2; ++J) {
    auto *Copy = cast<IntrinsicInst>(Stores[J]->clone());
    uint64_t Offset = J ? 73472 + 32 * K : 160 + 160 * K;
    Copy->setArgOperand(1, At.CreateGEP(Type::getInt8Ty(F.getContext()), A, At.getInt64(Offset)));
    Copy->insertBefore(PB->getIterator());
  }
  BranchInst::Create(Exit, PB); PB->eraseFromParent();
  Loop->dropAllReferences(); Loop->eraseFromParent();
  return 918;
}
struct CounterDomain { PHINode *P; unsigned Low, High; };
SmallVector<CounterDomain, 4> coefficientDomains(Function &F) {
  SmallVector<CounterDomain, 4> Out; bool Found[4] = {};
  constexpr unsigned Low[] = {0,0,1,2}, High[] = {24,6,6,6};
  for (Instruction &I : instructions(F)) if (MDNode *MD = I.getMetadata(DomainKey)) {
    if (MD->getNumOperands() != 1) reject("counter annotation shape differs", &I);
    auto *CM = dyn_cast<ConstantAsMetadata>(MD->getOperand(0));
    auto *C = CM ? dyn_cast<ConstantInt>(CM->getValue()) : nullptr; auto *P = dyn_cast<PHINode>(&I);
    if (!C || !P || C->getZExtValue() > 3 || Found[C->getZExtValue()]) reject("source counter cloned/missing/unknown", &I);
    unsigned K = C->getZExtValue(); Found[K] = true; Out.push_back({P, Low[K], High[K]});
  }
  if (Out.size() != 4) reject("four source-bound coefficient counter domains required");
  return Out;
}
void proveCoefficientAccess(Instruction *I, Value *P, AllocaInst *A, const DataLayout &DL,
                            DominatorTree &DT, ArrayRef<CounterDomain> Domains) {
  Assignment Empty;
  if (auto Offset = pointerOffset(P, A, Empty, DL)) {
    if (Offset->isNegative() || Offset->ugt(3968) || (Offset->getZExtValue() % 32))
      reject("constant coefficient address exceeds exact initialized4000B field", I);
    return;
  }
  SmallVector<std::pair<Value *, bool>, 8> Conditions;
  for (BasicBlock &BB : *I->getFunction()) if (auto *B = dyn_cast<BranchInst>(BB.getTerminator()))
    if (B->isConditional()) for (unsigned Edge = 0; Edge < 2; ++Edge)
      if (DT.dominates(BasicBlockEdge(&BB, B->getSuccessor(Edge)), I->getParent()))
        Conditions.push_back({B->getCondition(), Edge == 0});
  unsigned Accepted = 0; Assignment Values;
  // Finite over-approximation of only four proven scalar source loop domains.
  // Dominating scalar edges remove impossible stage0/stage1 read cases.
  auto Visit = [&](auto &&Self, unsigned D) -> void {
    if (D != Domains.size()) {
      auto Domain = Domains[D];
      for (unsigned K = Domain.Low; K <= Domain.High; ++K) {
        Values[Domain.P] = APInt(64, K); Self(Self, D + 1);
      }
      return;
    }
    for (auto Constraint : Conditions) if (auto V = integerValue(Constraint.first, Values))
      if ((!V->isZero()) != Constraint.second) return;
    auto Offset = pointerOffset(P, A, Values, DL);
    if (!Offset || Offset->isNegative() || Offset->ugt(3968) || (Offset->getZExtValue() % 32))
      reject("coefficient access not proved inside exact initialized4000B field", I);
    ++Accepted;
  };
  Visit(Visit, 0); if (!Accepted) reject("coefficient access has no proved source-domain execution", I);
}
void isolateCoefficient(Function &F, AllocaInst *ContextA, const DataLayout &DL) {
  auto Domains = coefficientDomains(F);
  MemSetInst *ContextZero = nullptr;
  for (Instruction &I : instructions(F)) if (auto *S = dyn_cast<MemSetInst>(&I)) {
    int64_t O = 0; auto *N = dyn_cast<ConstantInt>(S->getLength()); auto *Z = dyn_cast<ConstantInt>(S->getValue());
    if (GetPointerBaseWithConstantOffset(S->getRawDest(), O, DL) == ContextA && !O && N &&
        N->getZExtValue() == ContextBytes && Z && Z->isZero() && !S->isVolatile()) {
      if (ContextZero) reject("more than one complete context zero"); ContextZero = S;
    }
  }
  if (!ContextZero) reject("complete coefficient source initialization absent");
  SmallVector<GetElementPtrInst *, 32> Roots;
  for (Instruction &I : instructions(F)) if (auto *G = dyn_cast<GetElementPtrInst>(&I)) {
    int64_t O = 0;
    if (GetPointerBaseWithConstantOffset(G, O, DL) == ContextA && O >= 88192 && O < 92192) Roots.push_back(G);
  }
  if (Roots.empty()) reject("coefficient field pointer roots absent");
  IRBuilder<> Entry(&*F.getEntryBlock().getFirstInsertionPt());
  auto *Type = ArrayType::get(FixedVectorType::get(Entry.getInt32Ty(), 8), 125);
  auto *A = Entry.CreateAlloca(Type, nullptr, "closed.coefficient.private4000"); A->setAlignment(Align(64));
  IRBuilder<> Init(ContextZero); auto *Zero = Init.CreateMemSet(A, Init.getInt8(0), 4000, Align(64));
  for (GetElementPtrInst *G : Roots) {
    int64_t O = 0;
    // Earlier parent replacement can already move this constant child.
    if (GetPointerBaseWithConstantOffset(G, O, DL) != ContextA) continue;
    IRBuilder<> At(G); Value *New = At.CreateGEP(At.getInt8Ty(), A, At.getInt64(O - 88192));
    G->replaceAllUsesWith(New);
  }
  DominatorTree DT(F); SmallVector<Value *, 128> Pending{A}; SmallPtrSet<Value *, 32> Seen;
  unsigned Accesses = 0;
  while (!Pending.empty()) {
    Value *P = Pending.pop_back_val(); if (!Seen.insert(P).second) continue;
    for (User *U : P->users()) {
      if (auto *G = dyn_cast<GetElementPtrInst>(U)) {
        if (G->getPointerOperand() != P) reject("coefficient pointer used as index", G); Pending.push_back(G); continue;
      }
      Instruction *I = dyn_cast<Instruction>(U); if (!I) reject("coefficient pointer escapes outside instructions");
      if (I == Zero) continue;
      if (auto *L = dyn_cast<LoadInst>(I)) {
        if (L->getPointerOperand() != P || L->isAtomic() || L->isVolatile() || !words(L->getType()))
          reject("unsupported coefficient load", L);
      } else if (auto *S = dyn_cast<StoreInst>(I)) {
        if (S->getPointerOperand() != P || S->isAtomic() || S->isVolatile() || !words(S->getValueOperand()->getType()))
          reject("unsupported coefficient store/capture", S);
      } else if (auto *II = dyn_cast<IntrinsicInst>(I)) {
        bool Load = II->getIntrinsicID() == Intrinsic::masked_load;
        if ((!Load && II->getIntrinsicID() != Intrinsic::masked_store) || II->arg_size() != 3 ||
            II->getArgOperand(Load ? 0 : 1) != P ||
            !words(Load ? II->getType() : II->getArgOperand(0)->getType())) reject("coefficient callback/unknown intrinsic", II);
        auto *Mask = dyn_cast<FixedVectorType>(II->getArgOperand(Load ? 1 : 2)->getType());
        if (!Mask || Mask->getNumElements() != 8 || !Mask->getElementType()->isIntegerTy(1)) reject("coefficient mask ABI differs", II);
      } else reject("coefficient pointer escapes to unsupported operation", I);
      if (!DT.dominates(ContextZero, I) || !DT.dominates(Zero, I)) reject("full coefficient zero does not dominate memory", I);
      proveCoefficientAccess(I, P, A, DL, DT, Domains); ++Accesses;
    }
  }
  if (!Accesses) reject("isolated coefficient accesses absent");
  // Remove only unused original address calculations; no memory/FP operation.
  bool Changed;
  do { Changed = false;
    for (auto It = inst_begin(F), End = inst_end(F); It != End;) {
      Instruction *I = &*It++;
      if (isa<GetElementPtrInst>(I) && I->use_empty()) { I->eraseFromParent(); Changed = true; }
    }
  } while (Changed);
  errs() << "sirius-closed-transport-context: isolated_coefficient_bytes=4000 bounded_accesses=" << Accesses
         << " original_cache_RK_loops_and_masked_operations_unchanged=1 reset_zero_stores=918\n";
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
    markCoefficientDomains(M);
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
        if (!T || T->getNumElements() != 37592 ||
            !T->getElementType()->isIntegerTy(32) || A->getAlign() != Align(64) ||
            !isa<ConstantInt>(A->getArraySize()) || !cast<ConstantInt>(A->getArraySize())->isOne())
          reject("original transactional staging allocation differs", A);
        ++Staging;
      }
    if (Staging != 1) reject("one exact original transactional staging allocation required");
    Inventory Before = fpInventory(M);
    AllocaInst *Owned = ownedAllocation(*Entry, M.getDataLayout());
    if (expandReset(*Entry, Owned, M.getDataLayout()) != 918) reject("exact reset expansion incomplete");
    isolateCoefficient(*Entry, Owned, M.getDataLayout());
    if (fpInventory(M) != Before) reject("reset/cache isolation changed post-inline FP/math inventory");
    Context Object = inspect(*Entry, M.getDataLayout());
    unsigned Count = Object.Masked.size(); normalize(Object);
    if (fpInventory(M) != Before) reject("context normalization changed post-inline FP/math inventory");
    if (verifyModule(M, &errs())) reject("LLVM rejected context normalization output");
    errs() << "sirius-closed-transport-context: raw_definitions_inlined=" << RawDefinitions
           << " owned_contexts=1 context_bytes=92192 dominated_memory_accesses=" << Object.Memory.size()
           << " local_masked_operations_normalized=" << Count << " postinline_FP_inventory_unchanged=1\n";
    return PreservedAnalyses::none();
  }
};
} // namespace
extern "C" LLVM_ATTRIBUTE_WEAK PassPluginLibraryInfo llvmGetPassPluginInfo() {
  return {LLVM_PLUGIN_API_VERSION, "SiriusClosedTransportContext", LLVM_VERSION_STRING,
    [](PassBuilder &PB) { PB.registerPipelineParsingCallback(
      [](StringRef Name, ModulePassManager &MPM, ArrayRef<PassBuilder::PipelineElement>) {
        if (Name != "sirius-closed-transport-context") return false; MPM.addPass(ClosedContext()); return true;
      }); }, nullptr};
}
