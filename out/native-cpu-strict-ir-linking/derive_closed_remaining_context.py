"""Derive one pinned T/E context plugin source; never run a compiler or IR pass."""
from pathlib import Path
import ast
import hashlib
import json
import sys

R = Path(__file__).resolve().parents[2]
V = Path(__file__).resolve().parent
W = R / 'out/native-cpu-canonical-closed-remaining'
PARENT = V / 'derive_closed_dp_context.py'
PARENT_PIN = dict(bytes=22855, sha256='ffd94f50726aa20895ce91d6820dee67cb722f142bbd4798556d52b31ca0814d')
SOURCE = V / 'canonical_private_remaining_contract.json'
SOURCE_PIN = dict(bytes=58311, sha256='f2f9b0ec0fb0930485eaf0016d219327044f65602d22a9696b47e9450f98e1cc')
PROFILES = {
    'transport': dict(bytes=16012504, sha256='7d810098cf634829e40b4bfc3d4758ac276378a85114369b1cd95ca4ae0630a4',
        context_bytes=92192, fields=['ptr', '2295*v8i32', '459*v8i32', 'v8i32', '125*v8i32'],
        offsets=[0,32,73472,88160,88192], staged_words=37592, pointer_calls=2698, wrappers=8,
        value_contract=dict(bytes=12730, sha256='6854d06caef0a6b2204a9f7bff02bb7a95049ecc3bff1a89bbba3446ac7a1593')),
    'endpoint': dict(bytes=22601688, sha256='61321c43b452dcc98d0c3ed1a45623012437f2aedf6c47d3e8ebe799b3d02aaa',
        context_bytes=117792, fields=['ptr', '7*v8i32', 'v8i32', '3060*v8i32', '612*v8i32'],
        offsets=[0,32,256,288,98208], staged_words=30632, pointer_calls=3954, wrappers=9,
        value_contract=dict(bytes=12734, sha256='08ee08a9fbf42e3c34a50cd4d0c19f9069accd4d2a16eeb806995c7910f6f50f')),
}

TRANSPORT_CPP = r'''
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
'''


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def doc(path):
    return json.loads(path.read_text())


def terminal(path):
    row = doc(path)
    assert row['status'] == 'terminal' and row['passed'] is True and row['child_exit'] == 0
    assert row['source_unchanged'] and row['whole_inputs_unchanged']
    assert not row['remaining_owned_processes'] and not row['cleanup_errors']
    assert row['observed_births'] and all(x['absent'] for x in row['observed_births'])
    assert path.with_name('inputs-before.json').read_bytes() == path.with_name('inputs-after.json').read_bytes()


def replace_once(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new, 1)


def template(stage):
    assert pin(PARENT) == PARENT_PIN
    tree = ast.parse(PARENT.read_text())
    text = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == 'CPP' for t in n.targets))
    profile = PROFILES[stage]
    text = text.replace('closed DP', 'closed '+stage).replace('closed-dp', 'closed-'+stage)
    text = text.replace('sirius_closed_dp_packet', 'sirius_closed_'+stage+'_packet')
    text = text.replace('SiriusClosedDPContext', 'SiriusClosed'+stage.title()+'Context')
    text = text.replace('context_bytes=120832', 'context_bytes='+str(profile['context_bytes']))
    text = replace_once(text, 'ContextBytes = 120832', 'ContextBytes = '+str(profile['context_bytes']))
    start = text.index('bool contextType('); end = text.index('\nvoid bounds(', start)
    if stage == 'transport':
        types = '''
  if (!S || S->isOpaque() || S->isPacked() || S->getNumElements() != 5 ||
      !S->getElementType(0)->isPointerTy() || !words(S->getElementType(3))) return false;
  for (auto Pair : {std::pair<unsigned,unsigned>{1,2295}, {2,459}, {4,125}}) {
    auto *A = dyn_cast<ArrayType>(S->getElementType(Pair.first));
    if (!A || A->getNumElements() != Pair.second || !words(A->getElementType())) return false;
  }
'''
    else:
        types = '''
  if (!S || S->isOpaque() || S->isPacked() || S->getNumElements() != 5 ||
      !S->getElementType(0)->isPointerTy() || !words(S->getElementType(2))) return false;
  for (auto Pair : {std::pair<unsigned,unsigned>{1,7}, {3,3060}, {4,612}}) {
    auto *A = dyn_cast<ArrayType>(S->getElementType(Pair.first));
    if (!A || A->getNumElements() != Pair.second || !words(A->getElementType())) return false;
  }
'''
    layout = '\n'.join('    Layout->getElementOffset('+str(i)+') == '+str(offset)+(' &&' if i<4 else ';')
                       for i,offset in enumerate(profile['offsets']))
    new = '''bool contextType(Type *T, const DataLayout &DL) {
  auto *S = dyn_cast<StructType>(T);'''+types+'''  auto *Layout = DL.getStructLayout(S);
  return DL.getTypeAllocSize(S).getFixedValue() == ContextBytes &&
'''+layout+'\n}'
    text = text[:start]+new+text[end:]
    # The callback observer must unwind this final giant ordinary function.
    text = replace_once(text, '    unsigned MathDeclarations = 0;',
        '    Entry->setUWTableKind(UWTableKind::Async);\n    unsigned MathDeclarations = 0;')
    # Keep the original transactional staging allocation exact; it is not a
    # removable Context and is intentionally excluded from memory normalization.
    text = replace_once(text, '    Inventory Before = fpInventory(M);', '''    unsigned Staging = 0;
    for (Instruction &I : instructions(Entry)) if (auto *A = dyn_cast<AllocaInst>(&I))
      if (A->getName() == "staged_output") {
        auto *T = dyn_cast<ArrayType>(A->getAllocatedType());
        if (!T || T->getNumElements() != '''+str(profile['staged_words'])+''' ||
            !T->getElementType()->isIntegerTy(32) || A->getAlign() != Align(64) ||
            !isa<ConstantInt>(A->getArraySize()) || !cast<ConstantInt>(A->getArraySize())->isOne())
          reject("original transactional staging allocation differs", A);
        ++Staging;
      }
    if (Staging != 1) reject("one exact original transactional staging allocation required");
    Inventory Before = fpInventory(M);''')
    if stage == 'transport':
        text = text.replace('#include "llvm/ADT/DenseMap.h"',
            '#include "llvm/ADT/DenseMap.h"\n#include "llvm/ADT/MapVector.h"\n#include "llvm/IR/CFG.h"\n#include <optional>')
        text = replace_once(text, 'struct ClosedContext : PassInfoMixin<ClosedContext> {',
                            TRANSPORT_CPP+'\nstruct ClosedContext : PassInfoMixin<ClosedContext> {')
        text = replace_once(text, '    unsigned MathDeclarations = 0;',
                            '    markCoefficientDomains(M);\n    unsigned MathDeclarations = 0;')
        text = replace_once(text, '    Context Object = inspect(*Entry, M.getDataLayout());', '''    AllocaInst *Owned = ownedAllocation(*Entry, M.getDataLayout());
    if (expandReset(*Entry, Owned, M.getDataLayout()) != 918) reject("exact reset expansion incomplete");
    isolateCoefficient(*Entry, Owned, M.getDataLayout());
    if (fpInventory(M) != Before) reject("reset/cache isolation changed post-inline FP/math inventory");
    Context Object = inspect(*Entry, M.getDataLayout());''')
    return text


def main():
    assert len(sys.argv) == 3 and sys.argv[1] in PROFILES and sys.argv[2] == 'derive'
    stage = sys.argv[1]; profile = PROFILES[stage]
    B = R / ('bin/linux-gcc/native-cpu-canonical-closed-'+stage)
    input_path = B / (stage+'.value.bc')
    expected = dict(path=str(input_path.relative_to(R)), **{k:profile[k] for k in ('bytes','sha256')})
    assert pin(SOURCE) == SOURCE_PIN
    source = doc(SOURCE); proof = source['source_proofs'][stage]
    for key in ('source_inverse_to_publication_exact','source_inverse_to_AS_exact',
                'unchanged_outside_export_binding_and_EvaluateProgram','generic_fallback_removed_only_from_new_private_body',
                'canonical_tuple_records_match_exact_tables','complete_context_memset_and_full8_physical_packing_unchanged',
                'supported_callers_lane0','mixed_families_supported','homogeneous_all_active_guard_absent'):
        assert proof[key] is True, key
    assert pin(B/'ispc-source'/ (stage+'.ispc')) == source['derived_sources'][stage]
    owner = W / ('value-'+stage) / 'owner.json'; terminal(owner)
    result_path = W / (stage+'-value-result.json'); result = doc(result_path)
    assert result['phase'] == 'value-'+stage and result['output'] == expected
    assert pin(input_path) == {k:profile[k] for k in ('bytes','sha256')}
    assert result['accepted_artifact_owner'] == str(owner.relative_to(R))
    assert sum(result['pointer_calls'].values()) == profile['pointer_calls']
    assert sum(result['value_calls'].values()) == profile['wrappers']
    snap = doc(owner.with_name('inputs-before.json'))
    strict_input = result['input']
    raw = snap[strict_input['path']]
    assert {k:raw[k] for k in ('bytes','sha256')} == {k:strict_input[k] for k in ('bytes','sha256')}
    assert Path(raw['resolved_path']).resolve() == (R/strict_input['path']).resolve()
    assert pin(R/strict_input['path']) == {k:strict_input[k] for k in ('bytes','sha256')}
    value_path = V / ('closed_'+stage+'_value_contract.json'); value = doc(value_path)
    assert pin(value_path) == profile['value_contract'] == result['source_value_contract']
    assert value['phase'] == 'value' and value['stage'] == stage and value['source_inverse_exact']
    assert value['value_math_prefix'] == 'sirius_local_value_'
    assert pin(R/value['output']['path']) == {k:value['output'][k] for k in ('bytes','sha256')}
    inputs = dict(value['source_inputs']); inputs.update(source['source_inputs'])
    for name,identity in inputs.items(): assert pin(R/name) == identity, name
    for path in (Path(__file__).resolve(), PARENT, SOURCE, value_path, result_path, owner,
                 owner.with_name('inputs-before.json'),owner.with_name('inputs-after.json'),
                 input_path, B/'ispc-source'/(stage+'.ispc')):
        inputs[str(path.relative_to(R))] = pin(path)
    text = template(stage)
    for token,replacement in (('@INPUT@',json.dumps(str(input_path))),('@SHA@',json.dumps(profile['sha256'])),
                              ('@BYTES@',str(profile['bytes']))):
        text = replace_once(text,token,replacement)
    assert not any(token in text for token in ('@INPUT@','@SHA@','@BYTES@'))
    target = V / ('closed_'+stage+'_context.cpp'); contract = target.with_name('closed_'+stage+'_context_contract.json')
    assert not target.exists() and not contract.exists()
    with target.open('x') as stream: stream.write(text)
    record = dict(stage=stage,phase='context',source_only=True,pass_name='sirius-closed-'+stage+'-context',
        input=expected,output=dict(path=str(target.relative_to(R)),**pin(target)),source_inputs=inputs,
        private_scope=value['private_scope'],exact_context_layout=dict(bytes=profile['context_bytes'],
            fields=profile['fields'],offsets=profile['offsets'],all_physical_lanes=8),
        staged_output=dict(words=profile['staged_words'],bytes=profile['staged_words']*4,alignment=64,preserved=True),
        allowed_residual_alloca=(dict(role='coefficient',elements=125,vector='v8i32',bytes=4000,alignment=64)
                                if stage=='transport' else None),
        ordinary_entry='sirius_closed_'+stage+'_packet',ordinary_entry_async_uwtable=True,
        source_proof=dict(mixed_family_masks_and_order_preserved=True,live_endpoint11241_guard_preserved=True,
            all_root_validation_before_publication=True,original_allocation_and_transaction_contract_preserved=True,
            raw_only_inlining=True,external_math9_cc92_inaccessibleRW_unchanged=True,
            normalization_only_full_initialized_inbounds_private_context=True,
            postinline_FP_math_inventory_unchanged_by_memory_rewrite=True),
        reset=(dict(source_shape='one IV0/+3/exit459 self-loop, sixzero maskedstores, invariantmask, noFP/callback/liveouts',
                    action='same459 ordered valid/tag pairs,918 stores at exactconstantoffsets',
                    coefficient_RK_and_FP_loops_unchanged=True) if stage=='transport' else None),
        coefficient=(dict(source_field=[88192,92192],private_zero_initialized_bytes=4000,
             counter_domains={'cache':[0,24],'stage':[0,6],'inner':[1,6],'last':[2,6]},
             proof='Exact counted recurrences/terminating branches plus stage>1 loop admission; finite integer-address enumeration constrained by dominating scalar edges; rejectunknown addresses/captures. Keep all coefficient masked operations and source loops.')
                     if stage=='transport' else None),
        initialization='Full exact context zero dominates every normalized memory operation; all eight physical lanes initialized. Stop at GlobalParams pointer load, preserve original masks/passthroughs/alignment, no freezing.',
        pipeline='Pre-DCE masked duplicate; RAW AlwaysInliner once and DCE; preserve postinline FP/math inventory across private-memory-only changes; set ordinary AsyncUWTable. Root focused cleanup must eliminate main Context while preserving exact transactional staging and only optional4000B coefficient cache.',
        qualification='Source-only. No compile/native/runtime or numerical acceptance transfer; actual structural/native/FDE gates required before any genuine-job cohort.')
    with contract.open('x') as stream: json.dump(record,stream,indent=2); stream.write('\n')
    print(json.dumps(dict(stage=stage,source_only=True,compiler_run=False,output=record['output'],
        contract=dict(path=str(contract.relative_to(R)),**pin(contract)),pass_name=record['pass_name'])))


if __name__ == '__main__':
    main()
