// Exact four-input shared value ABI. No public/GCC object is changed.
// See shared_value_contract.json; reviewed local-bank logic is reused verbatim.
#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/STLExtras.h"
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/ADT/StringExtras.h"
#include "llvm/Analysis/ConstantFolding.h"
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
#include "llvm/Transforms/Utils/Cloning.h"
#include <array>
#include <map>
#include <memory>

#if LLVM_VERSION_MAJOR != 23 || LLVM_VERSION_MINOR != 1 || LLVM_VERSION_PATCH != 2
#error "Requires the pinned LLVM 23.1.2 API."
#endif
using namespace llvm;
namespace {
constexpr StringLiteral Directory = "/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-strict-ir-linking/";
constexpr StringLiteral Prefix = "sirius_scope_v4_";
constexpr StringLiteral SharedPrefix = "sirius_shared_value_";
constexpr StringLiteral MaskedSuffix = "___un_3C_unu_3E_unuun_3C_unu_3E_unuunuunu";
struct Profile {
  const char *Pass, *Input; uint64_t Bytes; const char *SHA, *Export;
  unsigned Mode, Calls, TagBanks;
};
constexpr std::array<Profile,4> Profiles{{
 {"sirius-shared-dp","dopri_phase.strict.sha256.bc",10451600,
  "bbdf247215b4f5f36f9872f55d202db879dbbb9d9301dbecdc863c27fa0c0979","sirius_scope_dp_packet",0,1581,1575},
 {"sirius-shared-transport","transport.strict.bc",16091084,
  "6e81247821df090e23b1859f726e756559e5e6b4487d95a6e71233bbb8b60d1f","sirius_scope_transport_packet",1,2701,2682},
 {"sirius-shared-endpoint","endpoint.strict.bc",22512336,
  "4dd9a944b89573a654fceddbeb47dd1f7c9c1833e3c5a4ef9282f9e9dfb6ed74","sirius_scope_endpoint_packet",2,3959,3924},
 {"sirius-shared-math","private-math.bc",344448,
  "f9108d507569ce44730aface0011c44247e181c7bdadf0b8466dbd3ab907c982","",3,0,0}}};
struct Recipe { const char *Suffix; bool Tagged; std::array<unsigned,3> Calls; };
constexpr std::array<Recipe,9> Recipes{{
 {"dense_add8",false,{3,8,13}}, {"dense_multiply8",false,{1,7,13}},
 {"dense_divide8",false,{1,3,7}}, {"dense_sqrt8",false,{1,1,2}},
 {"tagged_add8",true,{585,947,1607}}, {"tagged_subtract8",true,{284,10,197}},
 {"tagged_multiply8",true,{661,1482,1963}}, {"tagged_divide8",true,{45,240,154}},
 {"tagged_sqrt8",true,{0,3,3}}}};

[[noreturn]] void reject(const Twine &Why, const Instruction *I=nullptr) {
  errs()<<"sirius-shared-value: "<<Why<<'\n';
  if(I) { errs()<<"in "<<I->getFunction()->getName()<<": "; I->print(errs()); errs()<<'\n'; }
  report_fatal_error("exact private value ABI refused",false);
}
std::string digest(StringRef S) { SHA256 H;H.update(S);return toHex(H.final(),true); }
std::string fingerprint(const Module &M) {
  std::string S;raw_string_ostream O(S);M.print(O,nullptr);O.flush();return digest(S);
}
void requirePinned(const Module &M,const Profile &P) {
  auto Buffer=MemoryBuffer::getFile((Directory+P.Input).str());
  if(!Buffer || (*Buffer)->getBufferSize()!=P.Bytes || digest((*Buffer)->getBuffer())!=P.SHA)
    reject("exact input bitcode pin absent/changed");
  LLVMContext C;SMDiagnostic D;
  auto Reference=parseIR((*Buffer)->getMemBufferRef(),D,C);
  if(!Reference)reject("cannot independently parse pinned input");
  Reference->setModuleIdentifier(M.getModuleIdentifier());
  std::string Expected=fingerprint(*Reference);
  if(fingerprint(M)!=Expected)reject("full canonical module fingerprint differs");
  errs()<<P.Pass<<": input-sha256="<<P.SHA<<" module-fingerprint="<<Expected<<'\n';
}
bool words(Type *T) {
  auto *V=dyn_cast<FixedVectorType>(T);
  return V && V->getNumElements()==8 && V->getElementType()->isIntegerTy(32);
}
bool maskType(Type *T) {
  auto *V=dyn_cast<FixedVectorType>(T);
  return V && V->getNumElements()==8 && V->getElementType()->isIntegerTy(1);
}
uint64_t sizeOf(Type *T,const DataLayout &DL) {
  auto S=DL.getTypeStoreSize(T);if(S.isScalable())reject("scalable memory access");
  return S.getFixedValue();
}
int64_t offset(Value *P,AllocaInst *A,const DataLayout &DL) {
  int64_t O=0;
  if(GetPointerBaseWithConstantOffset(P,O,DL)!=A)reject("nonconstant/local-bank pointer");
  return O;
}
void bounds(Value *P,AllocaInst *A,uint64_t Bytes,uint64_t Width,const DataLayout &DL) {
  int64_t O=offset(P,A,DL);
  if(O<0 || uint64_t(O)>Bytes || Width>Bytes-uint64_t(O))reject("local-bank access outside exact allocation");
}

struct Bank {
  AllocaInst *A;
  uint64_t Bytes;
  CallInst *Math;
  SmallVector<StoreInst *,18> Stores;
  SmallVector<MemSetInst *,2> Zeros;
  SmallVector<IntrinsicInst *,12> Masked;
  SmallVector<Instruction *,2> Starts,Ends;
};

// This walks only one named fixed-size compiler-private allocation's uses.
// Captures, pointer arithmetic outside constant GEPs and external memory are
// rejected. The complete input fingerprint supplies the source proof domain.
Bank inspectBank(AllocaInst *A,unsigned Count,CallInst *Math,const DataLayout &DL) {
  auto *AT=dyn_cast<ArrayType>(A->getAllocatedType());
  auto *N=dyn_cast<ConstantInt>(A->getArraySize());
  if(!AT || AT->getNumElements()!=Count || !AT->getElementType()->isIntegerTy(32) ||
      !N || !N->isOne() || A->getAddressSpace()!=0)reject("unexpected private allocation",A);
  if(!A->getName().starts_with(Count==120 ? "Sirius" : "SiriusPackedTags"))
    reject("unnamed/unreviewed compiler-private allocation",A);
  Bank B{A,uint64_t(Count)*4,Math,{},{},{},{},{}};
  SmallVector<Value *,32> Pending{A};SmallPtrSet<Value *,32> Seen;
  while(!Pending.empty()) {
    Value *P=Pending.pop_back_val();if(!Seen.insert(P).second)continue;
    for(User *U:P->users()) {
      if(auto *G=dyn_cast<GetElementPtrInst>(U)) {
        if(G->getPointerOperand()!=P)reject("private pointer used as GEP index",G);
        bounds(G,A,B.Bytes,0,DL);Pending.push_back(G);continue;
      }
      if(auto *L=dyn_cast<LoadInst>(U)) {
        if(L->getPointerOperand()!=P || L->isVolatile() || L->isAtomic())reject("nonordinary private load",L);
        bounds(P,A,B.Bytes,sizeOf(L->getType(),DL),DL);continue;
      }
      if(auto *S=dyn_cast<StoreInst>(U)) {
        if(S->getPointerOperand()!=P || S->isVolatile() || S->isAtomic())reject("private pointer escape/store",S);
        bounds(P,A,B.Bytes,sizeOf(S->getValueOperand()->getType(),DL),DL);
        B.Stores.push_back(S);continue;
      }
      auto *C=dyn_cast<CallInst>(U);
      if(C && C==Math) {
        if(P!=A)reject("math receives interior private pointer",C);
        continue;
      }
      auto *II=dyn_cast<IntrinsicInst>(U);
      if(!II)reject("unsupported private pointer use",dyn_cast<Instruction>(U));
      if(auto *Z=dyn_cast<MemSetInst>(II)) {
        auto *Length=dyn_cast<ConstantInt>(Z->getLength());
        auto *Byte=dyn_cast<ConstantInt>(Z->getValue());
        if(Z->getRawDest()!=P || Z->isVolatile() || !Length || !Byte || !Byte->isZero())
          reject("nonzero/unknown private initializer",Z);
        bounds(P,A,B.Bytes,Length->getZExtValue(),DL);B.Zeros.push_back(Z);continue;
      }
      switch(II->getIntrinsicID()) {
      case Intrinsic::lifetime_start: B.Starts.push_back(II);break;
      case Intrinsic::lifetime_end: B.Ends.push_back(II);break;
      case Intrinsic::masked_load:
        if(II->arg_size()!=3 || II->getArgOperand(0)!=P || !words(II->getType()) ||
            !maskType(II->getArgOperand(1)->getType()))reject("unknown LLVM23 private masked load",II);
        bounds(P,A,B.Bytes,32,DL);B.Masked.push_back(II);break;
      case Intrinsic::masked_store:
        if(II->arg_size()!=3 || II->getArgOperand(1)!=P || !words(II->getArgOperand(0)->getType()) ||
            !maskType(II->getArgOperand(2)->getType()))reject("unknown LLVM23 private masked store",II);
        bounds(P,A,B.Bytes,32,DL);B.Masked.push_back(II);break;
      default:reject("unsupported call using original packed allocation",II);
      }
    }
  }
  return B;
}

bool zeroCovers(const Bank &B,unsigned Column,Instruction *At,DominatorTree &DT,const DataLayout &DL) {
  uint64_t Start=uint64_t(Column)*32;
  for(auto *Z:B.Zeros) {
    uint64_t O=uint64_t(offset(Z->getRawDest(),B.A,DL));
    uint64_t L=cast<ConstantInt>(Z->getLength())->getZExtValue();
    if(O<=Start && L>=32 && Start-O<=L-32 && DT.dominates(Z,At))return true;
  }
  return false;
}
bool columnInitialized(const Bank &B,unsigned Column,Instruction *At,DominatorTree &DT,
                       const DataLayout &DL,bool ZeroOnly=false) {
  if(zeroCovers(B,Column,At,DT,DL))return true;
  for(auto *S:B.Stores) {
    if(!words(S->getValueOperand()->getType()) || offset(S->getPointerOperand(),B.A,DL)!=int64_t(Column)*32 ||
        !DT.dominates(S,At))continue;
    if(!ZeroOnly)return true;
    auto *C=dyn_cast<Constant>(S->getValueOperand());if(C && C->isNullValue())return true;
  }
  return false;
}
void checkInitialization(const Bank &B,unsigned Destination,DominatorTree &DT,const DataLayout &DL) {
  for(unsigned C=0;C<B.Bytes/32;++C)
    if(!columnInitialized(B,C,B.Math,DT,DL))reject("private bank lacks dominating complete initialization",B.Math);
  unsigned Width=B.Bytes==480 ? 5:1;
  for(auto *S:B.Stores) {
    int64_t O=offset(S->getPointerOperand(),B.A,DL);
    if(O>=int64_t(Destination*Width)*32 && O<int64_t((Destination+1)*Width)*32 && DT.dominates(S,B.Math)) {
      auto *C=dyn_cast<Constant>(S->getValueOperand());
      if(!C || !C->isNullValue())reject("nonzero private destination overwrite before math",S);
    }
  }
  for(auto *II:B.Masked)if(II->getIntrinsicID()==Intrinsic::masked_store) {
    int64_t O=offset(II->getArgOperand(1),B.A,DL);
    if(O>=int64_t(Destination*Width)*32 && O<int64_t((Destination+1)*Width)*32 && DT.dominates(II,B.Math))
      reject("masked private destination overwrite before math",II);
  }
  for(unsigned C=0;C<Width;++C)
    if(!columnInitialized(B,Destination*Width+C,B.Math,DT,DL,true))
      reject("private old destination is not established zero",B.Math);
  if(!B.Starts.empty() && !llvm::any_of(B.Starts,[&](auto *S){return DT.dominates(S,B.Math);}))
    reject("private lifetime start does not dominate its math call",B.Math);
  if(llvm::any_of(B.Ends,[&](auto *E){return DT.dominates(E,B.Math);}))
    reject("private lifetime already ended before its math call",B.Math);
  for(auto *II:B.Masked) {
    unsigned P=II->getIntrinsicID()==Intrinsic::masked_load ? 0:1;
    unsigned Column=unsigned(offset(II->getArgOperand(P),B.A,DL)/32);
    if(offset(II->getArgOperand(P),B.A,DL)%32 || !columnInitialized(B,Column,II,DT,DL))
      reject("masked normalization would read an uninitialized full vector",II);
  }
}

void normalize(Bank &B) {
  for(auto *II:B.Masked) {
    bool Load=II->getIntrinsicID()==Intrinsic::masked_load;
    unsigned P=Load?0:1;
    IRBuilder<> Builder(II);
    Value *Pointer=II->getArgOperand(P);
    Type *T=Load?II->getType():II->getArgOperand(0)->getType();
    Align Alignment=II->getParamAlign(P).valueOrOne();
    auto *Old=Builder.CreateAlignedLoad(T,Pointer,Alignment,"value.old.private");
    Value *Mask=II->getArgOperand(Load?1:2);
    if(Load) {
      Value *Selected=Builder.CreateSelect(Mask,Old,II->getArgOperand(2),"value.masked.private");
      II->replaceAllUsesWith(Selected);
    } else {
      Value *Selected=Builder.CreateSelect(Mask,II->getArgOperand(0),Old,"value.merge.private");
      Builder.CreateAlignedStore(Selected,Pointer,Alignment);
    }
    II->eraseFromParent();
  }
}

Value *slot(IRBuilder<> &B,Value *P,unsigned Column) {
  return B.CreateConstGEP1_64(B.getInt8Ty(),P,uint64_t(Column)*32,"value.slot");
}
void foldFixedIndices(Function &F,const DataLayout &DL) {
  bool Changed;
  do {
    Changed=false;
    for(auto &BB:F)for(auto I=BB.begin(),End=BB.end();I!=End;) {
      Instruction *Inst=&*I++;
      if(auto *Phi=dyn_cast<PHINode>(Inst);Phi && Phi->getType()->isPointerTy()) {
        int64_t CommonOffset=0;
        auto *Base=dyn_cast<AllocaInst>(GetPointerBaseWithConstantOffset(
            Phi->getIncomingValue(0),CommonOffset,DL));
        bool Same=Base;
        for(Value *Incoming:Phi->incoming_values()) {
          int64_t ThisOffset=0;
          Same &= GetPointerBaseWithConstantOffset(Incoming,ThisOffset,DL)==Base && ThisOffset==CommonOffset;
        }
        if(Same) {
          // Both original edges denote the identical new local address.
          // Materialize it after initialization, dominating the cloned CFG;
          // this changes no memory operation, edge, mask or pointer choice.
          IRBuilder<> Entry(F.getEntryBlock().getTerminator());
          Value *Canonical=Entry.CreateConstGEP1_64(Entry.getInt8Ty(),Base,CommonOffset,"value.same.address");
          Phi->replaceAllUsesWith(Canonical);Phi->eraseFromParent();Changed=true;continue;
        }
      }
      if(!(isa<CastInst>(Inst) || isa<BinaryOperator>(Inst) || isa<GetElementPtrInst>(Inst) || isa<ICmpInst>(Inst)))continue;
      if(Inst->getType()->isFPOrFPVectorTy() || llvm::any_of(Inst->operands(),[](const Use &U){return U->getType()->isFPOrFPVectorTy();}))continue;
      if(auto *C=ConstantFoldInstruction(Inst,DL)) {
        Inst->replaceAllUsesWith(C);Inst->eraseFromParent();Changed=true;
      }
    }
  } while(Changed);
}

// Check the new local core bank too. No original bank address can enter an
// opaque core call. Only ordinary bounded accesses and constant memcpy/memset
// may use these NEW local pointers; the original alloca checker is stricter.
void checkCoreBank(AllocaInst *A,unsigned Count,const DataLayout &DL) {
  uint64_t Bytes=uint64_t(Count)*4;
  SmallVector<Value *,32> Pending{A};SmallPtrSet<Value *,32> Seen;
  while(!Pending.empty()) {
    Value *P=Pending.pop_back_val();if(!Seen.insert(P).second)continue;
    for(User *U:P->users()) {
      if(auto *G=dyn_cast<GetElementPtrInst>(U)) {bounds(G,A,Bytes,0,DL);Pending.push_back(G);continue;}
      if(auto *L=dyn_cast<LoadInst>(U)) {
        if(L->getPointerOperand()!=P || L->isVolatile() || L->isAtomic())reject("core local load not ordinary",L);
        bounds(P,A,Bytes,sizeOf(L->getType(),DL),DL);continue;
      }
      if(auto *S=dyn_cast<StoreInst>(U)) {
        if(S->getPointerOperand()!=P || S->isVolatile() || S->isAtomic())reject("core bank address escapes",S);
        bounds(P,A,Bytes,sizeOf(S->getValueOperand()->getType(),DL),DL);continue;
      }
      if(auto *MT=dyn_cast<MemTransferInst>(U)) {
        auto *N=dyn_cast<ConstantInt>(MT->getLength());
        if(!N || MT->isVolatile() || (MT->getRawDest()!=P && MT->getRawSource()!=P))reject("unknown core local memcpy",MT);
        bounds(P,A,Bytes,N->getZExtValue(),DL);continue;
      }
      if(auto *MS=dyn_cast<MemSetInst>(U)) {
        auto *N=dyn_cast<ConstantInt>(MS->getLength());
        if(!N || MS->isVolatile() || MS->getRawDest()!=P)reject("unknown core local memset",MS);
        bounds(P,A,Bytes,N->getZExtValue(),DL);continue;
      }
      reject("new value-core local bank still escapes",dyn_cast<Instruction>(U));
    }
  }
}

using Inventory=std::map<std::string,unsigned>;
Inventory fpInventory(const Module &M) {
  Inventory Result;
  for(const Function &F:M)for(const Instruction &I:instructions(F)) {
    if(auto *Op=dyn_cast<FPMathOperator>(&I))if(Op->getFastMathFlags().any())reject("fast math is outside exact input",&I);
    if(auto *II=dyn_cast<ConstrainedFPIntrinsic>(&I))++Result[II->getCalledFunction()->getName().str()];
    else if(I.getOpcode()==Instruction::FAdd || I.getOpcode()==Instruction::FSub ||
            I.getOpcode()==Instruction::FMul || I.getOpcode()==Instruction::FDiv ||
            I.getOpcode()==Instruction::FCmp || I.getOpcode()==Instruction::FPTrunc ||
            I.getOpcode()==Instruction::FPExt)reject("ordinary FP operation outside strict copy",&I);
  }
  return Result;
}

// All external incoming values are integer vectors, never addresses. Original
// bank/tag writes now address core-local arrays; all old nonintrinsic math
// callees are the exact source-bound internal slow/value helpers. FP effects
// remain inaccessible RW; readonly constants justify conservative Other Ref.
MemoryEffects coreEffects() {
  return MemoryEffects::otherMemOnly(ModRefInfo::Ref) |
         MemoryEffects::inaccessibleMemOnly(ModRefInfo::ModRef);
}
FunctionType *coreType(LLVMContext &C,const Recipe &R) {
  Type *V=FixedVectorType::get(Type::getInt32Ty(C),8);
  SmallVector<Type *,13> Params(10,V);if(R.Tagged){Params.push_back(V);Params.push_back(V);}Params.push_back(Type::getInt32Ty(C));
  SmallVector<Type *,6> Results(R.Tagged?6:5,V);
  return FunctionType::get(StructType::get(C,Results),Params,false);
}
void configureCore(Function &F) {
  F.setLinkage(GlobalValue::ExternalLinkage);F.setVisibility(GlobalValue::HiddenVisibility);
  F.setDSOLocal(true);F.setCallingConv(CallingConv::Fast);
  F.removeFnAttr(Attribute::AlwaysInline);F.addFnAttr(Attribute::NoInline);
  F.addFnAttr(Attribute::StrictFP);F.addFnAttr(Attribute::NoUnwind);
  F.setMemoryEffects(coreEffects());
}
Function *declareCore(Module &M,const Recipe &R) {
  std::string Name=(SharedPrefix+R.Suffix).str();
  if(M.getFunction(Name))reject("shared-core name collision");
  auto *F=Function::Create(coreType(M.getContext(),R),GlobalValue::ExternalLinkage,Name,M);
  configureCore(*F);return F;
}
void pointerABI(Function &F,const Recipe &R,bool Definition) {
  unsigned First=R.Tagged?2:1;
  if(F.isDeclaration()==Definition || F.arg_size()!=First+4 || !F.getReturnType()->isVoidTy() ||
     !F.doesNotThrow() || (Definition && !F.hasFnAttribute(Attribute::StrictFP)))
    reject("unexpected original pointer math declaration/definition");
  for(unsigned I=0;I<F.arg_size();++I) {
    Type *T=F.getArg(I)->getType();
    if(I<First ? (!T->isPointerTy() || cast<PointerType>(T)->getAddressSpace()!=0) : !T->isIntegerTy(32))
      reject("original pointer math ABI differs");
  }
}
void buildAdapter(Function &F,Function &Core,const Recipe &R) {
  if(!F.isDeclaration())reject("raw-mode adapter would overwrite a definition");
  LLVMContext &C=F.getContext();Type *V=FixedVectorType::get(Type::getInt32Ty(C),8);
  unsigned First=R.Tagged?2:1,Dest=R.Tagged?0:2,A=R.Tagged?1:0,BSlot=R.Tagged?2:1;
  F.setLinkage(GlobalValue::InternalLinkage);F.setVisibility(GlobalValue::DefaultVisibility);F.setDSOLocal(true);
  F.removeFnAttr(Attribute::NoInline);F.addFnAttr(Attribute::AlwaysInline);F.addFnAttr(Attribute::StrictFP);
  F.setMemoryEffects(MemoryEffects::argMemOnly()|coreEffects());
  IRBuilder<> B(BasicBlock::Create(C,"shared.adapter",&F));
  SmallVector<Value *,13> Args;
  for(unsigned S:{A,BSlot})for(unsigned I=0;I<5;++I)
    Args.push_back(B.CreateAlignedLoad(V,slot(B,F.getArg(0),S*5+I),Align(4)));
  if(R.Tagged)for(unsigned S:{1u,2u})Args.push_back(B.CreateAlignedLoad(V,slot(B,F.getArg(1),S),Align(4)));
  Args.push_back(F.getArg(First+3));
  auto *Call=B.CreateCall(&Core,Args,"shared.result");Call->setCallingConv(CallingConv::Fast);
  Call->addFnAttr(Attribute::StrictFP);Call->setDoesNotThrow();Call->setMemoryEffects(coreEffects());
  for(unsigned I=0;I<5;++I)B.CreateAlignedStore(B.CreateExtractValue(Call,I),slot(B,F.getArg(0),Dest*5+I),Align(4));
  if(R.Tagged)B.CreateAlignedStore(B.CreateExtractValue(Call,5),slot(B,F.getArg(1),0),Align(4));
  B.CreateRetVoid();
}

// Clone once in the math-only mode. There is no copy of this arithmetic CFG
// in any raw stage module, and external linkage fixes the cross-object ABI.
Function *cloneMath(Function &F,const Recipe &R,const DataLayout &DL) {
  Module &M=*F.getParent();LLVMContext &C=M.getContext();
  if(!F.use_empty())reject("standalone math entry has an unexpected caller");
  auto *Core=declareCore(M,R);
  Type *V=FixedVectorType::get(Type::getInt32Ty(C),8);
  auto *ReturnType=cast<StructType>(Core->getReturnType());
  BasicBlock *Entry=BasicBlock::Create(C,"shared.entry",Core);IRBuilder<> B(Entry);
  auto *Bank=B.CreateAlloca(ArrayType::get(Type::getInt32Ty(C),120),nullptr,"shared.local.bank");Bank->setAlignment(Align(32));
  AllocaInst *Tags=nullptr;
  if(R.Tagged){Tags=B.CreateAlloca(ArrayType::get(Type::getInt32Ty(C),24),nullptr,"shared.local.tags");Tags->setAlignment(Align(32));}
  unsigned Dest=R.Tagged?0:2,A=R.Tagged?1:0,BSlot=R.Tagged?2:1;
  for(unsigned I=0;I<5;++I) {
    B.CreateAlignedStore(Constant::getNullValue(V),slot(B,Bank,Dest*5+I),Align(32));
    B.CreateAlignedStore(Core->getArg(I),slot(B,Bank,A*5+I),Align(32));
    B.CreateAlignedStore(Core->getArg(5+I),slot(B,Bank,BSlot*5+I),Align(32));
  }
  if(R.Tagged) {
    B.CreateAlignedStore(Constant::getNullValue(V),slot(B,Tags,0),Align(32));
    B.CreateAlignedStore(Core->getArg(10),slot(B,Tags,1),Align(32));
    B.CreateAlignedStore(Core->getArg(11),slot(B,Tags,2),Align(32));
  }
  ValueToValueMapTy Map;unsigned First=R.Tagged?2:1;
  Map[F.getArg(0)]=Bank;if(R.Tagged)Map[F.getArg(1)]=Tags;
  Map[F.getArg(First)]=B.getInt32(Dest);Map[F.getArg(First+1)]=B.getInt32(A);Map[F.getArg(First+2)]=B.getInt32(BSlot);
  Map[F.getArg(First+3)]=Core->getArg(Core->arg_size()-1);
  SmallVector<ReturnInst *,4> Returns;
  CloneFunctionInto(Core,&F,Map,CloneFunctionChangeType::LocalChangesOnly,Returns);
  B.CreateBr(cast<BasicBlock>(Map[&F.getEntryBlock()]));
  Core->setAttributes(AttributeList::get(C,F.getAttributes().getFnAttrs(),AttributeSet{},{}));
  configureCore(*Core);
  if(F.hasPersonalityFn())Core->setPersonalityFn(F.getPersonalityFn());
  if(Returns.empty())reject("shared math core has no normal return");
  for(auto *Old:Returns) {
    IRBuilder<> Out(Old);Value *Result=PoisonValue::get(ReturnType);
    for(unsigned I=0;I<5;++I)Result=Out.CreateInsertValue(Result,Out.CreateAlignedLoad(V,slot(Out,Bank,Dest*5+I),Align(32)),I);
    if(R.Tagged)Result=Out.CreateInsertValue(Result,Out.CreateAlignedLoad(V,slot(Out,Tags,0),Align(32)),5);
    Out.CreateRet(Result);Old->eraseFromParent();
  }
  foldFixedIndices(*Core,DL);checkCoreBank(Bank,120,DL);if(Tags)checkCoreBank(Tags,24,DL);
  F.eraseFromParent();return Core;
}

void mathClosure(Module &M) {
  unsigned Slow=0;
  for(Function &F:M)if(!F.isDeclaration()) {
    if(!F.hasFnAttribute(Attribute::StrictFP) || !F.doesNotThrow())reject("math definition outside strict/nounwind source closure");
    auto ME=F.getMemoryEffects();
    if(isModSet(ME.getModRef(IRMemLocation::Other)) || isModSet(ME.getModRef(IRMemLocation::ErrnoMem)))
      reject("math source may write ordinary nonargument/errno memory");
    if(F.getName().contains("slow_add_triple") || F.getName().contains("slow_multiply_triple")) {
      if(!F.hasFnAttribute(Attribute::NoInline))reject("original slow function lost noinline");++Slow;
    }
    for(Instruction &I:instructions(F))if(auto *Call=dyn_cast<CallBase>(&I)) {
      Function *Callee=Call->getCalledFunction();
      if(!Callee || (!Callee->isIntrinsic() && Callee->isDeclaration()))reject("new opaque math callback",Call);
    }
  }
  if(Slow!=2)reject("standalone math slow-core inventory differs");
}

struct SharedValue : PassInfoMixin<SharedValue> {
  const Profile *P;
  explicit SharedValue(const Profile &P):P(&P){}
  PreservedAnalyses run(Module &M,ModuleAnalysisManager &) {
    requirePinned(M,*P);const DataLayout &DL=M.getDataLayout();
    Inventory Before=fpInventory(M);
    for(Function &F:M)if(F.hasFnAttribute("target-features")) {
      SmallVector<StringRef,32> Features;F.getFnAttribute("target-features").getValueAsString().split(Features,',');
      if(llvm::is_contained(Features,StringRef("-fma")))reject("FMA capability disable would also disable AVX512");
    }
    if(P->Mode==3) {
      mathClosure(M);
      for(const Recipe &R:Recipes) {
        Function *F=M.getFunction((Prefix+R.Suffix).str());
        if(!F)reject("missing standalone math body");pointerABI(*F,R,true);
      }
      for(const Recipe &R:Recipes)cloneMath(*M.getFunction((Prefix+R.Suffix).str()),R,DL);
      if(fpInventory(M)!=Before)reject("math clone/removal changed complete strict FP inventory");
      // The sole unused dense_subtract8 body is left for ordinary internalize
      // and DCE; nine shared names must be preserved by the producer.
      errs()<<P->Pass<<": shared_core_definitions=9 old_pointer_exports_removed=9 FP_inventory_unchanged=1\n";
    } else {
      unsigned Exports=0;
      for(Function &F:M)if(!F.isDeclaration() && !F.hasLocalLinkage()) {
        if(F.getName()!=P->Export && F.getName()!=(Twine(P->Export)+MaskedSuffix).str())reject("raw export set differs");
        if(F.getVisibility()!=GlobalValue::HiddenVisibility)reject("raw export visibility differs");++Exports;
      }
      if(Exports!=2)reject("raw must retain both ordinary and masked exports");
      SmallVector<Bank,8> Banks;SmallPtrSet<AllocaInst *,32> Unique;
      DenseMap<Function *,std::unique_ptr<DominatorTree>> Dominators;
      unsigned Calls=0,Tags=0,Masked=0,Adapters=0;
      for(const Recipe &R:Recipes) {
        Function *F=M.getFunction((Prefix+R.Suffix).str());unsigned Expected=R.Calls[P->Mode];
        if(!F){if(Expected)reject("missing referenced raw pointer declaration");continue;}
        pointerABI(*F,R,false);unsigned First=R.Tagged?2:1,Count=0;
        for(User *U:F->users()) {
          auto *Call=dyn_cast<CallInst>(U);
          if(!Call || Call->getCalledFunction()!=F || !Call->hasFnAttr(Attribute::StrictFP))reject("nondirect/unscoped pointer math use");
          std::array<unsigned,3> Indices=R.Tagged?std::array<unsigned,3>{0,1,2}:std::array<unsigned,3>{2,0,1};
          for(unsigned I=0;I<3;++I) {
            auto *K=dyn_cast<ConstantInt>(Call->getArgOperand(First+I));
            if(!K || K->getZExtValue()!=Indices[I])reject("raw common math tuple differs",Call);
          }
          auto *A=dyn_cast<AllocaInst>(Call->getArgOperand(0));auto *T=R.Tagged?dyn_cast<AllocaInst>(Call->getArgOperand(1)):nullptr;
          if(!A || (R.Tagged && (!T || T==A)))reject("raw private arguments are not separate exact allocas",Call);
          auto &DT=Dominators[Call->getFunction()];if(!DT)DT=std::make_unique<DominatorTree>(*Call->getFunction());
          if(!Unique.insert(A).second)reject("raw private bank shared between sites",Call);
          auto Data=inspectBank(A,120,Call,DL);checkInitialization(Data,Indices[0],*DT,DL);Masked+=Data.Masked.size();Banks.push_back(std::move(Data));
          if(T) {
            if(!Unique.insert(T).second)reject("raw private tag bank shared between sites",Call);
            auto Tag=inspectBank(T,24,Call,DL);checkInitialization(Tag,0,*DT,DL);Masked+=Tag.Masked.size();Banks.push_back(std::move(Tag));++Tags;
          }
          ++Count;
        }
        if(Count!=Expected)reject(Twine("raw exact call inventory differs for ")+R.Suffix);Calls+=Count;
      }
      if(Calls!=P->Calls || Tags!=P->TagBanks || Banks.size()!=Calls+Tags)reject("raw bank/call totals differ");
      for(Bank &B:Banks)normalize(B);
      for(const Recipe &R:Recipes) {
        Function *Core=declareCore(M,R);Function *F=M.getFunction((Prefix+R.Suffix).str());
        if(F && !F->use_empty()){buildAdapter(*F,*Core,R);++Adapters;}
      }
      if(fpInventory(M)!=Before)reject("raw adapters/normalization changed strict FP operations");
      errs()<<P->Pass<<": shared_core_declarations=9 adapters="<<Adapters<<" raw_calls="<<Calls<<" tags="<<Tags
            <<" initialized_local_masked_operations_normalized="<<Masked<<" FP_inventory_unchanged=1\n";
    }
    if(verifyModule(M,&errs()))reject("LLVM rejected shared value ABI output");
    return PreservedAnalyses::none();
  }
};
} // namespace
extern "C" LLVM_ATTRIBUTE_WEAK PassPluginLibraryInfo llvmGetPassPluginInfo() {
  return {LLVM_PLUGIN_API_VERSION,"SiriusSharedValue",LLVM_VERSION_STRING,
    [](PassBuilder &PB){PB.registerPipelineParsingCallback(
      [](StringRef Name,ModulePassManager &MPM,ArrayRef<PassBuilder::PipelineElement>){
        for(const Profile &P:Profiles)if(Name==P.Pass){MPM.addPass(SharedValue(P));return true;}return false;
      });},nullptr};
}
