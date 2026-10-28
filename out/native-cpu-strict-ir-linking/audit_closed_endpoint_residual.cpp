// Read-only audit of one retained Endpoint compiler artifact, LLVM23.1.2 only.
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/ADT/StringExtras.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/Config/llvm-config.h"
#include "llvm/IR/Dominators.h"
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
#include "llvm/Support/JSON.h"
#include "llvm/Support/MemoryBuffer.h"
#include "llvm/Support/ModRef.h"
#include "llvm/Support/SHA256.h"
#include "llvm/Support/SourceMgr.h"
#include "llvm/Support/raw_ostream.h"
#include <map>
#include <string>
#if LLVM_VERSION_MAJOR != 23 || LLVM_VERSION_MINOR != 1 || LLVM_VERSION_PATCH != 2
#error "Requires the pinned LLVM23.1.2 API."
#endif
using namespace llvm;
namespace {
constexpr StringLiteral BC="/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-canonical-closed-endpoint/endpoint.closed.bc";
constexpr StringLiteral LL="/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-canonical-closed-endpoint/endpoint.closed.ll";
constexpr StringLiteral BCSHA="2c794bcf6a0171171366e62c45ea0011e0480045abc6dff4ff58ac5bbf746075";
constexpr StringLiteral LLSHA="46bcf17beabf9b74d882bc70eca31b99f70669dea90172788b674b665acecbd8";
constexpr StringLiteral Prefix="sirius_local_value_";
[[noreturn]] void reject(const Twine &S,const Instruction *I=nullptr) {
  errs()<<"closed-endpoint-residual-audit: "<<S<<'\n';
  if(I) { I->print(errs());errs()<<'\n'; }
  report_fatal_error("Retained Endpoint artifact audit refused",false);
}
std::string digest(StringRef S) { SHA256 H;H.update(S);return toHex(H.final(),true); }
std::string fingerprint(const Module &M) {
  std::string S;raw_string_ostream O(S);M.print(O,nullptr);O.flush();return digest(S);
}
void reference(const Module &M,StringRef Path,uint64_t Bytes,StringRef SHA,const std::string &Expected) {
  auto Buffer=MemoryBuffer::getFile(Path);
  if(!Buffer||(*Buffer)->getBufferSize()!=Bytes||digest((*Buffer)->getBuffer())!=SHA) reject("complete retained input pin differs: "+Path);
  LLVMContext Context;SMDiagnostic Error;
  auto Parsed=parseIR((*Buffer)->getMemBufferRef(),Error,Context);
  if(!Parsed) reject("independent reference parse failed: "+Path);
  Parsed->setModuleIdentifier(M.getModuleIdentifier());
  if(fingerprint(*Parsed)!=Expected) reject("full canonical module differs from independent retained reference: "+Path);
}
bool words(Type *T) {
  auto *V=dyn_cast<FixedVectorType>(T);
  return V&&V->getNumElements()==8&&V->getElementType()->isIntegerTy(32);
}
bool vector32(Type *T) {
  auto *V=dyn_cast<FixedVectorType>(T);
  return V&&V->getNumElements()==8&&(V->getElementType()->isIntegerTy(32)||V->getElementType()->isFloatTy());
}
uint64_t bytes(Type *T,const DataLayout &DL) {
  auto S=DL.getTypeStoreSize(T);if(S.isScalable()) reject("scalable memory type");return S.getFixedValue();
}
using Counts=std::map<std::string,unsigned>;
json::Object object(const Counts &C) {json::Object O;for(auto &P:C)O[P.first]=int64_t(P.second);return O;}
struct Inventory {
  Counts Constrained,Math;unsigned FNeg=0,FAbs=0,Loads=0,Stores=0,MaskedLoads=0,MaskedStores=0;
  bool operator==(const Inventory &Other) const {
    return Constrained==Other.Constrained&&Math==Other.Math&&FNeg==Other.FNeg&&FAbs==Other.FAbs&&Loads==Other.Loads&&
      Stores==Other.Stores&&MaskedLoads==Other.MaskedLoads&&MaskedStores==Other.MaskedStores;
  }
};
std::string metadataString(Value *V) {
  auto *M=dyn_cast<MetadataAsValue>(V);auto *S=M?dyn_cast<MDString>(M->getMetadata()):nullptr;
  return S?S->getString().str():std::string();
}
Inventory inventory(Module &M) {
  Inventory Out;
  for(Function &F:M) for(Instruction &I:instructions(F)) {
    if(auto *FP=dyn_cast<FPMathOperator>(&I)) if(FP->getFastMathFlags().any())reject("fast math",&I);
    if(auto *CF=dyn_cast<ConstrainedFPIntrinsic>(&I)) {
      StringRef Name=CF->getCalledFunction()->getName();++Out.Constrained[Name.str()];
      if(metadataString(CF->getArgOperand(CF->arg_size()-1))!="fpexcept.strict")reject("non-strict FP exceptions",CF);
      if(!Name.contains(".fcmp")&&metadataString(CF->getArgOperand(CF->arg_size()-2))!="round.dynamic")reject("non-dynamic FP rounding",CF);
    } else if(I.getOpcode()==Instruction::FNeg)++Out.FNeg;
    else if(I.getOpcode()==Instruction::FAdd||I.getOpcode()==Instruction::FSub||I.getOpcode()==Instruction::FMul||
      I.getOpcode()==Instruction::FDiv||I.getOpcode()==Instruction::FRem||I.getOpcode()==Instruction::FCmp||
      I.getOpcode()==Instruction::FPTrunc||I.getOpcode()==Instruction::FPExt||I.getOpcode()==Instruction::SIToFP||
      I.getOpcode()==Instruction::UIToFP||I.getOpcode()==Instruction::FPToSI||I.getOpcode()==Instruction::FPToUI)
        reject("ordinary FP outside representation fneg",&I);
    if(isa<LoadInst>(I))++Out.Loads;if(isa<StoreInst>(I))++Out.Stores;
    if(auto *II=dyn_cast<IntrinsicInst>(&I)) {
      if(II->getIntrinsicID()==Intrinsic::masked_store)++Out.MaskedStores;
      if(II->getIntrinsicID()==Intrinsic::masked_load)++Out.MaskedLoads;
    }
    if(auto *C=dyn_cast<CallBase>(&I)) if(auto *Target=C->getCalledFunction()) if(Target->getName().starts_with("llvm.fabs."))++Out.FAbs;
    if(auto *C=dyn_cast<CallBase>(&I)) if(auto *Target=C->getCalledFunction()) if(Target->getName().starts_with(Prefix)) {
      if(C->getCallingConv()!=CallingConv::ID(92)||C->getFunctionType()!=Target->getFunctionType())reject("core call ABI drift",C);
      ++Out.Math[Target->getName().str()];
    }
  }
  return Out;
}
void mathDeclarations(Module &M) {
  const char *Suffixes[]={"dense_add8","dense_multiply8","dense_divide8","dense_sqrt8","tagged_add8", "tagged_subtract8","tagged_multiply8","tagged_divide8","tagged_sqrt8"};
  unsigned Found=0;
  for(Function &F:M) if(F.getName().starts_with(Prefix)) {
    bool Known=false;for(const char *S:Suffixes)if(F.getName()==(Prefix+S).str())Known=true;
    if(!Known||!F.isDeclaration()||F.getCallingConv()!=CallingConv::ID(92)||!F.hasFnAttribute(Attribute::StrictFP)||
       F.getMemoryEffects()!=MemoryEffects::inaccessibleMemOnly(ModRefInfo::ModRef)||F.isVarArg())reject("math9 declaration/effects drift");
    bool Tagged=F.getName().contains("_tagged_");unsigned N=Tagged?12:10;
    if(F.arg_size()!=N+1)reject("math9 parameter count differs");
    for(unsigned K=0;K<N;++K)if(!words(F.getArg(K)->getType()))reject("math9 vector parameter differs");
    if(!F.getArg(N)->getType()->isIntegerTy(32))reject("math9 mask differs");
    auto *T=dyn_cast<StructType>(F.getReturnType());
    if(!T||T->isOpaque()||T->getNumElements()!=(Tagged?6:5))reject("math9 return differs");
    for(Type *E:T->elements())if(!words(E))reject("math9 result vector differs");
    ++Found;
  }
  if(Found!=9)reject("nine fixed external math declarations required");
}
struct LocalCounts {unsigned Loads=0,Stores=0,MaskedStores=0,MaskedLoads=0,Memsets=0,GEPs=0,Lifetimes=0;};
LocalCounts local(AllocaInst *A,uint64_t Size,Function &F,const DataLayout &DL) {
  LocalCounts Counts;SmallVector<Value *,8> Pending{A};SmallPtrSet<Value *,16> Seen;
  SmallVector<Instruction *,32> Memory;SmallVector<Instruction *,4> Initializers;
  while(!Pending.empty()) {
    Value *P=Pending.pop_back_val();if(!Seen.insert(P).second)continue;
    int64_t O=0;if(GetPointerBaseWithConstantOffset(P,O,DL)!=A||O<0||uint64_t(O)>Size||O%32)reject("unbounded residual address",dyn_cast<Instruction>(P));
    for(User *U:P->users()) {
      if(auto *G=dyn_cast<GetElementPtrInst>(U)) {
        if(G->getPointerOperand()!=P||!G->hasAllConstantIndices())reject("dynamic/captured residual GEP",G);
        Pending.push_back(G);++Counts.GEPs;continue;
      }
      if(auto *L=dyn_cast<LoadInst>(U)) {
        if(L->getPointerOperand()!=P||L->isAtomic()||L->isVolatile()||L->getAlign()>Align(32)||!vector32(L->getType())||uint64_t(O)+32>Size)reject("residual load differs",L);
        Memory.push_back(L);++Counts.Loads;continue;
      }
      if(auto *S=dyn_cast<StoreInst>(U)) {
        if(S->getPointerOperand()!=P||S->isAtomic()||S->isVolatile()||S->getAlign()>Align(32)||!vector32(S->getValueOperand()->getType())||uint64_t(O)+32>Size)reject("residual store/capture differs",S);
        Memory.push_back(S);++Counts.Stores;if(!O&&Size==32)Initializers.push_back(S);continue;
      }
      auto *II=dyn_cast<IntrinsicInst>(U);if(!II)reject("residual pointer escapes",dyn_cast<Instruction>(U));
      if((II->getIntrinsicID()==Intrinsic::lifetime_start||II->getIntrinsicID()==Intrinsic::lifetime_end)&&
          II->arg_size()==1&&II->getArgOperand(0)==P&&!O) {++Counts.Lifetimes;continue;}
      if(auto *Set=dyn_cast<MemSetInst>(II)) {
        auto *N=dyn_cast<ConstantInt>(Set->getLength());auto *Z=dyn_cast<ConstantInt>(Set->getValue());
        if(Size!=160||O||Set->getRawDest()!=P||!N||N->getZExtValue()!=160||!Z||!Z->isZero()||Set->isVolatile())reject("residual memset differs",Set);
        Memory.push_back(Set);Initializers.push_back(Set);++Counts.Memsets;continue;
      }
      bool Load=II->getIntrinsicID()==Intrinsic::masked_load;
      if((!Load&&II->getIntrinsicID()!=Intrinsic::masked_store)||II->arg_size()!=3||II->getArgOperand(Load?0:1)!=P||
          !vector32(Load?II->getType():II->getArgOperand(0)->getType())||uint64_t(O)+32>Size)reject("residual callback/intrinsic differs",II);
      auto *Mask=dyn_cast<FixedVectorType>(II->getArgOperand(Load?1:2)->getType());
      if(!Mask||Mask->getNumElements()!=8||!Mask->getElementType()->isIntegerTy(1))reject("residual mask ABI differs",II);
      if(II->getParamAlign(Load?0:1).valueOrOne()>Align(32))reject("unproved residual alignment",II);
      Memory.push_back(II);if(Load)++Counts.MaskedLoads;else++Counts.MaskedStores;
    }
  }
  DominatorTree DT(F);unsigned FullInitializers=0;
  for(Instruction *Init:Initializers) {
    bool Dominates=true;for(Instruction *I:Memory)if(I!=Init&&!DT.dominates(Init,I))Dominates=false;
    if(Dominates)++FullInitializers;
  }
  if(FullInitializers!=1||Counts.Lifetimes!=2||Counts.MaskedLoads)reject("residual full initialization/lifetime proof differs",A);
  if(Size==160&&(Counts.Memsets!=1||Counts.GEPs!=4))reject("root must have sole160B zero and four constant fields",A);
  if(Size==32&&(Counts.Memsets||Counts.GEPs))reject("metadata vector has unexpected descendants",A);
  return Counts;
}
struct Audit:PassInfoMixin<Audit> {
  PreservedAnalyses run(Module &M,ModuleAnalysisManager &) {
    std::string Before=fingerprint(M);reference(M,BC,19497068,BCSHA,Before);reference(M,LL,119473814,LLSHA,Before);
    Function *Entry=M.getFunction("sirius_closed_endpoint_packet");unsigned Definitions=0;
    for(Function &F:M)if(!F.isDeclaration()){++Definitions;if(&F!=Entry)reject("unexpected defined raw/helper");}
    if(Definitions!=1||!Entry||Entry->arg_size()!=6||Entry->isVarArg()||!Entry->getReturnType()->isVoidTy()||
        Entry->getCallingConv()!=CallingConv::C||!Entry->hasFnAttribute(Attribute::StrictFP)||
        Entry->getUWTableKind()!=UWTableKind::Async||Entry->getVisibility()!=GlobalValue::HiddenVisibility)reject("sole ordinary C ABI/strict/Async entry differs");
    for(unsigned K=0;K<6;++K)if((K==0||K==2)?!Entry->getArg(K)->getType()->isPointerTy():!Entry->getArg(K)->getType()->isIntegerTy(32))reject("six-argument C ABI differs");
    mathDeclarations(M);Inventory Inv=inventory(M);const DataLayout &DL=M.getDataLayout();
    std::map<std::string,unsigned> Wanted={{"programBase_3.i",32},{"count_2.i",32},{"layers_2.i",32},{"prefix_5.i",32},{"component_4.i",32},{"root_1.i",160}};
    unsigned Total=0,Staging=0,Residual=0;LocalCounts Sum;json::Array Locals;
    for(Instruction &I:instructions(Entry))if(auto *A=dyn_cast<AllocaInst>(&I)) {
      ++Total;auto *N=dyn_cast<ConstantInt>(A->getArraySize());if(!N||!N->isOne()||A->getAddressSpace()!=0)reject("allocation count/address space differs",A);
      if(A->getName()=="staged_output") {
        auto *T=dyn_cast<ArrayType>(A->getAllocatedType());if(!T||T->getNumElements()!=30632||!T->getElementType()->isIntegerTy(32)||A->getAlign()!=Align(64))reject("transactional staging differs",A);
        ++Staging;continue;
      }
      auto It=Wanted.find(A->getName().str());if(It==Wanted.end()||!It->second)reject("unexpected private residual/Context",A);
      unsigned Size=It->second;It->second=0;
      if(bytes(A->getAllocatedType(),DL)!=Size||A->getAlign()!=Align(32))reject("residual size/alignment differs",A);
      if(Size==32) {if(!words(A->getAllocatedType()))reject("metadata vector type differs",A);}
      else {auto *S=dyn_cast<StructType>(A->getAllocatedType());if(!S||S->isOpaque()||S->isPacked()||S->getNumElements()!=5)reject("root struct type differs",A);
        for(unsigned K=0;K<5;++K)if(!vector32(S->getElementType(K))||DL.getStructLayout(S)->getElementOffset(K)!=K*32)reject("root struct fields differ",A);
        for(unsigned K=0;K<4;++K)if(!cast<FixedVectorType>(S->getElementType(K))->getElementType()->isFloatTy())reject("root FP representation differs",A);
        if(!words(S->getElementType(4)))reject("root validity word differs",A);}
      LocalCounts C=local(A,Size,*Entry,DL);Residual+=Size;
      Sum.Loads+=C.Loads;Sum.Stores+=C.Stores;Sum.MaskedStores+=C.MaskedStores;Sum.MaskedLoads+=C.MaskedLoads;Sum.Memsets+=C.Memsets;Sum.GEPs+=C.GEPs;Sum.Lifetimes+=C.Lifetimes;
      Locals.push_back(json::Object{{"name",A->getName().str()},{"bytes",int64_t(Size)},{"alignment",32},{"full_initialization_dominates_every_memory_use",true},{"no_capture_or_dynamic_address",true}});
    }
    if(Total!=7||Staging!=1||Residual!=320||Sum.Loads!=17||Sum.Stores!=11||Sum.MaskedStores!=16||Sum.MaskedLoads||Sum.Memsets!=1||Sum.GEPs!=4||Sum.Lifetimes!=12)reject("exact residual inventory differs");
    for(auto &P:Wanted)if(P.second)reject("required residual absent");
    if(Inv.Loads!=579||Inv.Stores!=11||Inv.MaskedStores!=572||Inv.MaskedLoads!=29)reject("whole-function memory inventory differs");
    Counts ExpectedMath={{"sirius_local_value_dense_multiply8",12},{"sirius_local_value_dense_add8",11},
      {"sirius_local_value_dense_divide8",6},{"sirius_local_value_dense_sqrt8",1},
      {"sirius_local_value_tagged_multiply8",1963},{"sirius_local_value_tagged_divide8",154},
      {"sirius_local_value_tagged_add8",1607},{"sirius_local_value_tagged_sqrt8",3},{"sirius_local_value_tagged_subtract8",197}};
    Counts ExpectedFP={{"llvm.experimental.constrained.fadd.v8f32",96},{"llvm.experimental.constrained.fsub.v8f32",13},
      {"llvm.experimental.constrained.fcmp.v8f32",451},{"llvm.experimental.constrained.fcmps.v8f32",28},
      {"llvm.experimental.constrained.fmul.v8f32",1}};
    if(Inv.Math!=ExpectedMath||Inv.Constrained!=ExpectedFP||Inv.FNeg!=16||Inv.FAbs!=312)reject("exact retained strict FP/math inventory differs");
    if(verifyModule(M,&errs()))reject("invalid retained module");
    if(fingerprint(M)!=Before||!(inventory(M)==Inv))reject("read-only audit changed module or FP/call inventory");
    json::Object Out{{"passed",true},{"read_only",true},{"input_sha256",BCSHA.str()},{"text_sha256",LLSHA.str()},
      {"canonical_module_fingerprint",Before},{"full_module_fingerprint_unchanged",true},{"FP_call_inventory_unchanged",true},
      {"definitions",1},{"owned_context_allocas",0},{"strict_dynamic_exceptions",true},{"ordinary_FP_other_than_fneg",0},
      {"ordinary_entry_async_uwtable",true},{"math9_fixed_cc92_inaccessibleRW",true},{"staging_bytes",122528},
      {"residual_bytes",320},{"residual_loads",17},{"residual_stores",11},{"residual_masked_stores",16},{"residual_masked_loads",0},
      {"residual_memory_operations",44},{"residual_root_zero_memsets",1},{"residuals",std::move(Locals)},
      {"constrained_calls",object(Inv.Constrained)},{"value_calls",object(Inv.Math)},{"ordinary_fneg",int64_t(Inv.FNeg)},{"fabs_calls",int64_t(Inv.FAbs)},{"value_call_total",3954},
      {"whole_loads",579},{"whole_stores",11},{"whole_masked_stores",572},{"whole_masked_loads",29},
      {"numerical_qualification_transferred",false}};
    outs()<<json::Value(std::move(Out))<<'\n';return PreservedAnalyses::all();
  }
};
}
extern "C" LLVM_ATTRIBUTE_WEAK PassPluginLibraryInfo llvmGetPassPluginInfo() {
  return {LLVM_PLUGIN_API_VERSION,"SiriusClosedEndpointResidualAudit",LLVM_VERSION_STRING,
    [](PassBuilder &PB){PB.registerPipelineParsingCallback([](StringRef Name,ModulePassManager &MPM,ArrayRef<PassBuilder::PipelineElement>){
      if(Name!="sirius-audit-closed-endpoint-residual")return false;MPM.addPass(Audit());return true;});},nullptr};
}
