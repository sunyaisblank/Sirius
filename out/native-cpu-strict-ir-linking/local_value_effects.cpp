// Narrow only ordinary memory effects; all strict floating effects survive.
#include "llvm/ADT/SmallVector.h"
#include "llvm/ADT/StringExtras.h"
#include "llvm/Config/llvm-config.h"
#include "llvm/IR/Constants.h"
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
#include <array>
#include <memory>
using namespace llvm;
#if LLVM_VERSION_MAJOR != 23 || LLVM_VERSION_MINOR != 1 || LLVM_VERSION_PATCH != 2
#error "Requires the pinned LLVM23.1.2 API."
#endif
namespace {
constexpr StringLiteral Directory="/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-register-value/";
constexpr StringLiteral Old="sirius_register_value_",New="sirius_local_value_";
struct Profile { const char *Name,*SHA;uint64_t Bytes;unsigned Calls; };
constexpr std::array<Profile,4> Profiles{{
 {"math","6e98d3a0e81e092fbc6adb9404a378ae8ea8ddc1d9ee355ee8e959e6813628e0",297520,0},
 {"dp","53df413939ca1fc894726e34f39daa30c804ecba2b5fbf933d2700bce91e5567",8644344,1581},
 {"transport","76baefdefa21317c002eef960bb756e7a8c471f02e1f96e19649b7a8da8a03fb",13019408,2701},
 {"endpoint","1bf9f678fdf60b93eb555cbe0873bf454f0b0b0e3795b804069f77e720b6aa5d",17593448,3959}
}};
[[noreturn]] void reject(const Twine &S) {errs()<<"sirius-local-value: "<<S<<'\n';report_fatal_error("closed local-memory effect conversion refused",false);}
std::string digest(StringRef S) {SHA256 H;H.update(S);return toHex(H.final(),true);}
std::string fingerprint(const Module &M) {std::string S;raw_string_ostream O(S);M.print(O,nullptr);O.flush();return digest(S);}
std::unique_ptr<Module> read(LLVMContext &C,const Profile &P) {
 auto Buffer=MemoryBuffer::getFile((Directory+P.Name+".bc").str());
 if(!Buffer || (*Buffer)->getBufferSize()!=P.Bytes || digest((*Buffer)->getBuffer())!=P.SHA)reject("input bitcode pin differs");
 SMDiagnostic D;auto M=parseIR((*Buffer)->getMemBufferRef(),D,C);if(!M)reject("cannot independently parse input");return M;
}
bool numeric(Type *T) {
 if(T->isVoidTy() || T->isIntegerTy() || T->isFloatingPointTy())return true;
 if(auto *V=dyn_cast<VectorType>(T))return numeric(V->getElementType());
 if(auto *S=dyn_cast<StructType>(T)) {for(Type *E:S->elements())if(!numeric(E))return false;return true;}
 return false;
}
// The exact-source proof in private_finite_closure.json binds closure pointer
// reloads to locally captured Triple/x and bounds every fixed bank access.
// This structural check rejects new channels that would invalidate that proof.
void closedMath(Module &M) {
 if(!M.global_empty() || !M.alias_empty() || !M.ifunc_empty() || M.hasModuleInlineAsm())reject("ordinary global memory channel");
 unsigned Definitions=0,PointerLoads=0,PointerStores=0,DynamicGEP=0;
 for(Function &F:M)if(!F.isDeclaration()) {
  ++Definitions;
  if(!numeric(F.getReturnType()) || !F.doesNotThrow() || !F.hasFnAttribute(Attribute::StrictFP))reject("unexpected math return/effect");
  if(!F.hasLocalLinkage())for(Argument &A:F.args())if(!numeric(A.getType()))reject("external pointer argument");
  for(Instruction &I:instructions(F)) {
   if(isa<IntToPtrInst>(I) || isa<PtrToIntInst>(I) || isa<VAArgInst>(I) || isa<AtomicRMWInst>(I) || isa<AtomicCmpXchgInst>(I) || isa<FenceInst>(I))reject("foreign address/atomic operation");
   if(auto *FI=dyn_cast<FreezeInst>(&I))if(FI->getType()->isPointerTy())reject("arbitrary frozen pointer");
   if(auto *LI=dyn_cast<LoadInst>(&I)) {if(LI->isVolatile() || LI->isAtomic())reject("volatile/atomic load");PointerLoads+=LI->getType()->isPointerTy();}
   if(auto *SI=dyn_cast<StoreInst>(&I)) {
    if(SI->isVolatile() || SI->isAtomic())reject("volatile/atomic store");
    PointerStores+=SI->getValueOperand()->getType()->isPointerTy();
   }
   if(auto *G=dyn_cast<GetElementPtrInst>(&I)) {APInt Offset(64,0);DynamicGEP+=!G->accumulateConstantOffset(M.getDataLayout(),Offset);}
   for(Use &U:I.operands())if(U->getType()->isPointerTy()) {
    if(isa<UndefValue>(U.get()) || isa<PoisonValue>(U.get()))reject("undefined pointer source");
    if(isa<GlobalValue>(U.get())) {
     auto *Call=dyn_cast<CallBase>(&I);if(!Call || U.get()!=Call->getCalledOperand())reject("global address used as data");
    }
   }
   if(auto *Call=dyn_cast<CallBase>(&I)) {
    auto *C=Call->getCalledFunction();if(!C)reject("opaque callback");
    if(!C->isIntrinsic()) {if(C->isDeclaration() || !C->hasLocalLinkage())reject("external callback");continue;}
    StringRef N=C->getName();
    if(N.starts_with("llvm.experimental.constrained."))continue;
    switch(C->getIntrinsicID()) {
     case Intrinsic::lifetime_start:case Intrinsic::lifetime_end:case Intrinsic::memcpy:case Intrinsic::memset:
     case Intrinsic::experimental_noalias_scope_decl:case Intrinsic::fabs:case Intrinsic::usub_sat:case Intrinsic::is_fpclass:break;
     default:reject("new intrinsic memory/observer channel");
    }
   }
  }
 }
 if(Definitions!=12 || PointerLoads!=2 || PointerStores!=30 || DynamicGEP)reject("closed math pointer inventory differs");
}
struct FunctionEdit {Function *F;std::string Before,After;MemoryEffects Effects;bool Core;};
struct CallEdit {CallBase *C;MemoryEffects Effects;};
struct LocalEffects : PassInfoMixin<LocalEffects> {
 const Profile *P;explicit LocalEffects(const Profile &P):P(&P){}
 PreservedAnalyses run(Module &M,ModuleAnalysisManager &) {
  LLVMContext C;auto Reference=read(C,*P);Reference->setModuleIdentifier(M.getModuleIdentifier());
  std::string Before=fingerprint(M);if(fingerprint(*Reference)!=Before)reject("complete module fingerprint differs");
  LLVMContext MC;auto Math=read(MC,Profiles[0]);closedMath(*Math);
  SmallVector<FunctionEdit,12> Functions;SmallVector<CallEdit,2000> Calls;unsigned Cores=0;
  auto Effects=MemoryEffects::inaccessibleMemOnly(ModRefInfo::ModRef);
  for(Function &F:M)if(F.getName().starts_with(Old)) {
   if(F.getCallingConv()!=CallingConv::X86_RegCall || !F.hasFnAttribute(Attribute::StrictFP) || !F.doesNotThrow() || !numeric(F.getReturnType()))reject("core interface changed");
   for(Argument &A:F.args())if(!numeric(A.getType()))reject("pointer operand at value boundary");
   auto ME=F.getMemoryEffects();
   if(ME.getModRef(IRMemLocation::Other)!=ModRefInfo::Ref || ME.getModRef(IRMemLocation::InaccessibleMem)!=ModRefInfo::ModRef)reject("old core effects differ");
   std::string Name=F.getName().str(),Replacement=Name;Replacement.replace(0,Old.size(),New.str());
   if(M.getNamedValue(Replacement))reject("name collision");
   Functions.push_back({&F,Name,Replacement,ME,true});++Cores;
   for(User *U:F.users()) {
    auto *Call=dyn_cast<CallBase>(U);if(!Call || Call->getCalledFunction()!=&F || Call->getCallingConv()!=CallingConv::X86_RegCall || !Call->hasFnAttr(Attribute::StrictFP))reject("unexpected core use");
    Calls.push_back({Call,Call->getMemoryEffects()});
   }
  }
  if(Cores!=(StringRef(P->Name)=="dp"?8u:9u) || Calls.size()!=P->Calls)reject("core inventory differs");
  for(Function &F:M)if(!F.isDeclaration() && !F.hasLocalLinkage() && !F.getName().starts_with(Old)) {
   if(!F.getName().starts_with("sirius_register_") || !F.getName().contains("_packet"))reject("unexpected raw export");
   std::string Name=F.getName().str(),Replacement=Name;Replacement.replace(0,StringRef("sirius_register_").size(),"sirius_local_");
   if(M.getNamedValue(Replacement))reject("raw name collision");Functions.push_back({&F,Name,Replacement,F.getMemoryEffects(),false});
  }
  auto Change=[&](bool Forward) {
   for(auto &E:Functions) {E.F->setName(Forward?E.After:E.Before);if(E.Core)E.F->setMemoryEffects(Forward?Effects:E.Effects);}
   for(auto &E:Calls)E.C->setMemoryEffects(Forward?Effects:E.Effects);
  };
  Change(true);if(verifyModule(M,&errs()))reject("converted module fails verifier");std::string After=fingerprint(M);
  Change(false);if(fingerprint(M)!=Before)reject("non-effect/rename difference detected");Change(true);if(fingerprint(M)!=After)reject("conversion not reproducible");
  errs()<<"sirius-local-"<<P->Name<<": core-count="<<Cores<<" calls="<<Calls.size()<<" only-ordinary-OtherRead-removed; strict/inaccessibleRW-preserved normalized-roundtrip=exact\n";
  return PreservedAnalyses::none();
 }
};
}
extern "C" LLVM_ATTRIBUTE_WEAK PassPluginLibraryInfo llvmGetPassPluginInfo() {
 return {LLVM_PLUGIN_API_VERSION,"SiriusLocalValueEffects",LLVM_VERSION_STRING,[](PassBuilder &PB) {
 PB.registerPipelineParsingCallback([](StringRef Name,ModulePassManager &MPM,ArrayRef<PassBuilder::PipelineElement>) {
  for(const auto &P:Profiles)if(Name==(Twine("sirius-local-")+P.Name).str()){MPM.addPass(LocalEffects(P));return true;}return false;
 });}};
}
