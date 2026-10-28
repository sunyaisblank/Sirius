// Change only the private nine-core calling convention, never its arithmetic.
#include "llvm/ADT/SmallVector.h"
#include "llvm/ADT/StringExtras.h"
#include "llvm/Config/llvm-config.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/IR/Module.h"
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
#include <array>
#include <memory>

#if LLVM_VERSION_MAJOR != 23 || LLVM_VERSION_MINOR != 1 || LLVM_VERSION_PATCH != 2
#error "Requires the pinned LLVM 23.1.2 API."
#endif
using namespace llvm;
namespace {
constexpr StringLiteral Directory = "/home/astra/.project/Sirius/bin/linux-gcc/native-cpu-strict-ir-linking/";
constexpr StringLiteral Old = "sirius_shared_value_", New = "sirius_register_value_";
struct Profile { const char *Pass, *File, *SHA, *Export; uint64_t Bytes; unsigned Mode; };
constexpr std::array<Profile,4> Profiles{{
 {"sirius-register-dp", "shared-dp.O2.bc", "09f18c6909944b9f23ba87a0937da50d92e68f55d5036412efba35d1cb470cb3", "sirius_scope_dp_packet", 8644716, 0},
 {"sirius-register-transport", "shared-transport.O2.bc", "922560a7126c6aaa15aa37e2933ea0947c4dad4e96809f4f3dd50b8f395e9f74", "sirius_scope_transport_packet", 13019744, 1},
 {"sirius-register-endpoint", "shared-endpoint.cleanup.O2.bc", "d5949b32641a9e2d58b02a5b9630352cb970366f18cf870db3ba118d227e4758", "sirius_scope_endpoint_packet", 17593768, 2},
 {"sirius-register-math", "shared-math.O2.bc", "c4b9d7a91cdb1c49165bfb43991b243523d0e27687fc099db3f53283f1262d1d", "", 297664, 3}
}};
struct Core { const char *Suffix; bool Tagged; std::array<unsigned,3> Calls; };
constexpr std::array<Core,9> Cores{{
 {"dense_add8",false,{3,8,13}}, {"dense_multiply8",false,{1,7,13}},
 {"dense_divide8",false,{1,3,7}}, {"dense_sqrt8",false,{1,1,2}},
 {"tagged_add8",true,{585,947,1607}}, {"tagged_subtract8",true,{284,10,197}},
 {"tagged_multiply8",true,{661,1482,1963}}, {"tagged_divide8",true,{45,240,154}},
 {"tagged_sqrt8",true,{0,3,3}}
}};
[[noreturn]] void reject(const Twine &Why) {
 errs()<<"sirius-register-value: "<<Why<<'\n';
 report_fatal_error("exact private calling-convention conversion refused",false);
}
std::string digest(StringRef S) { SHA256 H;H.update(S);return toHex(H.final(),true); }
std::string fingerprint(const Module &M) {
 std::string S;raw_string_ostream O(S);M.print(O,nullptr);O.flush();return digest(S);
}
void pinned(const Module &M,const Profile &P) {
 auto Buffer=MemoryBuffer::getFile((Directory+P.File).str());
 if(!Buffer || (*Buffer)->getBufferSize()!=P.Bytes || digest((*Buffer)->getBuffer())!=P.SHA)
  reject("exact optimized shared input pin absent/changed");
 LLVMContext C;SMDiagnostic D;auto Reference=parseIR((*Buffer)->getMemBufferRef(),D,C);
 if(!Reference)reject("cannot independently parse pinned input");
 Reference->setModuleIdentifier(M.getModuleIdentifier());
 if(fingerprint(*Reference)!=fingerprint(M))reject("full canonical module fingerprint differs");
}
bool vectorWords(Type *T) {
 auto *V=dyn_cast<FixedVectorType>(T);
 return V && V->getNumElements()==8 && V->getElementType()->isIntegerTy(32);
}
void signature(Function &F,const Core &C,bool Definition) {
 unsigned N=C.Tagged?12:10;
 auto *R=dyn_cast<StructType>(F.getReturnType());
 if(F.isDeclaration()==Definition || F.arg_size()!=N+1 || F.isVarArg() || !R ||
    R->isPacked() || R->getNumElements()!=(C.Tagged?6u:5u) ||
    F.getCallingConv()!=CallingConv::Fast || F.getVisibility()!=GlobalValue::HiddenVisibility ||
    !F.hasExternalLinkage() || !F.isDSOLocal() || !F.hasFnAttribute(Attribute::NoInline) ||
    !F.hasFnAttribute(Attribute::StrictFP) || !F.doesNotThrow())reject("core interface differs");
 for(unsigned I=0;I<N;++I)if(!vectorWords(F.getArg(I)->getType()))reject("core operand type differs");
 if(!F.getArg(N)->getType()->isIntegerTy(32))reject("core mask type differs");
 for(Type *T:R->elements())if(!vectorWords(T))reject("core result field type differs");
}
struct Edit {
 Function *F; std::string Original, Replacement; SmallVector<CallBase *,8> Calls;
};
void change(SmallVectorImpl<Edit> &Edits,bool Forward) {
 for(auto &E:Edits) {
  E.F->setName(Forward?E.Replacement:E.Original);
  if(!E.Calls.empty() || StringRef(E.Original).starts_with(Old)) {
   auto CC=Forward?CallingConv::X86_RegCall:CallingConv::Fast;
   E.F->setCallingConv(CC);for(auto *C:E.Calls)C->setCallingConv(CC);
  }
 }
}
struct RegisterValue : PassInfoMixin<RegisterValue> {
 const Profile *P; explicit RegisterValue(const Profile &P):P(&P){}
 PreservedAnalyses run(Module &M,ModuleAnalysisManager &) {
  pinned(M,*P);
  SmallVector<Edit,12> Edits;unsigned Total=0,Defined=0;
  for(const Core &C:Cores) {
   std::string Name=(Old+C.Suffix).str();auto *F=M.getFunction(Name);
   unsigned Expected=P->Mode==3?0:C.Calls[P->Mode];
   if(!F) { if(P->Mode==3 || Expected)reject("missing core");continue; }
   signature(*F,C,P->Mode==3);
   Edit E{F,Name,(New+C.Suffix).str(),{}};
   if(M.getNamedValue(E.Replacement))reject("new name collision");
   for(User *U:F->users()) {
    auto *Call=dyn_cast<CallBase>(U);
    if(!Call || Call->getCalledFunction()!=F || Call->getFunctionType()!=F->getFunctionType() ||
       Call->getCallingConv()!=CallingConv::Fast || !Call->hasFnAttr(Attribute::StrictFP) ||
       !Call->doesNotThrow())reject("unexpected core use");
    E.Calls.push_back(Call);
   }
   if(E.Calls.size()!=Expected)reject("core call multiplicity differs");
   Total+=Expected;Defined+=!F->isDeclaration();Edits.push_back(std::move(E));
  }
  if(P->Mode==3 ? Defined!=9 : Defined!=0)reject("unexpected core body ownership");
  unsigned Exports=0;
  for(Function &F:M)if(!F.isDeclaration() && !F.hasLocalLinkage() && !F.getName().starts_with(Old)) {
   if(P->Mode==3 || !F.getName().starts_with(P->Export))reject("unexpected raw export");
   std::string Name=F.getName().str(),Replacement=Name;
   Replacement.replace(0,StringRef("sirius_scope_").size(),"sirius_register_");
   if(M.getNamedValue(Replacement))reject("raw export name collision");
   Edits.push_back({&F,Name,Replacement,{}});++Exports;
  }
  if(P->Mode==3 ? Exports!=0 : Exports!=2)reject("raw export count differs");
  // An exact normalized round trip proves the entire module differs only in
  // these function/call conventions and private symbol names. Attributes,
  // instruction bodies, operand order, CFG, masks and memory effects survive.
  std::string Before=fingerprint(M);change(Edits,true);
  if(verifyModule(M,&errs()))reject("converted module fails LLVM verification");
  std::string After=fingerprint(M);change(Edits,false);
  if(fingerprint(M)!=Before)reject("non-ABI module change detected");
  change(Edits,true);if(fingerprint(M)!=After)reject("conversion is not reproducible");
  errs()<<P->Pass<<": calls="<<Total<<" core-definitions="<<Defined
        <<" before="<<Before<<" after="<<After<<" normalized-roundtrip=exact\n";
  return PreservedAnalyses::none();
 }
};
}
extern "C" LLVM_ATTRIBUTE_WEAK PassPluginLibraryInfo llvmGetPassPluginInfo() {
 return {LLVM_PLUGIN_API_VERSION,"SiriusRegisterValueABI",LLVM_VERSION_STRING,
 [](PassBuilder &PB) { PB.registerPipelineParsingCallback(
 [](StringRef Name,ModulePassManager &MPM,ArrayRef<PassBuilder::PipelineElement>) {
  for(const auto &P:Profiles)if(Name==P.Pass){MPM.addPass(RegisterValue(P));return true;}
  return false;
 }); }};
}
