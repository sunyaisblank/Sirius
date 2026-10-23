#include "product17.h"
#include <llvm/IR/Verifier.h>
#include <llvm/IRReader/IRReader.h>
#include <llvm/Support/SourceMgr.h>
#include <llvm/Support/FileSystem.h>
#include <llvm/Support/raw_ostream.h>

std::string siriusCanonicalModule(const llvm::Module &);
bool siriusVerifyInverse(llvm::StringRef, llvm::StringRef);

int main(int argc, char **argv) {
   if (argc != 3) return 2;
   llvm::LLVMContext context;
   llvm::SMDiagnostic error;
   auto module = llvm::parseIRFile(argv[1], error, context);
   if (!module) { error.print(argv[0], llvm::errs()); return 3; }
   if (llvm::verifyModule(*module, &llvm::errs())) return 4;
   const auto original = siriusCanonicalModule(*module);
   llvm::Function *owner = module->getFunction("cs_variant");
   if (!owner || owner->isDeclaration()) return 5;
   for (llvm::Function &f : *module)
      if (!f.isDeclaration() && &f != owner) return 6;
   const auto result = sirius_product17::outline(*owner);
   if (llvm::verifyModule(*module, &llvm::errs())) return 7;
   std::error_code ec;
   llvm::raw_fd_ostream output(argv[2], ec, llvm::sys::fs::OF_None);
   if (ec) return 8;
   module->print(output, nullptr);
   output.flush();
   if (output.has_error()) return 9;
   if (!siriusVerifyInverse(original, argv[2])) return 10;
   llvm::outs() << "{\"windows\":" << result.windows
                << ",\"exact_clone_check\":" << (result.exact_clone_check ? "true" : "false")
                << ",\"original_and_output_verifier_passed\":true"
                << ",\"serialized_inverse_matches_original_whole_module\":true}\n";
}
