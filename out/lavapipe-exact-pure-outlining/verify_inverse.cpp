// Independent whole-module inverse: decode the emitted helper's actual return
// values, expand every serialized call, then compare canonical original IR.
// It deliberately does not include the forward matcher's template or indices.
#include <llvm/IR/Constants.h>
#include <llvm/IR/Instructions.h>
#include <llvm/IR/IRBuilder.h>
#include <llvm/IR/Module.h>
#include <llvm/IR/Verifier.h>
#include <llvm/IRReader/IRReader.h>
#include <llvm/Support/SourceMgr.h>
#include <llvm/Support/raw_ostream.h>
#include <llvm/Transforms/Utils/Cloning.h>

std::string siriusCanonicalModule(const llvm::Module &source) {
   auto m = llvm::CloneModule(source);
   m->setModuleIdentifier("sirius_product17_roundtrip");
   for (llvm::Function &f : *m) {
      for (llvm::Argument &a : f.args()) a.setName("");
      for (llvm::BasicBlock &bb : f) {
         bb.setName("");
         for (llvm::Instruction &i : bb) i.setName("");
      }
   }
   std::string text;
   llvm::raw_string_ostream out(text);
   m->print(out, nullptr);
   return text;
}

bool siriusVerifyInverse(llvm::StringRef original, llvm::StringRef path) {
   using namespace llvm;
   LLVMContext context;
   SMDiagnostic error;
   auto m = parseIRFile(path, error, context);
   if (!m || verifyModule(*m, &errs())) return false;
   auto *helper = m->getFunction("__sirius_product17_v1");
   if (!helper)
      return siriusCanonicalModule(*m) == original;
   if (helper->size() != 1 || helper->arg_size() != 2 ||
       helper->getCallingConv() != CallingConv::PreserveAll ||
       !helper->getReturnType()->isStructTy()) return false;
   auto *return_type = cast<StructType>(helper->getReturnType());
   if (return_type->getNumElements() != 2 ||
       return_type->getElementType(0) != helper->getArg(0)->getType() ||
       return_type->getElementType(1) != helper->getArg(1)->getType()) return false;
   auto *terminator = dyn_cast<ReturnInst>(helper->front().getTerminator());
   if (!terminator) return false;
   Value *returned = terminator->getReturnValue();
   Value *fields[2] = {nullptr,nullptr};
   std::vector<Instruction *> body;
   std::vector<Instruction *> inserts;
   for (Instruction &i : helper->front()) {
      if (&i == terminator) break;
      if (isa<InsertValueInst>(i)) {
         inserts.push_back(&i);
      } else {
         if (!inserts.empty() || !isa<BinaryOperator>(i) ||
             !i.getType()->isFPOrFPVectorTy() ||
             i.getRawSubclassOptionalData() != 0 || i.hasMetadata()) return false;
         body.push_back(&i);
      }
   }
   if (inserts.size() != 2 || body.empty()) return false;
   for (unsigned n = 0; n < 2; ++n) {
      auto *insert = dyn_cast<InsertValueInst>(returned);
      if (!insert || insert->getNumIndices() != 1) return false;
      const unsigned index = *insert->idx_begin();
      if (index >= 2 || fields[index]) return false;
      fields[index] = insert->getInsertedValueOperand();
      returned = insert->getAggregateOperand();
   }
   if (!isa<PoisonValue>(returned) || !fields[0] || !fields[1]) return false;
   std::vector<CallInst *> calls;
   for (User *u : helper->users()) {
      auto *call = dyn_cast<CallInst>(u);
      if (!call || call->getCalledFunction() != helper || call->arg_size() != 2 ||
          call->getCallingConv() != CallingConv::PreserveAll) return false;
      calls.push_back(call);
   }
   if (calls.empty()) return false;
   for (CallInst *call : calls) {
      DenseMap<Value *,Value *> map;
      map[helper->getArg(0)] = call->getArgOperand(0);
      map[helper->getArg(1)] = call->getArgOperand(1);
      for (Instruction *i : body) {
         Instruction *copy = i->clone();
         for (unsigned j = 0; j < copy->getNumOperands(); ++j) {
            Value *v = copy->getOperand(j);
            if (Value *mapped = map.lookup(v)) copy->setOperand(j,mapped);
            else if (!isa<Constant>(v)) { copy->deleteValue(); return false; }
         }
         copy->insertBefore(call->getIterator());
         map[i] = copy;
      }
      Value *results[2] = {map.lookup(fields[0]), map.lookup(fields[1])};
      if (!results[0] || !results[1]) return false;
      std::vector<ExtractValueInst *> projections;
      for (User *u : call->users()) {
         auto *extract = dyn_cast<ExtractValueInst>(u);
         if (!extract || extract->getNumIndices() != 1 ||
             *extract->idx_begin() >= 2) return false;
         projections.push_back(extract);
      }
      for (ExtractValueInst *extract : projections) {
         extract->replaceAllUsesWith(results[*extract->idx_begin()]);
         extract->eraseFromParent();
      }
      call->eraseFromParent();
   }
   helper->eraseFromParent();
   if (verifyModule(*m, &errs())) return false;
   const auto restored = siriusCanonicalModule(*m);
   if (restored != original) {
      errs() << "serialized inverse differs from original whole-module IR\n";
      return false;
   }
   return true;
}
