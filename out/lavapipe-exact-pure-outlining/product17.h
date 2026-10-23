// Private LLVM 21 code-screen candidate. No floating operation is reassociated.
#pragma once

#include <array>
#include <iterator>
#include <vector>

#include <llvm/ADT/DenseMap.h>
#include <llvm/IR/Constants.h>
#include <llvm/IR/Function.h>
#include <llvm/IR/IRBuilder.h>
#include <llvm/IR/Instructions.h>
#include <llvm/IR/Module.h>
#include <llvm/Support/ModRef.h>

namespace sirius_product17 {
using namespace llvm;

struct Window {
   std::array<Instruction *, 17> body;
};

// Nonnegative operands name earlier definitions; negative operands are A,B,C.
// C is the exact binary32 split constant 4097, never a mask-derived value.
inline constexpr int operands[17][2] = {
   {-1,-3}, {0,-1}, {0,1}, {-1,2}, {-2,-3}, {4,-2}, {4,5}, {-2,6},
   {-1,-2}, {6,2}, {8,9}, {6,3}, {10,11}, {7,2}, {12,13}, {3,7}, {15,14}
};
inline constexpr unsigned opcodes[17] = {
   Instruction::FMul, Instruction::FSub, Instruction::FSub, Instruction::FSub,
   Instruction::FMul, Instruction::FSub, Instruction::FSub, Instruction::FSub,
   Instruction::FMul, Instruction::FMul, Instruction::FSub, Instruction::FMul,
   Instruction::FSub, Instruction::FMul, Instruction::FSub, Instruction::FMul,
   Instruction::FSub
};

inline bool match(Instruction *start, Window &w) {
   auto *vt = dyn_cast<FixedVectorType>(start->getType());
   if (!vt || vt->getNumElements() != 8 || !vt->getElementType()->isFloatTy())
      return false;
   Instruction *p = start;
   for (unsigned i = 0; i < 17; ++i) {
      if (!p || !isa<BinaryOperator>(p) || p->getOpcode() != opcodes[i] ||
          p->getType() != vt || p->getRawSubclassOptionalData() != 0 || p->hasMetadata())
         return false;
      w.body[i] = p;
      p = p->getNextNode();
   }
   Value *a = w.body[0]->getOperand(0), *b = w.body[4]->getOperand(0);
   // Alias-specialized/constant cases retain the original code. The first
   // representative therefore supplies an unambiguous two-argument map.
   if (a == b || isa<Constant>(a) || isa<Constant>(b))
      return false;
   // Both arguments must already dominate the replacement at the FIRST
   // instruction. In particular B may not be a private split of A above it.
   for (Instruction *inside : w.body)
      if (a == inside || b == inside)
         return false;
   auto *c = dyn_cast<Constant>(w.body[0]->getOperand(1));
   auto *split = c ? dyn_cast_or_null<ConstantFP>(c->getSplatValue()) : nullptr;
   if (!split || split->getValueAPF().bitcastToAPInt().getZExtValue() != 0x45800800u)
      return false;
   for (unsigned i = 0; i < 17; ++i) {
      for (unsigned j = 0; j < 2; ++j) {
         const int ref = operands[i][j];
         Value *expected;
         if (ref == -1) expected = a;
         else if (ref == -2) expected = b;
         else if (ref == -3) expected = c;
         else expected = w.body[ref];
         if (w.body[i]->getOperand(j) != expected)
            return false;
      }
      if (i == 8 || i == 16)
         continue;
      for (const User *user : w.body[i]->users()) {
         bool internal = false;
         for (Instruction *inside : w.body)
            internal |= user == inside;
         if (!internal)
            return false;
      }
   }
   return true;
}

struct Result {
   unsigned windows = 0;
   bool exact_clone_check = false;
};

// The provider caller must separately admit the exact original module, owner,
// target and cache namespace. This function changes only exact matched windows.
inline Result outline(Function &owner) {
   Module &m = *owner.getParent();
   if (m.getNamedValue("__sirius_product17_v1"))
      return {};
   std::vector<Window> windows;
   for (BasicBlock &bb : owner) {
      for (auto it = bb.begin(); it != bb.end();) {
         Window w;
         if (!match(&*it, w)) {
            ++it;
            continue;
         }
         windows.push_back(w);
         it = std::next(w.body.back()->getIterator());
      }
   }
   if (windows.empty())
      return {};

   auto &first = windows.front();
   Type *vt = first.body[0]->getType();
   auto *rt = StructType::get(m.getContext(), {vt,vt});
   auto *ft = FunctionType::get(rt, {vt,vt}, false);
   auto *helper = Function::Create(ft, GlobalValue::InternalLinkage,
                                   "__sirius_product17_v1", m);
   helper->setCallingConv(CallingConv::PreserveAll);
   helper->addFnAttr(Attribute::NoInline);
   helper->addFnAttr(Attribute::NoUnwind);
   helper->addFnAttr(Attribute::NoSync);
   helper->addFnAttr(Attribute::NoFree);
   helper->addFnAttr(Attribute::NoRecurse);
   helper->addFnAttr(Attribute::WillReturn);
   helper->setMemoryEffects(MemoryEffects::none());
   for (const char *name : {"target-cpu", "target-features", "denormal-fp-math",
                            "denormal-fp-math-f32"}) {
      if (owner.hasFnAttribute(name))
         helper->addFnAttr(owner.getFnAttribute(name));
   }
   if (owner.hasFnAttribute(Attribute::StrictFP))
      helper->addFnAttr(Attribute::StrictFP);

   DenseMap<Value *, Value *> mapping;
   mapping[first.body[0]->getOperand(0)] = helper->getArg(0);
   mapping[first.body[4]->getOperand(0)] = helper->getArg(1);
   auto *entry = BasicBlock::Create(m.getContext(), "entry", helper);
   std::array<Instruction *,17> cloned;
   for (unsigned i = 0; i < 17; ++i) {
      Instruction *copy = first.body[i]->clone();
      for (unsigned j = 0; j < copy->getNumOperands(); ++j)
         if (Value *mapped = mapping.lookup(copy->getOperand(j)))
            copy->setOperand(j, mapped);
      copy->insertInto(entry, entry->end());
      cloned[i] = copy;
      mapping[first.body[i]] = copy;
   }
   IRBuilder<> ret(entry);
   Value *pair = ret.CreateInsertValue(PoisonValue::get(rt), cloned[8], 0);
   pair = ret.CreateInsertValue(pair, cloned[16], 1);
   ret.CreateRet(pair);

   // Re-expand the actual helper against EVERY original region before the
   // first caller mutation. isIdenticalTo also compares optional flags, unlike
   // stock similarity's isSameOperationAs. No extra metadata is accepted.
   for (const Window &w : windows) {
      DenseMap<Value *, Value *> back;
      back[helper->getArg(0)] = w.body[0]->getOperand(0);
      back[helper->getArg(1)] = w.body[4]->getOperand(0);
      for (unsigned i = 0; i < 17; ++i) {
         Instruction *copy = cloned[i]->clone();
         for (unsigned j = 0; j < copy->getNumOperands(); ++j)
            if (Value *mapped = back.lookup(copy->getOperand(j)))
               copy->setOperand(j, mapped);
         const bool exact = copy->isIdenticalTo(w.body[i]) &&
                            copy->hasSameSubclassOptionalData(w.body[i]);
         copy->deleteValue();
         if (!exact) {
            helper->eraseFromParent();
            return {};
         }
         back[cloned[i]] = w.body[i];
      }
   }

   for (const Window &w : windows) {
      // Read current operands: a preceding replacement may have rewritten
      // an input to this window. The same replacement also rewrites PHI uses.
      Value *a = w.body[0]->getOperand(0), *b = w.body[4]->getOperand(0);
      IRBuilder<> call(w.body[0]);
      auto *invoke = call.CreateCall(helper, {a,b}, "product17");
      invoke->setCallingConv(helper->getCallingConv());
      w.body[8]->replaceAllUsesWith(call.CreateExtractValue(invoke, 0));
      w.body[16]->replaceAllUsesWith(call.CreateExtractValue(invoke, 1));
      for (auto it = w.body.rbegin(); it != w.body.rend(); ++it)
         (*it)->eraseFromParent();
   }
   return {static_cast<unsigned>(windows.size()), true};
}
} // namespace sirius_product17
