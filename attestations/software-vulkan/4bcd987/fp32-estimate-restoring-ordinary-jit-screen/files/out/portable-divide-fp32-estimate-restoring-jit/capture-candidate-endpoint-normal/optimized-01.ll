; ModuleID = 'jit_size_function'
source_filename = "jit_size_function"
target datalayout = "e-m:e-p270:32:32-p271:32:32-p272:64:64-i64:64-i128:128-f80:128-n8:16:32:64-S128"

declare ptr @coro_malloc(i32)

declare void @coro_free(ptr)

define { <8 x i32>, <8 x i32>, <8 x i32>, <8 x i32> } @size(i64 %0) {
entry:
  %1 = add i64 %0, 56
  %2 = inttoptr i64 %1 to ptr
  %3 = load i64, ptr %2, align 8
  %4 = add i64 %3, 288
  %5 = inttoptr i64 %4 to ptr
  %6 = load i64, ptr %5, align 8
  %7 = add i64 %6, 132168
  %8 = inttoptr i64 %7 to ptr
  %9 = load ptr, ptr %8, align 8
  %10 = call i64 %9(i64 %6, i64 %3, i32 1)
  %11 = inttoptr i64 %10 to ptr
  %12 = call { <8 x i32>, <8 x i32>, <8 x i32>, <8 x i32> } %11(i64 %0)
  ret { <8 x i32>, <8 x i32>, <8 x i32>, <8 x i32> } %12
}
