; ModuleID = 'jit_size_function'
source_filename = "jit_size_function"
target datalayout = "e-m:e-p270:32:32-p271:32:32-p272:64:64-i64:64-i128:128-f80:128-n8:16:32:64-S128"

declare ptr @coro_malloc(i32)

declare void @coro_free(ptr)

define { <8 x i32>, <8 x i32>, <8 x i32>, <8 x i32> } @size(i64 %0, <8 x i32> %1) {
entry:
  %2 = add i64 %0, 56
  %3 = inttoptr i64 %2 to ptr
  %4 = load i64, ptr %3, align 8
  %5 = add i64 %4, 288
  %6 = inttoptr i64 %5 to ptr
  %7 = load i64, ptr %6, align 8
  %8 = add i64 %7, 132168
  %9 = inttoptr i64 %8 to ptr
  %10 = load ptr, ptr %9, align 8
  %11 = call i64 %10(i64 %7, i64 %4, i32 0)
  %12 = inttoptr i64 %11 to ptr
  %13 = call { <8 x i32>, <8 x i32>, <8 x i32>, <8 x i32> } %12(i64 %0, <8 x i32> %1)
  ret { <8 x i32>, <8 x i32>, <8 x i32>, <8 x i32> } %13
}
