; ModuleID = 'jit_fetch_function'
source_filename = "jit_fetch_function"
target datalayout = "e-m:e-p270:32:32-p271:32:32-p272:64:64-i64:64-i128:128-f80:128-n8:16:32:64-S128"

declare ptr @coro_malloc(i32)

declare void @coro_free(ptr)

define { <8 x float>, <8 x float>, <8 x float>, <8 x float>, <8 x i32> } @fetch(i64 %0, i64 %1, <8 x i32> %2, <8 x i32> %3, <8 x i32> %4, <8 x i32> %5, <8 x i32> %6) {
entry:
  %7 = add i64 %0, 56
  %8 = inttoptr i64 %7 to ptr
  %9 = load i64, ptr %8, align 8
  %10 = add i64 %9, 288
  %11 = inttoptr i64 %10 to ptr
  %12 = load i64, ptr %11, align 8
  %13 = add i64 %12, 132144
  %14 = inttoptr i64 %13 to ptr
  %15 = load ptr, ptr %14, align 8
  %16 = call i64 %15(i64 %12, i64 %9, i32 36)
  %17 = inttoptr i64 %16 to ptr
  %18 = call { <8 x float>, <8 x float>, <8 x float>, <8 x float>, <8 x i32> } %17(i64 %0, i64 %1, <8 x i32> %2, <8 x i32> %3, <8 x i32> %4, <8 x i32> %5, <8 x i32> %6)
  ret { <8 x float>, <8 x float>, <8 x float>, <8 x float>, <8 x i32> } %18
}
