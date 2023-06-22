const std = @import("std");
pub const Error = std.mem.Allocator.Error || error{UndefinedMath};
//TODO:
// ? currently unused, maybe define all operations on example Error!Vector
//  + avoid needing try statements in the calculations
//  - error handling is no longer enforced by heart
//    could result in huge error traces, the error gets passed around until the next try by the user
//  ? use this only on inline operations add, mul, ...
//    do not use on deinit, cmp, LU, print, set, at ...
//   - inline operations might allocate new instances which later cannot be freed unless not using inline in the frist place
