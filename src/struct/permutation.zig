const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;
const Allocator = std.mem.Allocator;

///Permutation struct
/// 'Index' can be any unsigned integer type
pub fn PermutationType(comptime Index: type) type {
    return struct {
        const Permutation = @This();
        val: [*]Index,
        len: Index,

        /// allocate Permutation of size n with undefined values
        fn init(n: Index, allocator: Allocator) !Permutation {
            return Permutation{
                .val = (try allocator.alloc(Index, n)).ptr,
                .len = n,
            };
        }

        pub fn deinit(a: Permutation, allocator: Allocator) void {
            allocator.free(a.val[0..a.len]);
        }

        /// allocate identity Permutation of size n
        pub fn eye(n: Index, allocator: Allocator) !Permutation {
            var res = try init(n, allocator);
            for (0..res.len) |i| {
                res.val[i] = i;
            }
            return res;
        }

        /// allocate copy of Permutation a
        pub fn copy(a: Permutation, allocator: Allocator) !Permutation {
            return Permutation{
                .val = (try allocator.dupe(Index, a.val[0..a.len])).ptr,
                .len = a.len,
            };
        }

        /// return target at index
        pub inline fn at(a: Permutation, i: Index) Index {
            return a.val[i];
        }

        /// swaps the two elements at i and j in the Permutation.
        pub fn swap(a: Permutation, i: Index, j: Index) void {
            const value_i = a.val[i];
            const value_j = a.val[j];
            a.val[i] = value_j;
            a.val[j] = value_i;
        }

        /// a <- a(b(.))
        pub fn applyTo(a: Permutation, b: Permutation) void {
            assert(a.len == b.len);
            for (0..a.len) |i| {
                a.val[i] = b.at(a.val[i]);
            }
        }
    };
}

/// Permutation struct that keeps track of the inverse
/// 'Index' can be any unsigned integer type
pub fn BiPermutationType(comptime Index: type) type {
    return struct {
        const BiPermutation = @This();
        const Permutation = PermutationType(Index);
        val: [*]Index,
        inv: [*]Index,
        len: Index,

        /// allocate BiPermutation of size n with undefined values
        fn init(n: Index, allocator: Allocator) !BiPermutation {
            return BiPermutation{
                .val = (try allocator.alloc(Index, n)).ptr,
                .inv = (try allocator.alloc(Index, n)).ptr,
                .len = n,
            };
        }

        /// deinitialize BiPermutation
        pub fn deinit(a: BiPermutation, allocator: Allocator) void {
            allocator.free(a.val[0..a.len]);
            allocator.free(a.inv[0..a.len]);
        }

        pub fn deinitToPermutation(a: BiPermutation, allocator: Allocator) Permutation {
            allocator.free(a.inv[0..a.len]);
            return Permutation{
                .val = a.val,
                .len = a.len,
            };
        }

        pub fn deinitToInvPermutation(a: BiPermutation, allocator: Allocator) Permutation {
            allocator.free(a.val[0..a.len]);
            return Permutation{
                .val = a.inv,
                .len = a.len,
            };
        }

        /// allocate identity BiPermutation of size n
        pub fn eye(n: Index, allocator: Allocator) !BiPermutation {
            var res = try init(n, allocator);
            for (0..res.len) |i| {
                res.val[i] = i;
            }
            @memcpy(res.inv[0..n], res.val[0..n]);
            return res;
        }

        /// allocate copy of BiPermutation a
        pub fn copy(a: BiPermutation, allocator: Allocator) !BiPermutation {
            return BiPermutation{
                .val = (try allocator.dupe(Index, a.val[0..a.len])).ptr,
                .inv = (try allocator.dupe(Index, a.inv[0..a.len])).ptr,
                .len = a.len,
            };
        }

        /// allocte inverse permuatation
        pub fn inv(a: BiPermutation, allocator: Allocator) !BiPermutation {
            return BiPermutation{
                .val = (try allocator.dupe(Index, a.inv[0..a.len])).ptr,
                .inv = (try allocator.dupe(Index, a.val[0..a.len])).ptr,
                .len = a.len,
            };
        }

        /// return target at index
        pub inline fn at(a: BiPermutation, i: Index) Index {
            return a.val[i];
        }

        ///return source of index
        pub inline fn atInv(a: BiPermutation, i: Index) Index {
            return a.inv[i];
        }

        /// swaps the two elements at i and j in the BiPermutation.
        /// The operation is performed on a.
        pub fn swap(a: BiPermutation, i: Index, j: Index) void {
            const value_i = a.val[i];
            const value_j = a.val[j];
            a.val[i] = value_j;
            a.val[j] = value_i;
            a.inv[value_i] = j;
            a.inv[value_j] = i;
        }

        /// swaps the two elements at i and j in the initail configuration.
        /// The operation is performed on a.
        pub fn swapInv(a: BiPermutation, i: Index, j: Index) BiPermutation {
            const value_i = a.inv[i];
            const value_j = a.inv[j];
            a.inv[i] = value_j;
            a.inv[j] = value_i;
            a.val[value_i] = j;
            a.val[value_j] = i;
            return a;
        }

        /// a <- a(b(.))
        pub fn applyTo(a: BiPermutation, b: BiPermutation) void {
            assert(a.len == b.len);
            for (0..a.len) |i| {
                a.val[i] = b.at(a.val[i]);
                a.inv[a.val[i]] = i;
            }
        }

        /// a <- b(a(.))
        pub fn apply(a: BiPermutation, b: BiPermutation) void {
            assert(a.len == b.len);
            for (0..a.len) |i| {
                a.inv[i] = b.atInv(a.inv[i]);
                a.val[a.inv[i]] = i;
            }
        }
    };
}

test "permutation creation" {
    const ally = std.testing.allocator;
    const P = PermutationType(usize);

    const n = 100;
    const p = try P.eye(n, ally);
    defer p.deinit(ally);

    for (0..n) |i| {
        try testing.expectEqual(i, p.at(i));
    }
}

test "permutation operations" {
    const ally = std.testing.allocator;
    const P = PermutationType(usize);

    const n = 4;
    var a = try P.eye(n, ally);
    defer a.deinit(ally);
    a.swap(1, 2);
    a.swap(0, 1);
    try testing.expectEqual(@as(usize, 2), a.at(0));
    try testing.expectEqual(@as(usize, 0), a.at(1));
    try testing.expectEqual(@as(usize, 1), a.at(2));
    try testing.expectEqual(@as(usize, 3), a.at(3));

    var b = try P.eye(n, ally);
    defer b.deinit(ally);
    b.swap(0, 3);
    b.applyTo(a);
    try testing.expectEqual(@as(usize, 3), b.at(0));
    try testing.expectEqual(@as(usize, 0), b.at(1));
    try testing.expectEqual(@as(usize, 1), b.at(2));
    try testing.expectEqual(@as(usize, 2), b.at(3));
}

test "bipermutation creation" {
    const ally = std.testing.allocator;
    const P = BiPermutationType(usize);

    const n = 100;
    const p = try P.eye(n, ally);
    defer p.deinit(ally);

    for (0..n) |i| {
        try testing.expectEqual(i, p.at(i));
    }
}

test "bipermutation operations" {
    const ally = std.testing.allocator;
    const P = BiPermutationType(usize);

    const n = 4;
    var a = try P.eye(n, ally);
    defer a.deinit(ally);
    a.swap(1, 2);
    a.swap(0, 1);
    try testing.expectEqual(@as(usize, 2), a.at(0));
    try testing.expectEqual(@as(usize, 0), a.at(1));
    try testing.expectEqual(@as(usize, 1), a.at(2));
    try testing.expectEqual(@as(usize, 3), a.at(3));

    try testing.expectEqual(@as(usize, 1), a.atInv(0));
    try testing.expectEqual(@as(usize, 2), a.atInv(1));
    try testing.expectEqual(@as(usize, 0), a.atInv(2));
    try testing.expectEqual(@as(usize, 3), a.atInv(3));

    var b = try P.eye(n, ally);
    defer b.deinit(ally);
    b.swap(0, 3);
    b.apply(a);
    try testing.expectEqual(@as(usize, 2), b.at(0));
    try testing.expectEqual(@as(usize, 3), b.at(1));
    try testing.expectEqual(@as(usize, 1), b.at(2));
    try testing.expectEqual(@as(usize, 0), b.at(3));

    try testing.expectEqual(@as(usize, 3), b.atInv(0));
    try testing.expectEqual(@as(usize, 2), b.atInv(1));
    try testing.expectEqual(@as(usize, 0), b.atInv(2));
    try testing.expectEqual(@as(usize, 1), b.atInv(3));

    var c = try P.inv(b, ally);
    defer c.deinit(ally);
    c.applyTo(b);
    try testing.expectEqual(@as(usize, 0), c.at(0));
    try testing.expectEqual(@as(usize, 1), c.at(1));
    try testing.expectEqual(@as(usize, 2), c.at(2));
    try testing.expectEqual(@as(usize, 3), c.at(3));
}
