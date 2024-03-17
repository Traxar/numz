const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;
const Allocator = std.mem.Allocator;

/// permutation struct that keeps track of the inverse
/// 'Index' can be any unsigned integer type
pub fn PermutationType(comptime Index: type) type {
    return struct {
        const Permutation = @This();
        val: [*]Index,
        len: Index,

        /// allocate permutation of size n with undefined values
        pub fn init(n: Index, allocator: Allocator) !Permutation {
            return Permutation{
                .val = (try allocator.alloc(Index, 2 * n)).ptr,
                .len = n,
            };
        }

        ///allocate permutation with undefined values and same dimensions as a
        pub fn like(a: Permutation, allocator: Allocator) !Permutation {
            return init(a.len, allocator);
        }

        /// deinitialize permutation
        pub fn deinit(a: Permutation, allocator: Allocator) void {
            allocator.free(a.val[0 .. 2 * a.len]);
        }

        /// set permutation to identity
        pub fn eye(res: Permutation) void {
            var i: Index = 0;
            const n = res.len;
            while (i < n) : (i += 1) {
                res.val[i] = i;
            }
            @memcpy(res.val[n .. 2 * n], res.val[0..n]);
        }

        /// res <- a
        pub fn copy(a: Permutation, res: Permutation) void {
            assert(a.len == res.len);
            const n = res.len;
            @memcpy(res.val[0 .. 2 * n], a.val[0 .. 2 * n]);
        }

        /// res <- a^-1
        pub fn inv(a: Permutation, res: Permutation) void {
            assert(a.len == res.len);
            assert(a.val != res.val); //inplace
            const n = res.len;
            @memcpy(res.val[0..n], a.val[n .. 2 * n]);
            @memcpy(res.val[n .. 2 * n], a.val[0..n]);
        }

        /// return target at index
        pub inline fn at(a: Permutation, i: Index) Index {
            return a.val[i];
        }

        ///return source of index
        pub inline fn atInv(a: Permutation, i: Index) Index {
            return a.val[a.len + i];
        }

        /// swaps the two elements at i and j in the permutation.
        /// The operation is performed on a.
        pub fn swap(a: Permutation, i: Index, j: Index) void {
            const value_i = a.val[i];
            const value_j = a.val[j];
            a.val[i] = value_j;
            a.val[j] = value_i;
            a.val[a.len + value_i] = j;
            a.val[a.len + value_j] = i;
        }

        /// swaps the two elements at i and j in the initial configuration.
        /// The operation is performed on a.
        pub fn swapInv(a: Permutation, i: Index, j: Index) Permutation {
            const value_i = a.val[a.len + i];
            const value_j = a.val[a.len + j];
            a.val[a.len + i] = value_j;
            a.val[a.len + j] = value_i;
            a.val[value_i] = j;
            a.val[value_j] = i;
            return a;
        }

        /// a <- a(b(.))
        pub fn applyTo(a: Permutation, b: Permutation) void {
            assert(a.len == b.len);
            for (0..a.len) |i| {
                a.val[i] = b.at(a.val[i]);
                a.val[a.len + a.val[i]] = i;
            }
        }

        /// a <- b(a(.))
        pub fn apply(a: Permutation, b: Permutation) void {
            assert(a.len == b.len);
            for (0..a.len) |i| {
                const j = a.len + i;
                a.val[j] = b.atInv(a.val[j]);
                a.val[a.val[j]] = i;
            }
        }
    };
}

test "permutation creation" {
    const ally = std.testing.allocator;
    const P = PermutationType(usize);

    const n = 100;
    const p = try P.init(n, ally);
    defer p.deinit(ally);
    p.eye();

    for (0..n) |i| {
        try testing.expectEqual(i, p.at(i));
    }
}

test "permutation operations" {
    const ally = std.testing.allocator;
    const P = PermutationType(usize);

    const n = 4;
    var a = try P.init(n, ally);
    defer a.deinit(ally);
    a.eye();
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

    var b = try P.like(a, ally);
    defer b.deinit(ally);
    b.eye();
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

    var c = try P.like(b, ally);
    defer c.deinit(ally);
    b.inv(c);
    c.applyTo(b);
    try testing.expectEqual(@as(usize, 0), c.at(0));
    try testing.expectEqual(@as(usize, 1), c.at(1));
    try testing.expectEqual(@as(usize, 2), c.at(2));
    try testing.expectEqual(@as(usize, 3), c.at(3));
}
