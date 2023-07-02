const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;
const Allocator = std.mem.Allocator;

/// Permutation struct:
/// keeps track of the inverse
pub const Permutation = struct {
    val: [*]usize,
    inv: [*]usize,
    len: usize,

    //initializes Permutation of size n with undefined values
    fn init(n: usize, allocator: Allocator) !Permutation {
        return Permutation{
            .val = (try allocator.alloc(usize, n)).ptr,
            .inv = (try allocator.alloc(usize, n)).ptr,
            .len = n,
        };
    }

    /// deinitialize Permutation
    pub fn deinit(a: Permutation, allocator: Allocator) void {
        allocator.free(a.val[0..a.len]);
        allocator.free(a.inv[0..a.len]);
    }

    /// allocate identity permutation of n elements
    pub fn eye(n: usize, allocator: Allocator) !Permutation {
        var res = try init(n, allocator);
        for (0..res.len) |i| {
            res.val[i] = i;
        }
        @memcpy(res.inv[0..n], res.val[0..n]);
        return res;
    }

    /// allocate copy of permutation a
    pub fn copy(a: Permutation, allocator: Allocator) !Permutation {
        return Permutation{
            .val = (try allocator.dupe(usize, a.val[0..a.len])).ptr,
            .inv = (try allocator.dupe(usize, a.inv[0..a.len])).ptr,
            .len = a.len,
        };
    }

    /// allocte inverse permuatation
    pub fn inv(a: Permutation, allocator: Allocator) !Permutation {
        return Permutation{
            .val = (try allocator.dupe(usize, a.inv[0..a.len])).ptr,
            .inv = (try allocator.dupe(usize, a.val[0..a.len])).ptr,
            .len = a.len,
        };
    }

    /// return target at index
    pub fn at(a: Permutation, i: usize) usize {
        assert(i < a.len);
        return a.val[i];
    }

    ///return source of index
    pub fn atInv(a: Permutation, i: usize) usize {
        assert(i < a.len);
        return a.inv[i];
    }

    /// swaps the two elements at i and j in the permutation.
    /// The operation is performed on a.
    pub fn swap(a: Permutation, i: usize, j: usize) Permutation {
        assert(i < a.len);
        assert(j < a.len);
        const value_i = a.val[i];
        const value_j = a.val[j];
        a.val[i] = value_j;
        a.val[j] = value_i;
        a.inv[value_i] = j;
        a.inv[value_j] = i;
        return a;
    }

    /// swaps the two elements at i and j in the initail configuration.
    /// The operation is performed on a.
    pub fn swapInv(a: Permutation, i: usize, j: usize) Permutation {
        assert(i < a.len);
        assert(j < a.len);
        const value_i = a.inv[i];
        const value_j = a.inv[j];
        a.inv[i] = value_j;
        a.inv[j] = value_i;
        a.val[value_i] = j;
        a.val[value_j] = i;
        return a;
    }

    /// set a to: a applied to b
    pub fn applyTo(a: Permutation, b: Permutation) Permutation {
        assert(a.len == b.len);
        for (0..a.len) |i| {
            a.val[i] = b.at(a.val[i]);
            a.inv[a.val[i]] = i;
        }
        return a;
    }
};

test "creation" {
    const ally = std.testing.allocator;
    const P = Permutation;

    const n = 100;
    const p = try P.eye(n, ally);
    defer p.deinit(ally);

    for (0..n) |i| {
        try testing.expectEqual(i, p.at(i));
    }
}

test "manipulation" {
    const ally = std.testing.allocator;
    const P = Permutation;

    const n = 4;
    const a = (try P.eye(n, ally)).swap(1, 2).swap(0, 1);
    defer a.deinit(ally);
    try testing.expectEqual(@as(usize, 2), a.at(0));
    try testing.expectEqual(@as(usize, 0), a.at(1));
    try testing.expectEqual(@as(usize, 1), a.at(2));
    try testing.expectEqual(@as(usize, 3), a.at(3));

    try testing.expectEqual(@as(usize, 1), a.atInv(0));
    try testing.expectEqual(@as(usize, 2), a.atInv(1));
    try testing.expectEqual(@as(usize, 0), a.atInv(2));
    try testing.expectEqual(@as(usize, 3), a.atInv(3));

    const b_ = (try P.eye(n, ally)).swap(0, 3);
    defer b_.deinit(ally);
    const b = (try a.copy(ally)).applyTo(b_);
    defer b.deinit(ally);
    try testing.expectEqual(@as(usize, 2), b.at(0));
    try testing.expectEqual(@as(usize, 3), b.at(1));
    try testing.expectEqual(@as(usize, 1), b.at(2));
    try testing.expectEqual(@as(usize, 0), b.at(3));

    try testing.expectEqual(@as(usize, 3), b.atInv(0));
    try testing.expectEqual(@as(usize, 2), b.atInv(1));
    try testing.expectEqual(@as(usize, 0), b.atInv(2));
    try testing.expectEqual(@as(usize, 1), b.atInv(3));

    const c = (try P.inv(b, ally)).applyTo(b);
    defer c.deinit(ally);
    try testing.expectEqual(@as(usize, 0), c.at(0));
    try testing.expectEqual(@as(usize, 1), c.at(1));
    try testing.expectEqual(@as(usize, 2), c.at(2));
    try testing.expectEqual(@as(usize, 3), c.at(3));
}
