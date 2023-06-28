const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;
const Allocator = std.mem.Allocator;

// TODO:
// ? keep track of inverse permutations
//   maybe second datastructure ?

/// return struct that can create Permutations
pub const Permutation = struct {
    val: [*]usize,
    len: usize,

    fn init(n: usize, allocator: Allocator) !Permutation {
        return Permutation{
            .val = (try allocator.alloc(usize, n)).ptr,
            .len = n,
        };
    }

    /// deinitialize Permutation
    pub fn deinit(a: Permutation, allocator: Allocator) void {
        allocator.free(a.val[0..a.len]);
    }

    /// allocates a vector with n elements with value a
    pub fn eye(n: usize, allocator: Allocator) !Permutation {
        var res = try init(n, allocator);
        for (0..res.len) |i| {
            res.val[i] = i;
        }
        return res;
    }

    /// allocates a copy of permutation a
    pub fn copy(a: Permutation, allocator: Allocator) !Permutation {
        return Permutation{
            .val = (try allocator.dupe(usize, a.val[0..a.len])).ptr,
            .len = a.len,
        };
    }

    /// create inverse permuatation
    pub fn inv(a: Permutation, allocator: Allocator) !Permutation {
        var res = try init(a.len, allocator);
        for (0..res.len) |i| {
            res.val[a.at(i)] = i;
        }
        return res;
    }

    /// returns target at index
    pub fn at(a: Permutation, i: usize) usize {
        assert(i < a.len);
        return a.val[i];
    }

    /// swaps the two elements at i and j in the permutation.
    /// The operation is performed on a.
    pub fn swap(a: Permutation, i: usize, j: usize) Permutation {
        assert(i < a.len);
        assert(j < a.len);
        const h = a.val[i];
        a.val[i] = a.val[j];
        a.val[j] = h;
        return a;
    }

    /// set a to: a applied to b
    pub fn applyTo(a: Permutation, b: Permutation) Permutation {
        assert(a.len == b.len);
        for (0..a.len) |i| {
            a.val[i] = b.at(a.at(i));
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

    const n = 3;
    var a = (try P.eye(n, ally)).swap(1, 2).swap(0, 1);
    defer a.deinit(ally);
    try testing.expectEqual(@as(usize, 2), a.at(0));
    try testing.expectEqual(@as(usize, 0), a.at(1));
    try testing.expectEqual(@as(usize, 1), a.at(2));

    var b = (try P.eye(n, ally)).swap(0, 2).applyTo(a);
    defer b.deinit(ally);
    _ = a.swap(0, 2);
    try testing.expectEqual(a.at(0), b.at(0));
    try testing.expectEqual(a.at(1), b.at(1));
    try testing.expectEqual(a.at(2), b.at(2));

    const c = (try P.inv(a, ally)).applyTo(a);
    defer c.deinit(ally);
    try testing.expectEqual(@as(usize, 0), c.at(0));
    try testing.expectEqual(@as(usize, 1), c.at(1));
    try testing.expectEqual(@as(usize, 2), c.at(2));
}
