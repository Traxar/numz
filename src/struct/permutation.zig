const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;

/// TODO: keep track of inverse Permutation
/// return struct that can create Permutations
pub const PermutationType = struct {
    const Self = @This();
    const ArrayList = std.ArrayList(usize);
    allocator: std.mem.Allocator,

    fn init(self: Self, n: usize) !Permutation {
        return Permutation{
            .val = (try self.allocator.alloc(usize, n)).ptr,
            .len = n,
        };
    }

    /// deinitialize Permutation
    pub fn deinit(self: Self, a: Permutation) void {
        self.allocator.free(a.val[0..a.len]);
    }

    /// allocates a vector with n elements with value a
    pub fn eye(self: Self, n: usize) !Permutation {
        var res = try self.init(n);
        for (0..res.len) |i| {
            res.val[i] = i;
        }
        return res;
    }

    /// allocates a copy of permutation a
    pub fn copy(self: Self, a: Permutation) !Permutation {
        return Permutation{
            .val = (try self.allocator.dupe(usize, a.val[0..a.len])).ptr,
            .len = a.len,
        };
    }

    /// create inverse permuatation
    pub fn inv(self: Self, a: Permutation) !Permutation {
        var res = try self.init(a.len);
        for (0..res.len) |i| {
            res.val[a.at(i)] = i;
        }
        return res;
    }

    /// permutation with runtime length
    /// all functions in this struct are inplace and do not require memory allocation
    const Permutation = struct {
        val: [*]usize,
        len: usize,

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
};

test "creation" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();
    defer {
        const deinit_status = gpa.deinit();
        //fail test; can't try in defer as defer is executed after we return
        if (deinit_status == .leak) testing.expect(false) catch @panic("TEST FAIL");
    }

    const P = PermutationType{ .allocator = allocator };

    const n = 100;
    const p = try P.eye(n);
    defer P.deinit(p);

    for (0..n) |i| {
        try testing.expectEqual(i, p.at(i));
    }
}

test "manipulation" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();
    defer {
        const deinit_status = gpa.deinit();
        //fail test; can't try in defer as defer is executed after we return
        if (deinit_status == .leak) testing.expect(false) catch @panic("TEST FAIL");
    }

    const P = PermutationType{ .allocator = allocator };

    const n = 3;
    var a = (try P.eye(n)).swap(1, 2).swap(0, 1);
    defer P.deinit(a);
    try testing.expectEqual(@as(usize, 2), a.at(0));
    try testing.expectEqual(@as(usize, 0), a.at(1));
    try testing.expectEqual(@as(usize, 1), a.at(2));

    var b = (try P.eye(n)).swap(0, 2).applyTo(a);
    defer P.deinit(b);
    _ = a.swap(0, 2);
    try testing.expectEqual(a.at(0), b.at(0));
    try testing.expectEqual(a.at(1), b.at(1));
    try testing.expectEqual(a.at(2), b.at(2));

    const c = (try P.inv(a)).applyTo(a);
    defer P.deinit(c);
    try testing.expectEqual(@as(usize, 0), c.at(0));
    try testing.expectEqual(@as(usize, 1), c.at(1));
    try testing.expectEqual(@as(usize, 2), c.at(2));
}
