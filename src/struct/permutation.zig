const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;

/// return struct that can create Permutations
pub const PermutationType = struct {
    const Self = @This();
    const ArrayList = std.ArrayList(usize);
    allocator: std.mem.Allocator,

    /// create a permutation with n elements
    /// do not forget to call deinit() on result
    pub fn eye(self: Self, n: usize) !Permutation {
        var res = Permutation{ .p = try ArrayList.initCapacity(self.allocator, n) };
        for (0..n) |i| {
            res.p.appendAssumeCapacity(i);
        }
        return res;
    }

    /// vector with runtime length
    const Permutation = struct {
        p: ArrayList,

        /// deinitialize Vector
        pub fn deinit(a: Permutation) void {
            a.p.deinit();
        }

        /// create a copy of this permutation, using the same allocator
        pub fn copy(a: Permutation) !Permutation {
            return Permutation{ .p = try ArrayList.clone(a.p) };
        }

        /// create inverse permuatation, using the same allocator
        pub fn inv(a: Permutation) !Permutation {
            var res = Permutation{ .p = try ArrayList.initCapacity(a.p.allocator, a.size()) };
            res.p.items.len = a.size();
            for (a.p.items, 0..) |a_, i| {
                res.p.items[a_] = i;
            }
            return res;
        }

        /// returns target at index
        pub fn at(a: Permutation, i: usize) usize {
            return a.p.items[i];
        }

        /// return size
        pub fn size(a: Permutation) usize {
            return a.p.items.len;
        }

        /// swaps the two elements at i and j in the permutation.
        /// The operation is performed on a.
        pub fn swap(a: Permutation, i: usize, j: usize) Permutation {
            const h = a.p.items[i];
            a.p.items[i] = a.p.items[j];
            a.p.items[j] = h;
            return a;
        }

        /// set a to: a applied to b
        pub fn applyTo(a: Permutation, b: Permutation) Permutation {
            assert(a.size() == b.size());
            for (a.p.items, 0..) |a_, i| {
                a.p.items[i] = b.at(a_);
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
    defer p.deinit();

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
    defer a.deinit();
    try testing.expectEqual(@as(usize, 2), a.at(0));
    try testing.expectEqual(@as(usize, 0), a.at(1));
    try testing.expectEqual(@as(usize, 1), a.at(2));

    var b = (try P.eye(n)).swap(0, 2).applyTo(a);
    defer b.deinit();
    _ = a.swap(0, 2);
    try testing.expectEqual(a.at(0), b.at(0));
    try testing.expectEqual(a.at(1), b.at(1));
    try testing.expectEqual(a.at(2), b.at(2));

    const c = (try a.inv()).applyTo(a);
    defer c.deinit();
    try testing.expectEqual(@as(usize, 0), c.at(0));
    try testing.expectEqual(@as(usize, 1), c.at(1));
    try testing.expectEqual(@as(usize, 2), c.at(2));
}
