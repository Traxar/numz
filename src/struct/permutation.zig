const std = @import("std");
const testing = std.testing;

/// vector datastructure
pub fn Permutation(comptime size: usize) type {
    return struct {
        const Perm = @This();
        pub const n = size;
        p: [n]usize,

        /// Identity
        pub const eye = eye: {
            var e = Perm{ .p = undefined };
            for (0..Perm.n) |i| {
                e.p[i] = i;
            }
            break :eye e;
        };

        /// returns target at index
        pub fn at(a: Perm, i: usize) usize {
            return a.p[i];
        }

        /// swaps the two elements at i and j in the permutation.
        /// The operation is performed on a itself.
        pub fn swap(a: *Perm, i: usize, j: usize) void {
            const h = a.p[i];
            a.p[i] = a.p[j];
            a.p[j] = h;
        }

        /// returns b applied to a
        pub fn apply(a: Perm, b: Perm) Perm {
            var res = Perm{ .p = undefined };
            for (0..Perm.n) |i| {
                res.p[i] = a.p[b.p[i]];
            }
            return res;
        }

        /// returns inverse permuatation
        pub fn inv(a: Perm) Perm {
            var res = Perm{ .p = undefined };
            for (0..Perm.n) |i| {
                res.p[a.p[i]] = i;
            }
            return res;
        }

        /// returns true if equal
        pub fn cmp(a: Perm, b: Perm) bool {
            for (0..Perm.n - 1) |i| {
                if (a.p[i] != b.p[i]) {
                    return false;
                }
            }
            return true;
        }
    };
}

test "creation" {
    const n = 100;
    const P = Permutation(n);

    for (0..n) |i| {
        try testing.expectEqual(i, P.eye.at(i));
    }
}

test "manipulation" {
    const n = 3;
    const P = Permutation(n);

    var a = P.eye; //012
    a.swap(1, 2); //021
    a.swap(0, 1); //201
    try testing.expectEqual(@as(usize, 2), a.at(0));
    try testing.expectEqual(@as(usize, 0), a.at(1));
    try testing.expectEqual(@as(usize, 1), a.at(2));

    var b = P.eye;
    b.swap(0, 2); //210
    const c = a.apply(b); //102
    a.swap(0, 2); //102
    try testing.expectEqual(a.at(0), c.at(0));
    try testing.expectEqual(a.at(1), c.at(1));
    try testing.expectEqual(a.at(2), c.at(2));

    const d = c.apply(c.inv()); //012
    try testing.expect(P.eye.cmp(d));
}
