const std = @import("std");
const testing = std.testing;

/// vector datastructure
pub fn Vector(comptime scalar: type, comptime size: usize) type {
    return struct {
        const Vec = @This();
        const Scalar = scalar;
        pub const n = size;
        v: [n]Scalar,

        /// returns vector filled with s
        pub fn fill(s: Scalar) Vec {
            return Vec{ .v = [_]Scalar{s} ** size };
        }

        /// returns value at index
        pub fn at(a: Vec, i: usize) Scalar {
            return a.v[i];
        }

        /// sets the value at index
        pub fn set(a: *Vec, i: usize, b: Scalar) void {
            a.v[i] = b;
        }

        /// returns a + b
        pub fn add(a: Vec, b: Vec) Vec {
            var res = Vec{ .v = undefined };
            for (a.v, b.v, 0..) |a_, b_, i| {
                res.v[i] = a_.add(b_);
            }
            return res;
        }

        /// returns a - b
        pub fn sub(a: Vec, b: Vec) Vec {
            var res = Vec{ .v = undefined };
            for (a.v, b.v, 0..) |a_, b_, i| {
                res.v[i] = a_.sub(b_);
            }
            return res;
        }

        ///returns elementwise a * b
        pub fn mul(a: Vec, b: Vec) Vec {
            var res = Vec{ .v = undefined };
            for (a.v, b.v, 0..) |a_, b_, i| {
                res.v[i] = a_.mul(b_);
            }
            return res;
        }

        ///returns elementwise a / b
        pub fn div(a: Vec, b: Vec) Vec {
            var res = Vec{ .v = undefined };
            for (a.v, b.v, 0..) |a_, b_, i| {
                res.v[i] = a_.div(b_);
            }
            return res;
        }

        ///returns elementwise a / b
        pub fn dot(a: Vec, b: Vec) Scalar {
            var res = Scalar.zero;
            for (a.v, b.v) |a_, b_| {
                res = res.add(a_.mul(b_));
            }
            return res;
        }

        ///returns the sum of elements
        pub fn sum(a: Vec) Scalar {
            var res = Scalar.zero;
            for (a.v) |a_| {
                res = res.add(a_);
            }
            return res;
        }

        ///returns the euklidean norm
        pub fn norm(a: Vec) Scalar {
            return a.dot(a).sqrt();
        }

        ///returns the smallest element
        pub fn min(a: Vec) Scalar {
            var res = a.v[0];
            for (1..Vec.n) |i| {
                res = res.min(a.v[i]);
            }
            return res;
        }

        ///returns the largest element
        pub fn max(a: Vec) Scalar {
            var res = a.v[0];
            for (1..Vec.n) |i| {
                res = res.max(a.v[i]);
            }
            return res;
        }
    };
}

test "creation" {
    const n = 100;
    const F = @import("../scalar.zig").Float(f32);
    const V = Vector(F, n);

    const a = F.from(-3.14);
    const v = V.fill(a);

    try testing.expectEqual([_]F{a} ** n, v.v);
    try testing.expect(a.cmp(.equal, v.at(0)));
}

test "operators" {
    const F = @import("../scalar.zig").Float(f32);
    const V = Vector(F, 3);

    const a_ = F.from(-3.14);
    const b_ = F.from(5.27);
    var a = V.fill(a_);
    a.set(1, F.zero);
    var b = V.fill(b_);
    b.set(0, F.eye);

    try testing.expect(a.at(1).cmp(.equal, F.zero));

    const c = a.add(b);
    try testing.expect(c.at(0).cmp(.equal, a_.add(F.eye)));
    try testing.expect(c.at(1).cmp(.equal, b_));
    try testing.expect(c.at(2).cmp(.equal, a_.add(b_)));

    const d = a.sub(b);
    try testing.expect(d.at(0).cmp(.equal, a_.sub(F.eye)));
    try testing.expect(d.at(1).cmp(.equal, F.zero.sub(b_)));
    try testing.expect(d.at(2).cmp(.equal, a_.sub(b_)));

    const e = a.mul(b);
    try testing.expect(e.at(0).cmp(.equal, a_));
    try testing.expect(e.at(1).cmp(.equal, F.zero));
    try testing.expect(e.at(2).cmp(.equal, a_.mul(b_)));

    const f = a.div(b);
    try testing.expect(f.at(0).cmp(.equal, a_));
    try testing.expect(f.at(1).cmp(.equal, F.zero));
    try testing.expect(f.at(2).cmp(.equal, a_.div(b_)));

    try testing.expect(a.dot(b).cmp(.equal, a_.add(a_.mul(b_))));

    try testing.expect(b.sum().cmp(.equal, b_.add(b_).add(F.eye)));

    const a_2 = a_.mul(a_);
    try testing.expect(a.norm().cmp(.equal, a_2.add(a_2).sqrt()));

    try testing.expect(a.min().cmp(.equal, a_.min(F.zero)));
    try testing.expect(b.max().cmp(.equal, b_.max(F.eye)));
}
