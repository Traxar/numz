const std = @import("std");
const testing = std.testing;
const Relation = @import("relation.zig").Relation;

/// wrapper for the float datastructure
/// choose float from {f16,f32,f64,f80,f128}
/// TODO: f128 is broken due to https://github.com/ziglang/zig/issues/15611
pub fn Float(comptime float: type) type {
    return struct {
        const Scalar = @This();
        f: float,

        /// the 0 element
        pub const zero = Scalar{ .f = 0 };

        /// the 1 element
        pub const eye = Scalar{ .f = 1 };

        /// return element isomorph to f
        pub fn from(f: float) Scalar {
            return Scalar{ .f = f };
        }

        /// returns a + b
        pub fn add(a: Scalar, b: Scalar) Scalar {
            return Scalar{ .f = a.f + b.f };
        }

        /// returns a - b
        pub inline fn sub(a: Scalar, b: Scalar) Scalar {
            return Scalar{ .f = a.f - b.f };
        }

        /// returns a * b
        pub fn mul(a: Scalar, b: Scalar) Scalar {
            return Scalar{ .f = a.f * b.f };
        }

        /// returns a / b
        pub fn div(a: Scalar, b: Scalar) Scalar {
            return Scalar{ .f = a.f / b.f };
        }

        /// returns |a|
        pub fn abs(a: Scalar) Scalar {
            return Scalar{ .f = @fabs(a.f) };
        }

        /// returns |√a|
        pub fn sqrt(a: Scalar) Scalar {
            return Scalar{ .f = @sqrt(a.f) };
        }

        /// returns a R b
        pub fn cmp(a: Scalar, R: Relation, b: Scalar) bool {
            if (a.f < b.f) {
                return Relation.less.imp(R);
            } else if (a.f == b.f) {
                return Relation.equal.imp(R);
            } else {
                return Relation.greater.imp(R);
            }
        }

        ///returns the minimum of a and b
        pub fn min(a: Scalar, b: Scalar) Scalar {
            return Scalar{ .f = @min(a.f, b.f) };
        }

        ///returns the maximum of a and b
        pub fn max(a: Scalar, b: Scalar) Scalar {
            return Scalar{ .f = @max(a.f, b.f) };
        }
    };
}

test "creation" {
    const fTypes = [_]type{ f16, f32, f64, f80, f128 };
    inline for (fTypes) |f| {
        const F = Float(f);

        try testing.expectEqual(@as(f, 0), F.zero.f);
        try testing.expectEqual(@as(f, 1), F.eye.f);

        const x: f = -3.14;
        try testing.expectEqual(x, F.from(x).f);
    }
}

test "operators" {
    const fTypes = [_]type{ f16, f32, f64, f80, f128 };
    inline for (fTypes) |f| {
        const F = Float(f);
        const a_: f = -3.14;
        const b_: f = 5.27;
        const a = F.from(a_);
        const b = F.from(b_);
        try testing.expectEqual(F.from(a_ + b_), a.add(b));
        try testing.expectEqual(F.from(a_ - b_), a.sub(b));
        try testing.expectEqual(F.from(b_ - a_), b.sub(a));
        try testing.expectEqual(F.from(a_ * b_), a.mul(b));
        try testing.expectEqual(F.from(a_ / b_), a.div(b));
        try testing.expectEqual(F.from(b_ / a_), b.div(a));
        try testing.expectEqual(F.from(@fabs(a_)), a.abs());
        try testing.expectEqual(F.from(@sqrt(b_)), b.sqrt());
    }
}

test "comparisons" {
    const fTypes = [_]type{ f16, f32, f64, f80 }; //, f128 }; // TODO: uncomment when https://github.com/ziglang/zig/issues/15611 gets fixed
    inline for (fTypes) |f| {
        const F = Float(f);
        const a = F.from(-3.14);
        const b = F.from(5.27);
        try testing.expect(a.cmp(.less, b));
        try testing.expect(a.cmp(.lessEqual, b));
        try testing.expect(!a.cmp(.equal, b));
        try testing.expect(!a.cmp(.greaterEqual, b));
        try testing.expect(!a.cmp(.greater, b));

        try testing.expect(!b.cmp(.less, a));
        try testing.expect(!b.cmp(.lessEqual, a));
        try testing.expect(!b.cmp(.equal, a));
        try testing.expect(b.cmp(.greaterEqual, a));
        try testing.expect(b.cmp(.greater, a));

        try testing.expect(!a.cmp(.less, a));
        try testing.expect(a.cmp(.lessEqual, a));
        try testing.expect(a.cmp(.equal, a));
        try testing.expect(a.cmp(.greaterEqual, a));
        try testing.expect(!a.cmp(.greater, a));

        try testing.expect(a.min(b).cmp(.equal, a));
        try testing.expect(a.max(b).cmp(.equal, b));
    }
}
