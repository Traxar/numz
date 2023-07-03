const std = @import("std");
const testing = std.testing;
const assert = std.debug.assert;
const CompareOperator = std.math.CompareOperator;
const Allocator = std.mem.Allocator;

// TODO:
// ~ change from to initalize from a fraction

/// wrapper for the float datastructure
/// choose float from {f16,f32,f64,f80,f128}
pub fn Float(comptime float: type) type {
    return struct {
        const Scalar = @This();
        f: float,

        /// the 0 element
        pub const zero = Scalar{ .f = 0 };

        /// the 1 element
        pub const eye = Scalar{ .f = 1 };

        /// return element isomorph to f
        pub fn from(p: isize, q: usize) Scalar {
            assert(q != 0);
            return Scalar{ .f = @as(float, @floatFromInt(p)) / @as(float, @floatFromInt(q)) };
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

        /// returns truth value of: a r b
        pub fn cmp(a: Scalar, r: CompareOperator, b: Scalar) bool {
            return std.math.order(a.f, b.f).compare(r);
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

        try testing.expectEqual(@as(f, -3.14), F.from(-314, 100).f);
    }
}

test "operators" {
    const fTypes = [_]type{ f16, f32, f64, f80, f128 };
    inline for (fTypes) |f| {
        const F = Float(f);
        const q = 100;
        const a_ = -314;
        const b_ = 527;
        const a = F.from(a_, q);
        const b = F.from(b_, q);
        try testing.expectEqual(a.f + b.f, a.add(b).f);
        try testing.expectEqual(a.f - b.f, a.sub(b).f);
        try testing.expectEqual(b.f - a.f, b.sub(a).f);
        try testing.expectEqual(a.f * b.f, a.mul(b).f);
        try testing.expectEqual(a.f / b.f, a.div(b).f);
        try testing.expectEqual(b.f / a.f, b.div(a).f);
        try testing.expectEqual(@fabs(a.f), a.abs().f);
        try testing.expectEqual(@sqrt(b.f), b.sqrt().f);
    }
}

test "comparisons" {
    const fTypes = [_]type{ f16, f32, f64, f80, f128 };
    inline for (fTypes) |f| {
        const F = Float(f);
        const a = F.from(-314, 100);
        const b = F.from(527, 100);
        try testing.expect(a.cmp(.lt, b));
        try testing.expect(a.cmp(.lte, b));
        try testing.expect(!a.cmp(.eq, b));
        try testing.expect(!a.cmp(.gte, b));
        try testing.expect(!a.cmp(.gt, b));

        try testing.expect(!b.cmp(.lt, a));
        try testing.expect(!b.cmp(.lte, a));
        try testing.expect(!b.cmp(.eq, a));
        try testing.expect(b.cmp(.gte, a));
        try testing.expect(b.cmp(.gt, a));

        try testing.expect(!a.cmp(.lt, a));
        try testing.expect(a.cmp(.lte, a));
        try testing.expect(a.cmp(.eq, a));
        try testing.expect(a.cmp(.gte, a));
        try testing.expect(!a.cmp(.gt, a));

        try testing.expect(a.min(b).cmp(.eq, a));
        try testing.expect(a.max(b).cmp(.eq, b));
    }
}
