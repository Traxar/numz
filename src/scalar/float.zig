const std = @import("std");
const testing = std.testing;
const assert = std.debug.assert;
const CompareOperator = std.math.CompareOperator;
const Allocator = std.mem.Allocator;

pub fn Float(comptime float: type) type {
    return SIMDFloat(float, 1);
}

/// wrapper for the float datastructure
/// choose float from {f16,f32,f64,f80,f128}
pub fn SIMDFloat(comptime float: type, comptime n: comptime_int) type {
    return struct {
        const Scalar = @This();
        const len = if (n <= 0) std.simd.suggestVectorSize(float) else n;
        f: @Vector(len, float),

        pub fn SIMDType(comptime length: ?comptime_int) type {
            return SIMDFloat(float, length);
        }

        /// the 0 element
        pub const zero = Scalar{ .f = @splat(0) };

        /// the 1 element
        pub const eye = Scalar{ .f = @splat(1) };

        /// return element isomorph to f
        pub fn from(p: isize, q: usize) Scalar {
            assert(q != 0);
            return Scalar{ .f = @splat(@as(float, @floatFromInt(p)) / @as(float, @floatFromInt(q))) };
        }

        /// returns a + b
        pub inline fn add(a: Scalar, b: Scalar) Scalar {
            return Scalar{ .f = a.f + b.f };
        }

        /// returns a - b
        pub inline fn sub(a: Scalar, b: Scalar) Scalar {
            return Scalar{ .f = a.f - b.f };
        }

        /// returns a * b
        pub inline fn mul(a: Scalar, b: Scalar) Scalar {
            return Scalar{ .f = a.f * b.f };
        }

        /// returns a / b
        pub inline fn div(a: Scalar, b: Scalar) Scalar {
            return Scalar{ .f = a.f / b.f };
        }

        /// returns -a
        pub inline fn neg(a: Scalar) Scalar {
            return Scalar{ .f = -a.f };
        }

        /// returns |a|
        pub inline fn abs(a: Scalar) Scalar {
            return Scalar{ .f = @fabs(a.f) };
        }

        /// returns |√a|
        pub inline fn sqrt(a: Scalar) Scalar {
            assert(@reduce(.And, a.cmp(.gte, Scalar.zero)));
            return Scalar{ .f = @sqrt(a.f) };
        }

        /// return log2(a)
        pub inline fn log2(a: Scalar) Scalar {
            assert(@reduce(.And, a.cmp(.gt, Scalar.zero)));
            return Scalar{ .f = @log2(a.f) };
        }

        /// returns truth value of: a r b
        pub inline fn cmp(a: Scalar, r: CompareOperator, b: Scalar) @Vector(len, bool) {
            return switch (r) {
                .eq => a.f == b.f,
                .lt => a.f < b.f,
                .lte => a.f <= b.f,
                .gt => a.f > b.f,
                .gte => a.f >= b.f,
                .neq => a.f != b.f,
            };
        }

        pub inline fn cmpAll(a: Scalar, r: CompareOperator, b: Scalar) bool {
            return @reduce(.And, a.cmp(r, b));
        }

        ///returns the minimum of a and b
        pub inline fn min(a: Scalar, b: Scalar) Scalar {
            return Scalar{ .f = @min(a.f, b.f) };
        }

        ///returns the maximum of a and b
        pub inline fn max(a: Scalar, b: Scalar) Scalar {
            return Scalar{ .f = @max(a.f, b.f) };
        }
    };
}

test "creation" {
    const fTypes = [_]type{ f16, f32, f64, f80, f128 };
    inline for (fTypes) |f| {
        const F = Float(f);

        try testing.expectEqual(@as(f, 0), F.zero.f[0]);
        try testing.expectEqual(@as(f, 1), F.eye.f[0]);

        try testing.expectEqual(@as(f, -3.14), F.from(-314, 100).f[0]);
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
        try testing.expectEqual(@log2(b.f), b.log2().f);

        try testing.expectEqual(@min(a.f, b.f), a.min(b).f);
        try testing.expectEqual(@max(a.f, b.f), a.max(b).f);
    }
}

test "comparisons" {
    const fTypes = [_]type{ f16, f32, f64, f80, f128 };
    inline for (fTypes) |f| {
        const F = Float(f);
        const a = F.from(-314, 100);
        const b = F.from(527, 100);
        try testing.expect(a.cmpAll(.lt, b));
        try testing.expect(a.cmpAll(.lte, b));
        try testing.expect(!a.cmpAll(.eq, b));
        try testing.expect(!a.cmpAll(.gte, b));
        try testing.expect(!a.cmpAll(.gt, b));
        try testing.expect(a.cmpAll(.neq, b));

        try testing.expect(!b.cmpAll(.lt, a));
        try testing.expect(!b.cmpAll(.lte, a));
        try testing.expect(!b.cmpAll(.eq, a));
        try testing.expect(b.cmpAll(.gte, a));
        try testing.expect(b.cmpAll(.gt, a));
        try testing.expect(b.cmpAll(.neq, a));

        try testing.expect(!a.cmpAll(.lt, a));
        try testing.expect(a.cmpAll(.lte, a));
        try testing.expect(a.cmpAll(.eq, a));
        try testing.expect(a.cmpAll(.gte, a));
        try testing.expect(!a.cmpAll(.gt, a));
        try testing.expect(!a.cmpAll(.neq, a));
    }
}
