const std = @import("std");
const testing = std.testing;
const assert = std.debug.assert;
const CompareOperator = std.math.CompareOperator;
const Allocator = std.mem.Allocator;

/// Scalar data structure for computation.
/// Choose `float` from `{f16, f32, f64, f80, f128}`.
pub fn Float(comptime float: type) type {
    return SIMDFloat(1, float);
}

/// Scalar data structure for SIMD computation.
/// This is ment for internal use, for scalar computations in userspace use `Float()` instead.
/// Choose `float` from `{f16, f32, f64, f80, f128}`.
/// `size` determines the SIMD size. if `size` is set to `null` the SIMD size will be choosen based on a heuristic
pub fn SIMDFloat(comptime size: ?usize, comptime float: type) type {
    assert(size != 0);
    return struct {
        const Scalar = @This();
        pub const SIMDsize: usize = size orelse (std.simd.suggestVectorSize(float) orelse 1);
        f: @Vector(SIMDsize, float),

        /// return a SIMD version of `@This()` with SIMD size `size_`
        pub fn SIMDType(comptime size_: ?usize) type {
            return SIMDFloat(size_, float);
        }

        /// the 0 element
        pub const zero = Scalar{ .f = @splat(0) };

        /// the 1 element
        pub const eye = Scalar{ .f = @splat(1) };

        /// return scalar with undefined value
        pub inline fn init() Scalar {
            return Scalar{ .f = undefined };
        }

        /// return element isomorph to p/q
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

        ///returns the minimum of a and b
        pub inline fn min(a: Scalar, b: Scalar) Scalar {
            return Scalar{ .f = @min(a.f, b.f) };
        }

        ///returns the maximum of a and b
        pub inline fn max(a: Scalar, b: Scalar) Scalar {
            return Scalar{ .f = @max(a.f, b.f) };
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
            assert(a.cmp(.gte, Scalar.zero));
            return Scalar{ .f = @sqrt(a.f) };
        }

        /// return log2(a)
        pub inline fn log2(a: Scalar) Scalar {
            assert(a.cmp(.gt, Scalar.zero));
            return Scalar{ .f = @log2(a.f) };
        }

        /// return red(a)
        pub inline fn SIMDred(a: Scalar, comptime red: std.builtin.ReduceOp) SIMDType(1) {
            return switch (red) {
                .Add, .Mul, .Min, .Max => SIMDType(1){ .f = [_]float{@reduce(red, a.f)} },
                else => unreachable,
            };
        }

        /// returns truth value of a r b
        fn SIMDcmp(a: Scalar, r: CompareOperator, b: Scalar) @Vector(SIMDsize, bool) {
            return switch (r) {
                .eq => a.f == b.f,
                .lt => a.f < b.f,
                .lte => a.f <= b.f,
                .gt => a.f > b.f,
                .gte => a.f >= b.f,
                .neq => a.f != b.f,
            };
        }

        /// returns true iff a r b holds for all SIMD elements
        pub inline fn cmp(a: Scalar, r: CompareOperator, b: Scalar) bool {
            return @reduce(.And, a.SIMDcmp(r, b));
        }
    };
}

test "float creation" {
    const fTypes = [_]type{ f16, f32, f64, f80, f128 };
    inline for (fTypes) |f| {
        const F = Float(f);

        try testing.expectEqual(@as(f, 0), F.zero.f[0]);
        try testing.expectEqual(@as(f, 1), F.eye.f[0]);

        try testing.expectEqual(@as(f, -3.14), F.from(-314, 100).f[0]);
    }
}

test "float operators" {
    const fTypes = [_]type{ f16, f32, f64, f80, f128 };
    inline for (fTypes) |f| {
        const F = Float(f);
        const a = F.from(-314, 100);
        const b = F.from(527, 100);
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

test "float comparisons" {
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
        try testing.expect(a.cmp(.neq, b));

        try testing.expect(!b.cmp(.lt, a));
        try testing.expect(!b.cmp(.lte, a));
        try testing.expect(!b.cmp(.eq, a));
        try testing.expect(b.cmp(.gte, a));
        try testing.expect(b.cmp(.gt, a));
        try testing.expect(b.cmp(.neq, a));

        try testing.expect(!a.cmp(.lt, a));
        try testing.expect(a.cmp(.lte, a));
        try testing.expect(a.cmp(.eq, a));
        try testing.expect(a.cmp(.gte, a));
        try testing.expect(!a.cmp(.gt, a));
        try testing.expect(!a.cmp(.neq, a));
    }
}

test "float SIMD" {
    const fTypes = [_]type{ f16, f32, f64, f80, f128 };
    inline for (fTypes) |f| {
        const F = SIMDFloat(null, f);
        const a = F.from(-314, 100);
        const b = F.from(527, 100);
        const c = F{ .f = a.f + b.f };
        try testing.expect(c.cmp(.eq, a.add(b)));
        _ = a.SIMDred(.Add);
    }
}
