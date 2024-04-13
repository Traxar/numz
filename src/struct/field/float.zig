const std = @import("std");
const testing = std.testing;
const assert = std.debug.assert;

/// Scalar data structure for computation.
/// Choose `float` from `{f16, f32, f64, f80, f128}`.
pub fn FloatType(comptime Float: type) type {
    return SimdFloatType(1, Float);
}

/// Scalar data structure for SIMD computation.
/// This is ment for internal use, for scalar computations in userspace use `FloatType()` instead.
/// Choose `float` from `{f16, f32, f64, f80, f128}`.
/// `size` determines the SIMD size. if `size` is set to `null` the SIMD size will be choosen based on a heuristic
fn SimdFloatType(comptime size: ?usize, comptime Float: type) type {
    comptime {
        if (@typeInfo(Float) != .Float) @compileError("provided Type must be a float");
        if (size == 0) @compileError("SIMD size = 0 is not allowed");
    }
    return struct {
        const Element = @This();
        pub const simd_size = size orelse (std.simd.suggestVectorLength(Float) orelse 1);

        f: @Vector(simd_size, Float),

        /// return a version of `@This()` with SIMD size `sz`
        pub fn SimdType(comptime sz: ?usize) type {
            return SimdFloatType(sz, Float);
        }

        /// return element isomorph to p/q
        pub fn from(p: isize, q: usize) Element {
            comptime {
                if (simd_size != 1) @compileError("can not create SimdElement from fraction");
            }
            return Element{ .f = @splat(@as(Float, @floatFromInt(p)) / @as(Float, @floatFromInt(q))) };
        }

        /// the 0 element
        /// this is the default element in a sparse struct
        pub const zero = Element{ .f = @splat(0) };

        /// the 1 element
        pub const one = Element{ .f = @splat(1) };

        /// result type of comparisons
        const Boolean = if (simd_size == 1) bool else @Vector(simd_size, bool);

        pub inline fn all(boolean: Boolean) bool {
            return if (simd_size == 1) boolean else @reduce(.And, boolean);
        }

        /// a == b
        pub inline fn eq(a: Element, b: Element) Boolean {
            return @bitCast(a.f == b.f);
        }

        /// a < b
        pub inline fn lt(a: Element, b: Element) Boolean {
            return @bitCast(a.f < b.f);
        }

        /// a <= b
        pub inline fn lte(a: Element, b: Element) Boolean {
            return @bitCast(a.f <= b.f);
        }

        /// a > b
        pub inline fn gt(a: Element, b: Element) Boolean {
            return @bitCast(a.f > b.f);
        }

        /// a >= b
        pub inline fn gte(a: Element, b: Element) Boolean {
            return @bitCast(a.f >= b.f);
        }

        /// a != b
        pub inline fn neq(a: Element, b: Element) Boolean {
            return @bitCast(a.f != b.f);
        }

        /// comparators defined on the field
        pub const Comparator = enum {
            eq,
            lt,
            lte,
            gt,
            gte,
            neq,

            /// <- if (T) TypeOf(op) else op
            inline fn fT(comptime cmp: Comparator, comptime T: bool) if (T) type else cmp.fT(true) {
                const f_ = switch (cmp) {
                    .eq => eq,
                    .lt => lt,
                    .lte => lte,
                    .gt => gt,
                    .gte => gte,
                    .neq => neq,
                };
                return if (T) @TypeOf(f_) else f_;
            }

            /// <- function assigned to op
            pub inline fn f(comptime cmp: Comparator) cmp.fT(true) {
                return cmp.fT(false);
            }

            // check domain and range of comparators:
            comptime {
                for (@typeInfo(Comparator).Enum.fields) |comparator| {
                    const cmp: Comparator = @enumFromInt(comparator.value);
                    const info = @typeInfo(@TypeOf(cmp.f())).Fn;
                    if (info.return_type.? != Boolean) @compileError("comparators must return 'Boolean'");
                    for (info.params) |param| {
                        if (param.type != Element) @compileError("comparators may only take 'Element's as arguments");
                    }
                }
            }
        };

        /// a
        pub inline fn id(a: Element) Element {
            return a;
        }

        /// a + b
        pub inline fn add(a: Element, b: Element) Element {
            return .{ .f = a.f + b.f };
        }

        /// -a
        pub inline fn neg(a: Element) Element {
            return .{ .f = -a.f };
        }

        /// a - b
        pub inline fn sub(a: Element, b: Element) Element {
            return .{ .f = a.f - b.f };
        }

        /// a * b
        pub inline fn mul(a: Element, b: Element) Element {
            return .{ .f = a.f * b.f };
        }

        //TODO: checkout why compiler segfaults on inlining this
        /// a^-1
        pub fn inv(a: Element) !Element {
            if (!Element.all(a.neq(Element.zero))) return error.DivisionByZero;
            return .{ .f = @as(@TypeOf(a.f), @splat(1.0)) / a.f };
        }

        /// a / b
        pub fn div(a: Element, b: Element) !Element {
            if (!Element.all(b.neq(Element.zero))) return error.DivisionByZero;
            return .{ .f = a.f / b.f };
        }

        /// min(a,b)
        pub inline fn min(a: Element, b: Element) Element {
            return .{ .f = @min(a.f, b.f) };
        }

        /// max(a,b)
        pub inline fn max(a: Element, b: Element) Element {
            return .{ .f = @max(a.f, b.f) };
        }

        /// |a|
        pub inline fn abs(a: Element) Element {
            return .{ .f = @abs(a.f) };
        }

        /// √a
        pub fn sqrt(a: Element) !Element {
            if (!Element.all(a.gte(Element.zero))) return error.SqrtOfNegative;
            return .{ .f = @sqrt(a.f) };
        }

        /// operators defined on the field
        pub const Operator = enum {
            id,
            add,
            neg,
            sub,
            mul,
            inv,
            div,
            min,
            max,
            abs,
            sqrt,

            /// <- if (T) TypeOf(op) else op
            inline fn fT(comptime op: Operator, comptime T: bool) if (T) type else op.fT(true) {
                const f_ = switch (op) {
                    .id => id,
                    .add => add,
                    .neg => neg,
                    .sub => sub,
                    .mul => mul,
                    .inv => inv,
                    .div => div,
                    .min => min,
                    .max => max,
                    .abs => abs,
                    .sqrt => sqrt,
                };
                return if (T) @TypeOf(f_) else f_;
            }

            /// <- function assigned to op
            pub inline fn f(comptime op: Operator) op.fT(true) {
                return op.fT(false);
            }

            /// op defined for all inputs
            pub inline fn ErrorSet(comptime op: Operator) ?type {
                comptime {
                    const info = @typeInfo(@TypeOf(op.f())).Fn;
                    if (info.return_type.? == Element) return null;
                    const ret_info = @typeInfo(info.return_type.?).ErrorUnion;
                    return ret_info.error_set;
                }
            }

            // check domain and range of operators:
            comptime {
                for (@typeInfo(Operator).Enum.fields) |operator| {
                    const op: Operator = @enumFromInt(operator.value);
                    const info = @typeInfo(@TypeOf(op.f())).Fn;
                    if (info.return_type.? != Element) {
                        const ret_info = @typeInfo(info.return_type.?).ErrorUnion;
                        if (ret_info.payload != Element) {
                            @compileError("operators must return 'Element' or '!Element'");
                        }
                    }
                    for (info.params) |param| {
                        if (param.type != Element) @compileError("operators may only take 'Element's as arguments");
                    }
                }
            }

            /// arg[i]: bool
            /// <- args => (op(args) == 0)
            /// This lets the sparse matrix know what happens to the majority of its elements
            pub inline fn zero(comptime op: Operator, args: [@typeInfo(op.fT(true)).Fn.params.len]bool) bool {
                return switch (op) {
                    .id => args[0],
                    .add => args[0] and args[1],
                    .neg => args[0],
                    .sub => args[0] and args[1],
                    .mul => args[0] or args[1],
                    .inv => false,
                    .div => args[0],
                    .min => args[0] and args[1],
                    .max => args[0] and args[1],
                    .abs => args[0],
                    .sqrt => args[0],
                };
            }
        };

        /// splat Element a of size 1 to simd_size
        pub inline fn simdSplat(a: SimdType(1)) Element {
            return Element{ .f = @splat(@bitCast(a.f)) };
        }

        /// reduce to size 1 using op
        pub inline fn simdReduce(a: Element, comptime op: Operator) SimdType(1) {
            return switch (op) {
                .add => .{ .f = @bitCast(@reduce(.Add, a.f)) },
                .mul => .{ .f = @bitCast(@reduce(.Mul, a.f)) },
                .min => .{ .f = @bitCast(@reduce(.Min, a.f)) },
                .max => .{ .f = @bitCast(@reduce(.Max, a.f)) },
                else => @compileError("Operator is not commutative"),
            };
        }

        /// returns the element at position i
        pub inline fn simdAt(a: Element, i: usize) SimdType(1) {
            if (simd_size == 1) {
                return a;
            } else {
                return .{ .f = @bitCast(a.f[i]) };
            }
        }

        /// sets the element at position i to b
        pub inline fn simdSet(a: *Element, i: usize, b: SimdType(1)) void {
            a.f[i] = @bitCast(b.f);
        }

        /// returns simd element with entries choosen according to b
        /// if b[i] is true
        ///     choose a[i]
        /// else
        ///     choose c[i]
        pub inline fn simdSelect(a: Element, b: @Vector(simd_size, bool), c: Element) Element {
            return Element{ .f = @select(Float, b, a.f, c.f) };
        }
    };
}

test "float 1" {
    const floatTypes = .{ f16, f32, f64, f80, f128 };
    inline for (floatTypes) |fx| {
        const F = FloatType(fx);
        try testing.expect(F.zero.eq(F.from(0, 1)));
        try testing.expect(F.one.eq(F.from(1, 1)));
        try testing.expect(F.zero.neq(F.one));

        const two = F.one.add(F.one);
        // 2 - 1 = 1
        try testing.expect(two.sub(F.one).eq(F.one));
        // (-1) * (-1) = 1
        try testing.expect(F.one.neg().mul(F.one.neg()).eq(F.one));
        // 2^-1 + (1 / 2) = 1
        try testing.expect((try two.inv()).add(try F.one.div(two)).eq(F.one));
        // √(2 * 2) = |-2|
        try testing.expect((try two.mul(two).sqrt()).eq(two.neg().abs()));
    }
}

test "float simd" {
    const simd_sizes = .{ 1, null };
    inline for (simd_sizes) |simd| {
        const float_types = .{ f16, f32, f64, f80 }; //TODO: add f128 when fixed upstream
        inline for (float_types) |fx| {
            const F = FloatType(fx);
            const SimdF = F.SimdType(simd);

            const a = SimdF.simdSplat(F.from(-1, 1));
            var b = SimdF.simdSplat(F.from(3, 1));
            b.simdSet(0, F.zero);
            const c = SimdF.simdSplat(F.from(2, 1));
            var pred: @Vector(SimdF.simd_size, bool) = @splat(false);
            pred[0] = true;

            try testing.expect(SimdF.all(a.add(b).eq(a.simdSelect(pred, c))));
            try testing.expect(a.add(b).simdReduce(.min).eq(F.from(-1, 1)));
            if (SimdF.simd_size > 1)
                try testing.expect(a.add(b).simdReduce(.max).eq(F.from(2, 1)))
            else
                try testing.expect(a.add(b).simdReduce(.max).eq(F.from(-1, 1)));
        }
    }
}

test "float zero operations" {
    const F = FloatType(f32);
    try testing.expect(F.Operator.add.zero(.{ true, true }));
    try testing.expect(!F.Operator.add.zero(.{ true, false }));
    try testing.expect(F.Operator.mul.zero(.{ true, false }));
    try testing.expect(!F.Operator.mul.zero(.{ false, false }));
    try testing.expect(F.Operator.sqrt.zero(.{true}));
    try testing.expect(!F.Operator.abs.zero(.{false}));
}

test "float enum calls" {
    const floatTypes = .{ f16, f32, f64, f80, f128 };
    inline for (floatTypes) |fx| {
        const F = FloatType(fx);
        try testing.expect(@call(
            .auto,
            F.Comparator.eq.f(),
            .{
                @call(
                    .auto,
                    F.Operator.add.f(),
                    .{ F.one, F.one },
                ),
                F.from(2, 1),
            },
        ));
    }
}
