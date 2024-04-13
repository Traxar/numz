const std = @import("std");
const testing = std.testing;
const assert = std.debug.assert;

/// Complex data structure for computation.
/// `Scalar` must already be a field.
pub fn ComplexType(comptime Scalar: type) type {
    return SimdComplexType(1, Scalar);
}

/// Complex data structure for SIMD computation.
/// This is ment for internal use, for scalar computations in userspace use `ComplexType()` instead.
/// `size` determines the SIMD size. if `size` is set to `null` the SIMD size will be choosen based on a heuristic
fn SimdComplexType(comptime size: ?usize, comptime Scalar: type) type {
    comptime {
        if (size == 0) @compileError("SIMD size = 0 is not allowed");
    }
    return struct {
        const Element = @This();

        pub const simd_size: usize = size orelse Scalar.SimdType(null).simd_size;
        const SimdScalar = Scalar.SimdType(simd_size);
        re: SimdScalar, // real part
        im: SimdScalar, // imaginary part

        /// return a version of `@This()` with SIMD size `sz`
        pub fn SimdType(comptime sz: ?usize) type {
            return SimdComplexType(sz, Scalar);
        }

        /// return element isomorph to a + ib
        pub fn from(a: SimdScalar, b: SimdScalar) Element {
            return Element{ .re = a, .im = b };
        }

        /// the 0 element
        /// this is the default element in a sparse struct
        pub const zero = Element{ .re = SimdScalar.zero, .im = SimdScalar.zero };

        /// the 1 element
        pub const one = Element{ .re = SimdScalar.one, .im = SimdScalar.zero };

        /// imaginary unit
        pub const i = Element{ .re = SimdScalar.zero, .im = SimdScalar.one };

        /// result type of comparisons
        const Boolean = if (simd_size == 1) bool else @Vector(simd_size, bool);

        pub fn all(boolean: Boolean) bool {
            return if (simd_size == 1) boolean else @reduce(.And, boolean);
        }

        /// a == b
        pub fn eq(a: Element, b: Element) Boolean {
            //FIXME: https://github.com/ziglang/zig/issues/14306
            const re: @Vector(simd_size, u1) = @bitCast(a.re.eq(b.re));
            const im: @Vector(simd_size, u1) = @bitCast(a.im.eq(b.im));
            return @bitCast(re & im);
        }

        /// a != b
        pub fn neq(a: Element, b: Element) Boolean {
            //FIXME: https://github.com/ziglang/zig/issues/14306
            const re: @Vector(simd_size, u1) = @bitCast(a.re.neq(b.re));
            const im: @Vector(simd_size, u1) = @bitCast(a.im.neq(b.im));
            return @bitCast(re | im);
        }

        /// comparators defined on the field
        pub const Comparator = enum {
            eq,
            neq,

            /// <- if (T) TypeOf(op) else op
            inline fn fT(comptime cmp: Comparator, comptime T: bool) if (T) type else cmp.fT(true) {
                const f_ = switch (cmp) {
                    .eq => eq,
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
        pub fn id(a: Element) Element {
            return a;
        }

        /// a'
        pub fn conj(a: Element) Element {
            return .{
                .re = a.re,
                .im = a.im.neg(),
            };
        }

        /// a + b
        pub fn add(a: Element, b: Element) Element {
            return .{
                .re = a.re.add(b.re),
                .im = a.im.add(b.im),
            };
        }

        /// -a
        pub fn neg(a: Element) Element {
            return .{
                .re = a.re.neg(),
                .im = a.im.neg(),
            };
        }

        /// a - b
        pub fn sub(a: Element, b: Element) Element {
            return .{
                .re = a.re.sub(b.re),
                .im = a.im.sub(b.im),
            };
        }

        /// a * b
        pub fn mul(a: Element, b: Element) Element {
            return .{
                .re = a.re.mul(b.re).sub(a.im.mul(b.im)),
                .im = a.re.mul(b.im).add(a.im.mul(b.re)),
            };
        }

        /// a^-1
        pub fn inv(a: Element) !Element {
            const q = a.re.mul(a.re).add(a.im.mul(a.im));
            return .{
                .re = try a.re.div(q),
                .im = a.im.neg().div(q) catch unreachable,
            };
        }

        /// a / b
        pub fn div(a: Element, b: Element) !Element {
            const q = b.re.mul(b.re).add(b.im.mul(b.im));
            return .{
                .re = try a.re.mul(b.re).add(a.im.mul(b.im)).div(q),
                .im = a.im.mul(b.re).sub(a.re.mul(b.im)).div(q) catch unreachable,
            };
        }

        /// |a|
        pub fn abs(a: Element) Element {
            return .{
                .re = a.re.mul(a.re).add(a.im.mul(a.im)).sqrt() catch unreachable,
                .im = SimdScalar.zero,
            };
        }

        /// √a
        pub fn sqrt(a: Element) Element {
            const r = a.abs();
            var b = a.add(r);
            b = Element.i.simdSelect(@bitCast(b.eq(Element.zero)), b);
            const s = r.re.sqrt() catch unreachable;
            const q = s.div(b.abs().re) catch unreachable;
            return .{
                .re = b.re.mul(q),
                .im = b.im.mul(q),
            };
        }

        /// operators defined on the field
        pub const Operator = enum {
            id,
            conj,
            add,
            neg,
            sub,
            mul,
            inv,
            div,
            abs,
            sqrt,

            /// <- if (T) TypeOf(op) else op
            inline fn fT(comptime op: Operator, comptime T: bool) if (T) type else op.fT(true) {
                const f_ = switch (op) {
                    .id => id,
                    .conj => conj,
                    .add => add,
                    .neg => neg,
                    .sub => sub,
                    .mul => mul,
                    .inv => inv,
                    .div => div,
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

            // arg[i]: bool
            // <- args => (op(args) == 0)
            pub fn zero(comptime op: Operator, args: anytype) bool {
                comptime {
                    const info = @typeInfo(op.fT(true)).Fn;
                    if (args.len != info.params.len) @compileError("number of args must match operator");
                    for (0..args.len) |j| {
                        if (@TypeOf(args[j]) != bool) @compileError("all args must be boolean");
                    }
                }
                return switch (op) {
                    .id => args[0],
                    .conj => args[0],
                    .add => args[0] and args[1],
                    .neg => args[0],
                    .sub => args[0] and args[1],
                    .mul => args[0] or args[1],
                    .inv => false,
                    .div => args[0],
                    .abs => args[0],
                    .sqrt => args[0],
                };
            }
        };

        /// splat Element a of size 1 to simd_size
        pub fn simdSplat(a: SimdType(1)) Element {
            return Element{
                .re = SimdScalar.simdSplat(a.re),
                .im = SimdScalar.simdSplat(a.im),
            };
        }

        /// reduce to size 1 using op
        pub fn simdReduce(a: Element, comptime op: Operator) SimdType(1) {
            switch (op) {
                .add => return .{
                    .re = a.re.simdReduce(.add),
                    .im = a.im.simdReduce(.add),
                },
                .mul => {
                    var res = a.simdAt(0);
                    inline for (1..simd_size) |j| {
                        res = res.mul(a.simdAt(j));
                    }
                    return res;
                },
                else => @compileError("Operator is not commutative"),
            }
        }

        /// returns the element at position j
        pub fn simdAt(a: Element, j: usize) SimdType(1) {
            if (simd_size == 1) {
                return a;
            } else {
                return .{
                    .re = a.re.simdAt(j),
                    .im = a.im.simdAt(j),
                };
            }
        }

        /// sets the element at position j to b
        pub fn simdSet(a: *Element, j: usize, b: SimdType(1)) void {
            a.re.simdSet(j, b.re);
            a.im.simdSet(j, b.im);
        }

        /// returns simd element with entries choosen according to c
        /// if b[i] is true
        ///     choose a[i]
        /// else
        ///     choose c[i]
        pub fn simdSelect(a: Element, b: @Vector(simd_size, bool), c: Element) Element {
            return .{
                .re = a.re.simdSelect(b, c.re),
                .im = a.im.simdSelect(b, c.im),
            };
        }
    };
}

test "complex" {
    const F = @import("float.zig").FloatType(f32);
    const C = ComplexType(F);

    try testing.expect(C.one.neg().sqrt().eq(C.i));
    try testing.expect(C.from(F.zero, F.from(2, 1)).sqrt().eq(C.one.add(C.i)));
}
