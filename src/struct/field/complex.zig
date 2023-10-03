const std = @import("std");
const testing = std.testing;
const assert = std.debug.assert;
const CompareOperator = std.math.CompareOperator;

// TODO sqrt
// TODO log

/// Scalar data structure for complex computation.
pub fn Complex(comptime Scalar: type) type {
    return SIMDComplex(1, Scalar);
}

/// Scalar data structure for SIMD computation.
/// This is ment for internal use, for scalar computations in userspace use `Complex()` instead.
/// `size` determines the SIMD size. if `size` is set to `null` the SIMD size will be choosen based on a heuristic
pub fn SIMDComplex(comptime size: ?usize, comptime Scalar: type) type {
    assert(Scalar.SIMDsize == 1);
    assert(size != 0);
    return struct {
        const Element = @This();
        pub const SIMDsize: usize = size orelse Scalar.SIMDType(null).SIMDsize;
        const SIMDScalar = Scalar.SIMDType(SIMDsize);
        re: SIMDScalar, // real part
        im: SIMDScalar, // imaginary part

        /// return a SIMD version of `@This()` with SIMD size `size_`
        pub fn SIMDType(comptime size_: ?usize) type {
            return SIMDComplex(size_, Scalar);
        }

        /// the 0 element
        pub const zero = Element{ .re = SIMDScalar.zero, .im = SIMDScalar.zero };

        /// the 1 element
        pub const eye = Element{ .re = SIMDScalar.eye, .im = SIMDScalar.zero };

        /// return scalar with undefined value
        pub inline fn init() Element {
            return Element{ .re = undefined, .im = undefined };
        }

        /// return element isomorph to a + i*b
        pub fn from(a: SIMDScalar, b: SIMDScalar) Element {
            return Element{ .re = a, .im = b };
        }

        /// returns a + b
        pub inline fn add(a: Element, b: Element) Element {
            return Element{ .re = a.re.add(b.re), .im = a.im.add(b.im) };
        }

        /// returns a - b
        pub inline fn sub(a: Element, b: Element) Element {
            return Element{ .re = a.re.sub(b.re), .im = a.im.sub(b.im) };
        }

        /// returns a * b
        pub inline fn mul(a: Element, b: Element) Element {
            return Element{ .re = a.re.mul(b.re).sub(a.im.mul(b.im)), .im = a.re.mul(b.im).add(a.im.mul(b.re)) };
        }

        /// returns a / b
        pub inline fn div(a: Element, b: Element) Element {
            const q = b.re.mul(b.re).add(b.im.mul(b.im));
            return Element{
                .re = a.re.mul(b.re).add(a.im.mul(b.im)).div(q),
                .im = a.im.mul(b.re).sub(a.re.mul(b.im)).div(q),
            };
        }

        /// returns -a
        pub inline fn neg(a: Element) Element {
            return Element{
                .re = a.re.neg(),
                .im = a.im.neg(),
            };
        }

        /// return a^-1
        pub inline fn inv(a: Element) Element {
            const q = a.re.mul(a.re).add(a.im.mul(a.im));
            return Element{
                .re = a.re.div(q),
                .im = a.im.neg().div(q),
            };
        }

        /// returns |a|
        pub inline fn abs(a: Element) Scalar {
            return a.re.mul(a.re).add(a.im.mul(a.im)).sqrt();
        }

        /// returns true iff a r b holds for all SIMD elements
        pub inline fn cmp(a: Element, r: CompareOperator, b: Element) bool {
            return @reduce(.And, a.SIMDcompare(r, b));
        }

        /// returns SIMD element with all entries set to a
        pub fn SIMDsplat(a: SIMDType(1)) Element {
            return Element{
                .re = SIMDScalar.SIMDsplat(a.re),
                .im = SIMDScalar.SIMDsplat(a.im),
            };
        }

        /// returns the element at position i
        pub fn SIMDat(a: Element, i: usize) SIMDType(1) {
            return SIMDType(1){
                .re = a.re.SIMDat(i),
                .im = a.im.SIMDat(i),
            };
        }

        /// sets the element at position i to b
        pub fn SIMDset(a: *Element, i: usize, b: SIMDType(1)) void {
            a.re.SIMDset(i, b.re);
            a.im.SIMDset(i, b.im);
        }

        /// return red(a)
        pub inline fn SIMDreduce(a: Element, comptime red: std.builtin.ReduceOp) SIMDType(1) {
            return switch (red) {
                .Add => SIMDType(1){
                    .re = a.re.SIMDreduce(red),
                    .im = a.re.SIMDreduce(red),
                },
                .Mul => unreachable, //TODO
                else => unreachable,
            };
        }

        /// returns SIMD element with entries choosen according to c
        /// if c[i] is true choose a[i]
        /// else choose b[i]
        pub fn SIMDselect(a: Element, b: Element, c: @Vector(SIMDsize, bool)) Element {
            return Element{
                .re = a.re.SIMDselect(b.re, c),
                .im = a.im.SIMDselect(b.im, c),
            };
        }

        /// returns truth value of a r b
        inline fn SIMDcompare(a: Element, r: CompareOperator, b: Element) @Vector(SIMDsize, bool) {
            return switch (r) {
                .eq => a.re.SIMDcompare(r, b.re) & a.im.SIMDcompare(r, b.im),
                .neq => a.re.SIMDcompare(r, b.re) | a.im.SIMDcompare(r, b.im),
                else => unreachable,
            };
        }
    };
}

test "complex creation" {
    const F = @import("float.zig").Float(f32);
    const C = Complex(F);

    const zero = C.zero;
    try testing.expectEqual(F.zero, zero.re);
    try testing.expectEqual(F.zero, zero.im);

    const one = C.eye;
    try testing.expectEqual(F.eye, one.re);
    try testing.expectEqual(F.zero, one.im);

    const a = F.from(-314, 100);
    const c = C.from(a, F.eye);

    try testing.expectEqual(a, c.re);
    try testing.expectEqual(F.eye, c.im);
}

test "complex operators" {
    const F = @import("float.zig").Float(f32);
    const C = Complex(F);

    const a = C.eye;
    const b = C.from(F.eye.neg(), F.eye);

    try testing.expectEqual(F.zero, a.add(b).re);
    try testing.expectEqual(F.eye, a.add(b).im);

    try testing.expectEqual(F.eye.add(F.eye), a.sub(b).re);
    try testing.expectEqual(F.eye.neg(), a.sub(b).im);

    try testing.expectEqual(F.eye.add(F.eye).neg(), b.sub(a).re);
    try testing.expectEqual(F.eye, b.sub(a).im);

    try testing.expectEqual(b, a.mul(b));

    try testing.expectEqual(F.from(-1, 2), a.div(b).re);
    try testing.expectEqual(F.from(-1, 2), a.div(b).re);

    try testing.expectEqual(b, b.div(a));

    try testing.expectEqual(F.eye.add(F.eye).sqrt(), b.abs());
}
