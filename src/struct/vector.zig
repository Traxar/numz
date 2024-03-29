const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;
const Allocator = std.mem.Allocator;

/// struct for Vector based computations
pub fn VectorType(comptime Element: type) type {
    assert(Element.SIMDsize == 1);
    const SIMDElement = Element.SIMDType(null);
    const SIMDsize = SIMDElement.SIMDsize;
    return struct {
        const Vector = @This();
        val: []SIMDElement,
        len: usize,

        /// deinitialize vector
        pub fn deinit(a: Vector, allocator: Allocator) void {
            allocator.free(a.val);
        }

        ///allocate vector with undefined values
        pub fn init(n: usize, allocator: Allocator) !Vector {
            const n_SIMD = 1 + @divFloor(n - 1, SIMDsize); //divCeil
            const res = Vector{
                .val = try allocator.alloc(SIMDElement, n_SIMD),
                .len = n,
            };
            res.SIMDsetTail(SIMDElement.zero);
            return res;
        }

        ///allocate vector with undefined values and same dimensions as a
        pub fn like(a: Vector, allocator: Allocator) !Vector {
            return init(a.len, allocator);
        }

        /// res <- a
        pub fn copy(a: Vector, res: Vector) void {
            assert(a.len == res.len);
            if (a.val.ptr != res.val.ptr) {
                @memcpy(res.val, a.val);
            }
        }

        /// set all elements of a to b
        pub fn fill(res: Vector, a: Element) void {
            const a_SIMD = SIMDElement.SIMDsplat(a);
            for (0..res.val.len) |i| {
                res.val[i] = a_SIMD;
            }
            res.SIMDsetTail(SIMDElement.zero);
        }

        /// set tail elements of a to b
        pub fn SIMDsetTail(a: Vector, b: SIMDElement) void {
            const i = a.val.len - 1;
            const tail_size = a.len - i * SIMDsize;
            if (tail_size == 0) return;
            const uSIMDsize = std.meta.Int(.unsigned, SIMDsize);
            //bitcast from int to @vector(bool) reverses the order
            const pred: @Vector(SIMDsize, bool) = @bitCast((@as(uSIMDsize, 1) << @truncate(tail_size)) - 1);

            a.val[i] = a.val[i].SIMDselect(b, pred);
        }

        /// return element at index i
        pub fn at(a: Vector, i: usize) Element {
            assert(i < a.len);
            const i_SIMD = @divFloor(i, SIMDsize);
            const a_SIMD = a.val[i_SIMD];
            const i_sub = i - i_SIMD * SIMDsize;
            return a_SIMD.SIMDat(i_sub);
        }

        /// set element at index i to b
        pub fn set(a: Vector, i: usize, b: Element) void {
            assert(i < a.len);
            const i_SIMD = @divFloor(i, SIMDsize);
            const i_sub = i - i_SIMD * SIMDsize;
            a.val[i_SIMD].SIMDset(i_sub, b);
        }

        /// res <- a + b
        pub fn add(a: Vector, b: Vector, res: Vector) void {
            assert(res.len == a.len);
            assert(a.len == b.len);
            for (0..res.val.len) |i_SIMD| {
                res.val[i_SIMD] = a.val[i_SIMD].add(b.val[i_SIMD]);
            }
        }

        /// res <- a - b
        pub fn sub(a: Vector, b: Vector, res: Vector) void {
            assert(res.len == a.len);
            assert(a.len == b.len);
            for (0..res.val.len) |i_SIMD| {
                res.val[i_SIMD] = a.val[i_SIMD].sub(b.val[i_SIMD]);
            }
        }

        /// res <- a * b
        pub fn mulE(a: Vector, b: Element, res: Vector) void {
            assert(res.len == a.len);
            const c = SIMDElement.SIMDsplat(b);
            for (0..res.val.len) |i_SIMD| {
                res.val[i_SIMD] = a.val[i_SIMD].mul(c);
            }
        }

        /// res <- a + b * c
        pub fn mulEAdd(a: Vector, b: Element, c: Vector, res: Vector) void {
            assert(res.len == a.len);
            const b_ = SIMDElement.SIMDsplat(b);
            for (0..res.val.len) |i_SIMD| {
                res.val[i_SIMD] = a.val[i_SIMD].add(b_.mul(c.val[i_SIMD]));
            }
        }

        /// res <- a * b
        pub fn mul(a: Vector, b: Vector, res: Vector) void {
            assert(res.len == a.len);
            assert(a.len == b.len);
            for (0..res.val.len) |i_SIMD| {
                res.val[i_SIMD] = a.val[i_SIMD].mul(b.val[i_SIMD]);
            }
        }

        /// res <- a + b * c
        pub fn mulAdd(a: Vector, b: Vector, c: Vector, res: Vector) void {
            assert(res.len == a.len);
            assert(res.len == b.len);
            assert(res.len == c.len);
            for (0..res.val.len) |i_SIMD| {
                res.val[i_SIMD] = a.val[i_SIMD].add(b.val[i_SIMD].mul(c.val[i_SIMD]));
            }
        }

        /// res <- a / b
        pub fn divE(a: Vector, b: Element, res: Vector) void {
            assert(res.len == a.len);
            const c = SIMDElement.SIMDsplat(b);
            for (0..res.val.len) |i_SIMD| {
                res.val[i_SIMD] = a.val[i_SIMD].div(c);
            }
        }

        /// res <- a / b
        pub fn div(a: Vector, b: Vector, res: Vector) void {
            assert(res.len == a.len);
            assert(a.len == b.len);
            for (0..res.val.len) |i_SIMD| {
                res.val[i_SIMD] = a.val[i_SIMD].div(b.val[i_SIMD]);
            }
        }

        /// return sum of elements
        pub fn sum(a: Vector) Element {
            var res = SIMDElement.zero;
            for (0..a.val.len) |i_SIMD| {
                res = res.add(a.val[i_SIMD]);
            }
            return res.SIMDreduce(.Add);
        }

        ///return a . b
        pub fn dot(a: Vector, b: Vector) Element {
            var res = SIMDElement.zero;
            for (0..a.val.len) |i_SIMD| {
                res = res.add(a.val[i_SIMD].mul(b.val[i_SIMD]));
            }
            return res.SIMDreduce(.Add);
        }

        ///return the euclidean norm
        pub fn norm(a: Vector) Element {
            return a.dot(a).sqrt();
        }

        ///return the smallest element
        pub fn min(a: Vector) Element {
            a.SIMDsetTail(SIMDElement.SIMDsplat(a.at(0)));
            var res = a.val[0];
            for (1..a.val.len) |i_SIMD| {
                res = res.min(a.val[i_SIMD]);
            }
            a.SIMDsetTail(SIMDElement.zero);
            return res.SIMDreduce(.Min);
        }

        ///return the largest element
        pub fn max(a: Vector) Element {
            a.SIMDsetTail(SIMDElement.SIMDsplat(a.at(0)));
            var res = a.val[0];
            for (1..a.val.len) |i_SIMD| {
                res = res.max(a.val[i_SIMD]);
            }
            a.SIMDsetTail(SIMDElement.zero);
            return res.SIMDreduce(.Max);
        }

        pub fn debugprint(self: Vector) void {
            std.debug.print("\n", .{});
            for (0..self.len) |i| {
                std.debug.print("{}\n", .{self.at(i).f});
            }
        }
    };
}

test "vector creation" {
    const ally = std.testing.allocator;
    const n = 101;
    const F = @import("field.zig").Float(f32);
    try testing.expect(F.SIMDsize == 1);
    const V = VectorType(F);

    const a = F.from(-314, 100);
    var v = try V.init(n, ally);
    defer v.deinit(ally);
    v.fill(a);

    try testing.expectEqual(@as(usize, n), v.len);

    for (0..v.len) |i| {
        try testing.expect(a.cmp(.eq, v.at(i)));
    }
}

test "vector operators" {
    const ally = std.testing.allocator;
    const n = 3;
    const F = @import("field.zig").Float(f32);
    const V = VectorType(F);

    const a_ = F.from(-314, 100);
    const b_ = F.from(527, 100);
    var a = try V.init(n, ally);
    defer a.deinit(ally);
    a.fill(a_);
    a.set(1, F.zero);
    try testing.expect(a.at(0).cmp(.eq, a_));
    try testing.expect(a.at(1).cmp(.eq, F.zero));
    try testing.expect(a.at(2).cmp(.eq, a_));

    var b = try V.init(n, ally);
    defer b.deinit(ally);
    b.fill(b_);
    b.set(0, F.eye);
    try testing.expect(b.at(0).cmp(.eq, F.eye));
    try testing.expect(b.at(1).cmp(.eq, b_));
    try testing.expect(b.at(2).cmp(.eq, b_));

    const c = try V.init(n, ally);
    a.add(b, c);
    defer c.deinit(ally);
    try testing.expect(c.at(0).cmp(.eq, a_.add(F.eye)));
    try testing.expect(c.at(1).cmp(.eq, b_));
    try testing.expect(c.at(2).cmp(.eq, a_.add(b_)));

    a.sub(b, c);
    try testing.expect(c.at(0).cmp(.eq, a_.sub(F.eye)));
    try testing.expect(c.at(1).cmp(.eq, b_.neg()));
    try testing.expect(c.at(2).cmp(.eq, a_.sub(b_)));

    a.mul(b, c);
    try testing.expect(c.at(0).cmp(.eq, a_));
    try testing.expect(c.at(1).cmp(.eq, F.zero));
    try testing.expect(c.at(2).cmp(.eq, a_.mul(b_)));

    a.div(b, c);
    try testing.expect(c.at(0).cmp(.eq, a_));
    try testing.expect(c.at(1).cmp(.eq, F.zero));
    try testing.expect(c.at(2).cmp(.eq, a_.div(b_)));

    try testing.expect(a.dot(b).cmp(.eq, a_.add(a_.mul(b_))));

    try testing.expect(b.sum().cmp(.eq, b_.add(b_).add(F.eye)));

    const a_2 = a_.mul(a_);
    try testing.expect(a.norm().cmp(.eq, a_2.add(a_2).sqrt()));

    try testing.expect(a.min().cmp(.eq, a_.min(F.zero)));
    try testing.expect(a.max().cmp(.eq, a_.max(F.zero)));
    try testing.expect(b.min().cmp(.eq, b_.min(F.eye)));
    try testing.expect(b.max().cmp(.eq, b_.max(F.eye)));
}
