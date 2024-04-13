const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;
const Allocator = std.mem.Allocator;

/// struct for Vector based computations
pub fn VectorType(comptime Element: type) type {
    assert(Element.simd_size == 1);
    const SimdElement = Element.SimdType(null);
    const simd_size = SimdElement.simd_size;
    return struct {
        const Vector = @This();
        val: []SimdElement,
        len: usize,

        /// deinitialize vector
        pub fn deinit(a: Vector, allocator: Allocator) void {
            allocator.free(a.val);
        }

        ///allocate vector with undefined values
        pub fn init(n: usize, allocator: Allocator) !Vector {
            const n_SIMD = 1 + @divFloor(n - 1, simd_size); //divCeil
            const res = Vector{
                .val = try allocator.alloc(SimdElement, n_SIMD),
                .len = n,
            };
            res.simdSetTail(SimdElement.zero);
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
            const a_SIMD = SimdElement.simdSplat(a);
            for (0..res.val.len) |i| {
                res.val[i] = a_SIMD;
            }
            res.simdSetTail(SimdElement.zero);
        }

        /// set tail elements of a to b
        pub fn simdSetTail(a: Vector, b: SimdElement) void {
            const i = a.val.len - 1;
            const tail_size = a.len - i * simd_size;
            if (tail_size == 0) return;
            const usimd_size = std.meta.Int(.unsigned, simd_size);
            //bitcast from int to @vector(bool) reverses the order
            const pred: @Vector(simd_size, bool) = @bitCast((@as(usimd_size, 1) << @truncate(tail_size)) - 1);

            a.val[i] = a.val[i].simdSelect(pred, b);
        }

        /// return element at index i
        pub fn at(a: Vector, i: usize) Element {
            assert(i < a.len);
            const i_SIMD = @divFloor(i, simd_size);
            const a_SIMD = a.val[i_SIMD];
            const i_sub = i - i_SIMD * simd_size;
            return a_SIMD.simdAt(i_sub);
        }

        /// set element at index i to b
        pub fn set(a: Vector, i: usize, b: Element) void {
            assert(i < a.len);
            const i_SIMD = @divFloor(i, simd_size);
            const i_sub = i - i_SIMD * simd_size;
            a.val[i_SIMD].simdSet(i_sub, b);
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
            const c = SimdElement.simdSplat(b);
            for (0..res.val.len) |i_SIMD| {
                res.val[i_SIMD] = a.val[i_SIMD].mul(c);
            }
        }

        /// res <- a + b * c
        pub fn mulEAdd(a: Vector, b: Element, c: Vector, res: Vector) void {
            assert(res.len == a.len);
            const b_ = SimdElement.simdSplat(b);
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
            const c = SimdElement.simdSplat(b);
            for (0..res.val.len) |i_SIMD| {
                res.val[i_SIMD] = a.val[i_SIMD].div(c);
            }
        }

        /// res <- a / b
        pub fn div(a: Vector, b: Vector, res: Vector) !void {
            assert(res.len == a.len);
            assert(a.len == b.len);
            b.simdSetTail(SimdElement.one);
            for (0..res.val.len) |i_SIMD| {
                res.val[i_SIMD] = try a.val[i_SIMD].div(b.val[i_SIMD]);
            }
            b.simdSetTail(SimdElement.zero);
        }

        /// return sum of elements
        pub fn sum(a: Vector) Element {
            var res = SimdElement.zero;
            for (0..a.val.len) |i_SIMD| {
                res = res.add(a.val[i_SIMD]);
            }
            return res.simdReduce(.add);
        }

        ///return a . b
        pub fn dot(a: Vector, b: Vector) Element {
            var res = SimdElement.zero;
            for (0..a.val.len) |i_SIMD| {
                res = res.add(a.val[i_SIMD].mul(b.val[i_SIMD]));
            }
            return res.simdReduce(.add);
        }

        ///return the euclidean norm
        pub fn norm(a: Vector) Element {
            return a.dot(a).sqrt() catch unreachable;
        }

        ///return the smallest element
        pub fn min(a: Vector) Element {
            a.simdSetTail(SimdElement.simdSplat(a.at(0)));
            var res = a.val[0];
            for (1..a.val.len) |i_SIMD| {
                res = res.min(a.val[i_SIMD]);
            }
            a.simdSetTail(SimdElement.zero);
            return res.simdReduce(.min);
        }

        ///return the largest element
        pub fn max(a: Vector) Element {
            a.simdSetTail(SimdElement.simdSplat(a.at(0)));
            var res = a.val[0];
            for (1..a.val.len) |i_SIMD| {
                res = res.max(a.val[i_SIMD]);
            }
            a.simdSetTail(SimdElement.zero);
            return res.simdReduce(.max);
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
    try testing.expect(F.simd_size == 1);
    const V = VectorType(F);

    const a = F.from(-314, 100);
    var v = try V.init(n, ally);
    defer v.deinit(ally);
    v.fill(a);

    try testing.expectEqual(@as(usize, n), v.len);

    for (0..v.len) |i| {
        try testing.expectEqual(a, v.at(i));
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
    try testing.expectEqual(a_, a.at(0));
    try testing.expectEqual(F.zero, a.at(1));
    try testing.expectEqual(a_, a.at(2));

    var b = try V.init(n, ally);
    defer b.deinit(ally);
    b.fill(b_);
    b.set(0, F.one);
    try testing.expectEqual(F.one, b.at(0));
    try testing.expectEqual(b_, b.at(1));
    try testing.expectEqual(b_, b.at(2));

    const c = try V.init(n, ally);
    a.add(b, c);
    defer c.deinit(ally);
    try testing.expectEqual(a_.add(F.one), c.at(0));
    try testing.expectEqual(b_, c.at(1));
    try testing.expectEqual(a_.add(b_), c.at(2));

    a.sub(b, c);
    try testing.expectEqual(a_.sub(F.one), c.at(0));
    try testing.expectEqual(b_.neg(), c.at(1));
    try testing.expectEqual(a_.sub(b_), c.at(2));

    a.mul(b, c);
    try testing.expectEqual(a_, c.at(0));
    try testing.expectEqual(F.zero, c.at(1));
    try testing.expectEqual(a_.mul(b_), c.at(2));

    try a.div(b, c);
    try testing.expectEqual(a_, c.at(0));
    try testing.expectEqual(F.zero, c.at(1));
    try testing.expectEqual(a_.div(b_), c.at(2));

    try testing.expectEqual(a_.add(a_.mul(b_)), a.dot(b));

    try testing.expectEqual(b_.add(b_).add(F.one), b.sum());

    const a_2 = a_.mul(a_);
    try testing.expectEqual(a_2.add(a_2).sqrt(), a.norm());

    try testing.expectEqual(a_.min(F.zero), a.min());
    try testing.expectEqual(a_.max(F.zero), a.max());
    try testing.expectEqual(b_.min(F.one), b.min());
    try testing.expectEqual(b_.max(F.one), b.max());
}
