const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;
const Allocator = std.mem.Allocator;

/// return struct that can allocate and free vectors
/// operations by this struct create new vectors
pub fn VectorType(comptime Scalar: type) type {
    const SIMDScalar = Scalar.SIMDType(null);
    const SIMDsize = SIMDScalar.SIMDsize;

    return struct {
        const Vector = @This();
        val: std.MultiArrayList(Scalar).Slice,

        ///allocate vector with undefined values
        pub fn init(n: usize, allocator: Allocator) !Vector {
            var multi_array_list = std.MultiArrayList(Scalar){};
            try multi_array_list.ensureTotalCapacity(allocator, n);
            multi_array_list.len = n;
            var slice = multi_array_list.slice();
            return Vector{ .val = slice };
        }

        /// allocate vector with n elements with value a
        pub fn rep(a: Scalar, n: usize, allocator: Allocator) !Vector {
            var res = try Vector.init(n, allocator);
            const ii = n - @mod(n, SIMDsize);
            var i: usize = 0;
            if (ii >= 1) {
                const SIMDa = a.SIMDsplat(SIMDScalar);
                while (i < ii) : (i += SIMDsize) {
                    SIMDa.toVec(res.val, i);
                }
            }
            while (i < n) : (i += 1) {
                a.toVec(res.val, i);
            }
            return res;
        }

        /// deinitialize vector
        pub fn deinit(a: *Vector, allocator: Allocator) void {
            a.val.deinit(allocator);
        }

        /// TODO SIMDify
        /// allocate copy of vector a
        pub fn copy(a: Vector, allocator: Allocator) !Vector {
            return Vector{
                .val = (try allocator.dupe(Scalar, a.val[0..a.len])).ptr,
                .len = a.len,
            };
        }

        /// return element at index i
        pub fn at(a: Vector, i: usize) Scalar {
            return Scalar.fromVec(a.val, i);
        }

        /// set element at index i to b
        pub fn set(a: Vector, i: usize, b: Scalar) void {
            b.toVe(a.val, i);
        }

        /// TODO SIMDify
        /// a <- a + b
        /// return a
        pub fn add(a: Vector, b: Vector) Vector {
            assert(a.len == b.len);
            for (0..a.len) |i| {
                a.set(i, a.at(i).add(b.at(i)));
            }
            return a;
        }

        /// TODO SIMDify
        /// a <- a - b
        /// return a
        pub fn sub(a: Vector, b: Vector) Vector {
            assert(a.len == b.len);
            for (0..a.len) |i| {
                a.set(i, a.at(i).sub(b.at(i)));
            }
            return a;
        }

        /// TODO SIMDify
        /// a <- b - a
        /// return a
        pub fn sub_(a: Vector, b: Vector) Vector { //TODO: change name
            assert(a.len == b.len);
            for (0..a.len) |i| {
                a.set(i, b.at(i).sub(a.at(i)));
            }
            return a;
        }

        /// TODO SIMDify
        /// a <- a * b
        /// return a
        pub fn mulS(a: Vector, b: Scalar) Vector {
            for (0..a.len) |i| {
                a.set(i, a.at(i).mul(b));
            }
            return a;
        }

        /// TODO SIMDify
        /// a <- a * b
        /// elemtentwise
        /// return a
        pub fn mul(a: Vector, b: Vector) Vector {
            assert(a.len == b.len);
            for (0..a.len) |i| {
                a.set(i, a.at(i).mul(b.at(i)));
            }
            return a;
        }

        /// TODO SIMDify
        /// a <- a / b
        /// return a
        pub fn divS(a: Vector, b: Scalar) Vector {
            for (0..a.len) |i| {
                a.set(i, a.at(i).div(b));
            }
            return a;
        }

        /// TODO SIMDify
        /// a <- a / b
        /// elementwise
        /// return a
        pub fn div(a: Vector, b: Vector) Vector {
            assert(a.len == b.len);
            for (0..a.len) |i| {
                a.set(i, a.at(i).div(b.at(i)));
            }
            return a;
        }

        /// TODO SIMDify
        /// a <- b / a
        /// elementwise
        /// return a
        pub fn div_(a: Vector, b: Vector) Vector { //TODO: change name
            assert(a.len == b.len);
            for (0..a.len) |i| {
                a.set(i, b.at(i).div(a.at(i)));
            }
            return a;
        }

        /// TODO SIMDify
        /// return sum of elements
        pub fn sum(a: Vector) Scalar {
            var res = Scalar.zero;
            for (0..a.len) |i| {
                res = res.add(a.at(i));
            }
            return res;
        }

        /// TODO SIMDify
        ///return the scalar product of a and b
        pub fn dot(a: Vector, b: Vector) Scalar {
            var res = Scalar.zero;
            for (0..a.len) |i| {
                res = res.add(a.at(i).mul(b.at(i)));
            }
            return res;
        }

        /// TODO SIMDify
        ///return the euklidean norm
        pub fn norm(a: Vector) Scalar {
            return a.dot(a).sqrt();
        }

        /// TODO SIMDify
        ///return the smallest element
        pub fn min(a: Vector) Scalar {
            var res = a.at(0);
            for (1..a.len) |i| {
                res = res.min(a.at(i));
            }
            return res;
        }

        /// TODO SIMDify
        ///return the largest element
        pub fn max(a: Vector) Scalar {
            var res = a.at(0);
            for (1..a.len) |i| {
                res = res.max(a.at(i));
            }
            return res;
        }

        pub fn debugprint(self: Vector) void {
            std.debug.print("\n", .{});
            for (0..self.len) |i| {
                std.debug.print("{}\n", .{self.at(i).f});
            }
        }
    };
}

test "creation" {
    const ally = std.testing.allocator;
    const n = 100;
    const F = @import("scalar.zig").Float(f32);
    try testing.expect(F.SIMDsize == 1);
    const V = VectorType(F);

    const a = F.from(-314, 100);
    var v = try V.rep(a, n, ally);
    //var v = try V.init(n, ally);
    defer v.deinit(ally);

    try testing.expectEqual(@as(usize, n), v.val.len);
    for (0..n) |i| {
        try testing.expect(a.cmp(.eq, v.at(i)));
    }
}

// test "operators" {
//     const ally = std.testing.allocator;
//     const n = 3;
//     const F = @import("../scalar.zig").Float(f32);
//     const V = VectorType(F);

//     const a_ = F.from(-314, 100);
//     const b_ = F.from(527, 100);
//     var a = try V.rep(a_, n, ally);
//     defer a.deinit(ally);
//     a.set(1, F.zero);
//     var b = try V.rep(b_, n, ally);
//     defer b.deinit(ally);
//     b.set(0, F.eye);

//     try testing.expect(a.at(1).cmp(.eq, F.zero));

//     const c = (try a.copy(ally)).add(b);
//     defer c.deinit(ally);
//     try testing.expect(c.at(0).cmp(.eq, a_.add(F.eye)));
//     try testing.expect(c.at(1).cmp(.eq, b_));
//     try testing.expect(c.at(2).cmp(.eq, a_.add(b_)));

//     const d = (try a.copy(ally)).add(b);
//     defer d.deinit(ally);
//     try testing.expect(d.at(0).cmp(.eq, a_.add(F.eye)));
//     try testing.expect(d.at(1).cmp(.eq, b_));
//     try testing.expect(d.at(2).cmp(.eq, a_.add(b_)));

//     const e = (try a.copy(ally)).mul(b);
//     defer e.deinit(ally);
//     try testing.expect(e.at(0).cmp(.eq, a_));
//     try testing.expect(e.at(1).cmp(.eq, F.zero));
//     try testing.expect(e.at(2).cmp(.eq, a_.mul(b_)));

//     const f = (try a.copy(ally)).div(b);
//     defer f.deinit(ally);
//     try testing.expect(f.at(0).cmp(.eq, a_));
//     try testing.expect(f.at(1).cmp(.eq, F.zero));
//     try testing.expect(f.at(2).cmp(.eq, a_.div(b_)));

//     try testing.expect(a.dot(b).cmp(.eq, a_.add(a_.mul(b_))));

//     try testing.expect(b.sum().cmp(.eq, b_.add(b_).add(F.eye)));

//     const a_2 = a_.mul(a_);
//     try testing.expect(a.norm().cmp(.eq, a_2.add(a_2).sqrt()));

//     try testing.expect(a.min().cmp(.eq, a_.min(F.zero)));
//     try testing.expect(b.max().cmp(.eq, b_.max(F.eye)));
// }
