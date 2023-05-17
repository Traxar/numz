const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;

/// return struct that can allocate and free vectors
/// operations by this struct result in a newly allocated vector
/// operation by the vector struct are inplace
pub fn VectorType(comptime Scalar: type) type {
    return struct {
        const Self = @This();
        allocator: std.mem.Allocator,

        fn init(self: Self, n: usize) !Vector {
            return Vector{
                .val = (try self.allocator.alloc(Scalar, n)).ptr,
                .len = n,
            };
        }

        /// deinitialize Vector
        pub fn deinit(self: Self, a: Vector) void {
            self.allocator.free(a.val[0..a.len]);
        }

        /// allocates a vector with n elements with value a
        pub fn rep(self: Self, a: Scalar, n: usize) !Vector {
            var res = try self.init(n);
            for (0..n) |i| {
                res.set(i, a);
            }
            return res;
        }

        /// allocates a copy of vector a
        pub fn copy(self: Self, a: Vector) !Vector {
            return Vector{
                .val = (try self.allocator.dupe(Scalar, a.val[0..a.len])).ptr,
                .len = a.len,
            };
        }

        /// vector with runtime length
        /// all functions in this struct are inplace and do not require memory allocation
        pub const Vector = struct {
            val: [*]Scalar,
            len: usize,

            /// return element at index i
            pub fn at(a: Vector, i: usize) Scalar {
                assert(i < a.len);
                return a.val[i];
            }

            /// set element at index i to b
            pub fn set(a: Vector, i: usize, b: Scalar) void {
                assert(i < a.len);
                a.val[i] = b;
            }

            /// a <- a + b
            /// return a
            pub fn add(a: Vector, b: Vector) Vector {
                assert(a.len == b.len);
                for (0..a.len) |i| {
                    a.set(i, a.at(i).add(b.at(i)));
                }
                return a;
            }

            /// a <- a - b
            /// return a
            pub fn sub(a: Vector, b: Vector) Vector {
                assert(a.len == b.len);
                for (0..a.len) |i| {
                    a.set(i, a.at(i).sub(b.at(i)));
                }
                return a;
            }

            /// a <- b - a
            /// return a
            pub fn sub_(a: Vector, b: Vector) Vector { //TODO: change name
                assert(a.len == b.len);
                for (0..a.len) |i| {
                    a.set(i, b.at(i).sub(a.at(i)));
                }
                return a;
            }

            /// a <- a * b
            /// return a
            pub fn mulS(a: Vector, b: Scalar) Vector {
                for (0..a.len) |i| {
                    a.set(i, a.at(i).mul(b));
                }
                return a;
            }

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

            /// a <- a / b
            /// return a
            pub fn divS(a: Vector, b: Scalar) Vector {
                for (0..a.len) |i| {
                    a.set(i, a.at(i).div(b));
                }
                return a;
            }

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

            /// return sum of elements
            pub fn sum(a: Vector) Scalar {
                var res = Scalar.zero;
                for (0..a.len) |i| {
                    res = res.add(a.at(i));
                }
                return res;
            }

            ///return the scalar product of a and b
            pub fn dot(a: Vector, b: Vector) Scalar {
                var res = Scalar.zero;
                for (0..a.len) |i| {
                    res = res.add(a.at(i).mul(b.at(i)));
                }
                return res;
            }

            ///return the euklidean norm
            pub fn norm(a: Vector) Scalar {
                return a.dot(a).sqrt();
            }

            ///return the smallest element
            pub fn min(a: Vector) Scalar {
                var res = a.at(0);
                for (1..a.len) |i| {
                    res = res.min(a.at(i));
                }
                return res;
            }

            ///return the largest element
            pub fn max(a: Vector) Scalar {
                var res = a.at(0);
                for (1..a.len) |i| {
                    res = res.max(a.at(i));
                }
                return res;
            }
        };
    };
}

test "creation" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();
    defer {
        const deinit_status = gpa.deinit();
        //fail test; can't try in defer as defer is executed after we return
        if (deinit_status == .leak) testing.expect(false) catch @panic("TEST FAIL");
    }

    const n = 100;
    const F = @import("../scalar.zig").Float(f32);
    const V = VectorType(F){ .allocator = allocator };

    const a = F.from(-3.14);
    var v = try V.rep(a, n);
    defer V.deinit(v);

    try testing.expectEqual(@as(usize, n), v.len);
    for (0..n) |i| {
        try testing.expect(a.cmp(.equal, v.at(i)));
    }
}

test "operators" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();
    defer {
        const deinit_status = gpa.deinit();
        //fail test; can't try in defer as defer is executed after we return
        if (deinit_status == .leak) testing.expect(false) catch @panic("TEST FAIL");
    }

    const n = 3;
    const F = @import("../scalar.zig").Float(f32);
    const V = VectorType(F){ .allocator = allocator };

    const a_ = F.from(-3.14);
    const b_ = F.from(5.27);
    var a = try V.rep(a_, n);
    defer V.deinit(a);
    a.set(1, F.zero);
    var b = try V.rep(b_, n);
    defer V.deinit(b);
    b.set(0, F.eye);

    try testing.expect(a.at(1).cmp(.equal, F.zero));

    const c = (try V.copy(a)).add(b);
    defer V.deinit(c);
    try testing.expect(c.at(0).cmp(.equal, a_.add(F.eye)));
    try testing.expect(c.at(1).cmp(.equal, b_));
    try testing.expect(c.at(2).cmp(.equal, a_.add(b_)));

    const d = (try V.copy(a)).add(b);
    defer V.deinit(d);
    try testing.expect(d.at(0).cmp(.equal, a_.add(F.eye)));
    try testing.expect(d.at(1).cmp(.equal, b_));
    try testing.expect(d.at(2).cmp(.equal, a_.add(b_)));

    const e = (try V.copy(a)).mul(b);
    defer V.deinit(e);
    try testing.expect(e.at(0).cmp(.equal, a_));
    try testing.expect(e.at(1).cmp(.equal, F.zero));
    try testing.expect(e.at(2).cmp(.equal, a_.mul(b_)));

    const f = (try V.copy(a)).div(b);
    defer V.deinit(f);
    try testing.expect(f.at(0).cmp(.equal, a_));
    try testing.expect(f.at(1).cmp(.equal, F.zero));
    try testing.expect(f.at(2).cmp(.equal, a_.div(b_)));

    try testing.expect(a.dot(b).cmp(.equal, a_.add(a_.mul(b_))));

    try testing.expect(b.sum().cmp(.equal, b_.add(b_).add(F.eye)));

    const a_2 = a_.mul(a_);
    try testing.expect(a.norm().cmp(.equal, a_2.add(a_2).sqrt()));

    try testing.expect(a.min().cmp(.equal, a_.min(F.zero)));
    try testing.expect(b.max().cmp(.equal, b_.max(F.eye)));
}
