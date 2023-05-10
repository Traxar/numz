const std = @import("std");
const testing = std.testing;

/// return struct that can create Vectors
pub fn VectorType(comptime Scalar: type) type {
    return struct {
        const Self = @This();
        const ArrayList = std.ArrayList(Scalar);
        allocator: std.mem.Allocator,

        /// create a vector with n elements with value a
        /// do not forget to call deinit() on result
        pub fn rep(self: Self, a: Scalar, n: usize) !Vector {
            var res = Vector{ .v = try ArrayList.initCapacity(self.allocator, n) };
            res.v.appendNTimesAssumeCapacity(a, n);
            return res;
        }

        /// vector with runtime length
        const Vector = struct {
            v: ArrayList,

            /// deinitialize Vector
            pub fn deinit(a: Vector) void {
                a.v.deinit();
            }

            /// create a copy of this vector, using the same allocator
            /// do not forget to call deinit() on result
            pub fn copy(a: Vector) !Vector {
                return Vector{ .v = try ArrayList.clone(a.v) };
            }

            /// return length of vector
            pub fn len(a: Vector) usize {
                return a.v.items.len;
            }

            /// return element at index i
            pub fn at(a: Vector, i: usize) Scalar {
                return a.v.items[i];
            }

            /// set element at index i to b
            pub fn set(a: Vector, i: usize, b: Scalar) void {
                a.v.items[i] = b;
            }

            /// add b to a
            /// return a
            pub fn add(a: Vector, b: Vector) Vector {
                for (a.v.items, b.v.items, 0..) |a_, b_, i| {
                    a.set(i, a_.add(b_));
                }
                return a;
            }

            /// subtract b from a
            /// return a
            pub fn sub(a: Vector, b: Vector) Vector {
                for (a.v.items, b.v.items, 0..) |a_, b_, i| {
                    a.set(i, a_.sub(b_));
                }
                return a;
            }

            /// elemtentwise multiply b onto a
            /// return a
            pub fn mul(a: Vector, b: Vector) Vector {
                for (a.v.items, b.v.items, 0..) |a_, b_, i| {
                    a.set(i, a_.mul(b_));
                }
                return a;
            }

            /// elementwise divide a by b
            /// return a
            pub fn div(a: Vector, b: Vector) Vector {
                for (a.v.items, b.v.items, 0..) |a_, b_, i| {
                    a.set(i, a_.div(b_));
                }
                return a;
            }

            /// return sum of elements
            pub fn sum(a: Vector) Scalar {
                var res = Scalar.zero;
                for (a.v.items) |a_| {
                    res = res.add(a_);
                }
                return res;
            }

            ///return elementwise a / b
            pub fn dot(a: Vector, b: Vector) Scalar {
                var res = Scalar.zero;
                for (a.v.items, b.v.items) |a_, b_| {
                    res = res.add(a_.mul(b_));
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
                for (1..a.len()) |i| {
                    res = res.min(a.at(i));
                }
                return res;
            }

            ///return the largest element
            pub fn max(a: Vector) Scalar {
                var res = a.at(0);
                for (1..a.len()) |i| {
                    res = res.max(a.at(i));
                }
                return res;
            }

            /// append b to a
            pub fn append(a: *Vector, b: Scalar) !void {
                try a.v.append(b);
            }

            /// concat b onto a
            pub fn concat(a: *Vector, b: Vector) !void {
                try a.v.appendSlice(b.v.items);
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
    defer v.deinit();

    try testing.expectEqual(@as(usize, n), v.len());
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
    defer a.deinit();
    a.set(1, F.zero);
    var b = try V.rep(b_, n);
    defer b.deinit();
    b.set(0, F.eye);

    try testing.expect(a.at(1).cmp(.equal, F.zero));

    const c = (try a.copy()).add(b);
    defer c.deinit();
    try testing.expect(c.at(0).cmp(.equal, a_.add(F.eye)));
    try testing.expect(c.at(1).cmp(.equal, b_));
    try testing.expect(c.at(2).cmp(.equal, a_.add(b_)));

    const d = (try a.copy()).add(b);
    defer d.deinit();
    try testing.expect(d.at(0).cmp(.equal, a_.add(F.eye)));
    try testing.expect(d.at(1).cmp(.equal, b_));
    try testing.expect(d.at(2).cmp(.equal, a_.add(b_)));

    const e = (try a.copy()).mul(b);
    defer e.deinit();
    try testing.expect(e.at(0).cmp(.equal, a_));
    try testing.expect(e.at(1).cmp(.equal, F.zero));
    try testing.expect(e.at(2).cmp(.equal, a_.mul(b_)));

    const f = (try a.copy()).div(b);
    defer f.deinit();
    try testing.expect(f.at(0).cmp(.equal, a_));
    try testing.expect(f.at(1).cmp(.equal, F.zero));
    try testing.expect(f.at(2).cmp(.equal, a_.div(b_)));

    try testing.expect(a.dot(b).cmp(.equal, a_.add(a_.mul(b_))));

    try testing.expect(b.sum().cmp(.equal, b_.add(b_).add(F.eye)));

    const a_2 = a_.mul(a_);
    try testing.expect(a.norm().cmp(.equal, a_2.add(a_2).sqrt()));

    try testing.expect(a.min().cmp(.equal, a_.min(F.zero)));
    try testing.expect(b.max().cmp(.equal, b_.max(F.eye)));

    try a.append(F.eye);
    try testing.expect(a.at(3).cmp(.equal, F.eye));

    try a.concat(b);
    try testing.expect(a.at(6).cmp(.equal, b_));
}
