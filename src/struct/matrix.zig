const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;

/// return struct that can create Vectors
pub fn MatrixType(comptime Scalar: type) type {
    return struct {
        const Self = @This();
        const ArrayHashMap = std.AutoArrayHashMap(struct { i: usize, j: usize }, Scalar);
        allocator: std.mem.Allocator,

        /// create an empty matrix with n rows and m columns
        /// do not forget to call deinit() on result
        pub fn zero(self: Self, n: usize, m: usize) Matrix {
            var res = Matrix{
                .entries = ArrayHashMap.init(self.allocator),
                .n = n,
                .m = m,
            };
            return res;
        }

        /// create the identity matrix of size n
        /// do not forget to call deinit() on result
        pub fn eye(self: Self, n: usize) !Matrix {
            var res = self.zero(n, n);
            for (0..n) |i| {
                try res.set(i, i, Scalar.eye);
            }
            return res;
        }

        /// vector with runtime length
        const Matrix = struct {
            entries: ArrayHashMap,
            n: usize, //amount of rows
            m: usize, //amount of columns

            /// deinitialize Vector
            pub fn deinit(a: *Matrix) void {
                a.entries.deinit();
            }

            /// return element at row i and column j
            pub fn at(a: Matrix, i: usize, j: usize) Scalar {
                assert(i < a.n);
                assert(j < a.m);
                return a.entries.get(.{ .i = i, .j = j }) orelse Scalar.zero;
            }

            // set element at row i and column j to b
            pub fn set(a: *Matrix, i: usize, j: usize, b: Scalar) !void {
                assert(i < a.n);
                assert(j < a.m);
                if (b.cmp(.equal, Scalar.zero)) {
                    _ = a.entries.swapRemove(.{ .i = i, .j = j });
                } else {
                    try a.entries.put(.{ .i = i, .j = j }, b);
                }
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

    const F = @import("../scalar.zig").Float(f32);
    const M = MatrixType(F){ .allocator = allocator };

    const n = 5;
    const m = 8;
    var a = M.zero(n, m);
    defer a.deinit();
    for (0..n) |i| {
        for (0..m) |j| {
            try testing.expect(a.at(i, j).cmp(.equal, F.zero));
        }
    }

    var b = try M.eye(n);
    defer b.deinit();
    for (0..n) |i| {
        for (0..n) |j| {
            if (i == j) {
                try testing.expect(b.at(i, j).cmp(.equal, F.eye));
            } else {
                try testing.expect(b.at(i, j).cmp(.equal, F.zero));
            }
        }
    }
}
