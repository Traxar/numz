const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;

//TODO:
//transpose
//add
//sub
//sub_
//mul
// matrix
//LU
// solve

/// return struct that can create new matrices
/// operations by this struct result in a newly allocated matrix
pub fn MatrixType(comptime Scalar: type) type {
    return struct {
        const Self = @This();
        const VectorType = @import("vector.zig").VectorType(Scalar);
        const Vector = VectorType.Vector;
        allocptr: *const std.mem.Allocator,

        ///allocate empty matrix
        fn init(allocptr: *const std.mem.Allocator, rows: usize, cols: usize) !Matrix {
            var res = Matrix{
                .val = (try allocptr.alloc(Matrix.Row, rows)).ptr,
                .rows = rows,
                .cols = cols,
                .allocptr = allocptr,
            };
            for (0..res.rows) |i| {
                res.val[i] = Matrix.Row{};
            }
            return res;
        }
        /// allocate empty matrix with n rows and m columns
        pub fn zero(self: Self, rows: usize, cols: usize) !Matrix {
            return init(self.allocptr, rows, cols);
        }

        /// allocate the identity matrix of size n
        pub fn eye(self: Self, size: usize) !Matrix {
            var res = try self.zero(size, size);
            for (0..res.rows) |i| {
                try res.set(i, i, Scalar.eye);
            }
            return res;
        }

        //TODO: create matrix with 2 vectors
        //TODO: convert vector to matrix
        //TODO: convert scalar to matrix

        /// sparse matrix with runtime size
        const Matrix = struct {
            const Row = std.MultiArrayList(struct { col: usize, val: Scalar });
            val: [*]Row,
            rows: usize,
            cols: usize,
            allocptr: *const std.mem.Allocator,

            /// deinitialize matrix
            pub fn deinit(a: Matrix) void {
                for (0..a.rows) |i| {
                    a.val[i].deinit(a.allocptr.*);
                }
                a.allocptr.free(a.val[0..a.rows]);
            }

            /// b := a
            /// return b
            pub fn copy(a: Matrix) !Matrix {
                var res = try init(a.allocptr, a.rows, a.cols);
                for (0..a.rows) |i| {
                    res.val[i] = a.val[i].clone(a.allocptr.*);
                }
                return res;
            }

            /// b := a^T
            /// return b
            /// O(m*n)
            pub fn transpose(a: Matrix) !Matrix {
                var res = try init(a.allocptr, a.cols, a.rows);
                for (0..a.rows) |i| {
                    for (0..a.entAt(i)) |j| {
                        try res.set(a.colAt(i, j), i, a.valAt(i, j, true));
                    }
                }
                return res;
            }

            /// return index of col in row
            /// performs binary search
            /// O(log(m)), O(1) if dense or last element of row
            fn indAt(a: Matrix, row: usize, col: usize) struct { ind: usize, ex: bool } {
                const r = a.val[row].items(.col);
                if (r.len == 0 or col > r[r.len - 1]) { // end of row
                    return .{ .ind = r.len, .ex = false };
                } else if (col == r[r.len - 1]) { // last of row
                    return .{ .ind = r.len - 1, .ex = true };
                } else { // binary search
                    var min = r.len -| (a.cols - col);
                    var max = @min(r.len - 1, col); // - 1 because entry is not last of row
                    while (min < max) {
                        const pivot = @divFloor(min + max, 2);
                        if (col <= r[pivot]) {
                            max = pivot;
                        } else {
                            min = pivot + 1;
                        }
                    }
                    return .{ .ind = min, .ex = (min < r.len and col == r[min]) };
                }
            }

            /// return number of entries in the row
            fn entAt(a: Matrix, row: usize) usize {
                return a.val[row].len;
            }

            /// return element at row and columnIndex
            /// O(1)
            fn valAt(a: Matrix, row: usize, i: usize, exists: bool) Scalar {
                assert(row < a.rows);
                assert(i <= a.val[row].len);
                if (!exists) {
                    return Scalar.zero;
                } else {
                    return a.val[row].items(.val)[i];
                }
            }

            /// return column at row and columnIndex i
            /// O(1)
            fn colAt(a: Matrix, row: usize, i: usize) usize {
                assert(row < a.rows);
                assert(i < a.val[row].len);
                return a.val[row].items(.col)[i];
            }

            /// set element at row and columnIndex i
            /// O(1)
            fn setAt(a: Matrix, row: usize, i: usize, exists: bool, col: usize, b: Scalar) !void {
                assert(row < a.rows);
                assert(col < a.cols);
                assert(i <= a.val[row].len);
                if (b.cmp(.equal, Scalar.zero)) {
                    if (exists) {
                        a.val[row].orderedRemove(i);
                    } // else do nothing
                } else {
                    if (exists) {
                        a.val[row].items(.val)[i] = b;
                    } else {
                        try a.val[row].insert(a.allocptr.*, i, .{ .col = col, .val = b });
                    }
                }
            }

            /// return element at row and column
            /// O(log(m))
            pub fn at(a: Matrix, row: usize, col: usize) Scalar {
                assert(row < a.rows);
                assert(col < a.cols);
                const i = a.indAt(row, col);
                if (i.ex) {
                    return a.val[row].items(.val)[i.ind];
                } else {
                    return Scalar.zero;
                }
            }

            // set element at row and column to b
            // O(m)
            // O(1) dense of last column in row
            pub fn set(a: Matrix, row: usize, col: usize, b: Scalar) !void {
                assert(row < a.rows);
                assert(col < a.cols);
                const i = a.indAt(row, col);
                try a.setAt(row, i.ind, i.ex, col, b);
                return;
            }

            /// a <- a * b
            /// return a
            /// O(n*m)
            pub fn mulS(a: Matrix, b: Scalar) Matrix {
                for (0..a.rows) |i| {
                    const e = a.entAt(i);
                    for (1..e + 1) |j_| {
                        const j = e - j_;
                        a.setAt(i, j, true, 0, a.valAt(i, j, true).mul(b)) catch unreachable;
                    }
                }
                return a;
            }

            /// c := a * b
            /// return c
            /// O(n*m)
            pub fn mulV(a: Matrix, b: Vector) !Vector {
                assert(a.cols == b.len);
                var res = try (VectorType{ .allocptr = b.allocptr }).rep(Scalar.zero, a.rows);
                for (0..a.rows) |i| {
                    for (0..a.entAt(i)) |j| {
                        res.set(i, res.at(i).add(a.valAt(i, j, true).mul(b.at(a.colAt(i, j)))));
                    }
                }
                return res;
            }

            /// c := a * b
            /// return c
            /// O(n*m^2)
            pub fn mul(a: Matrix, b: Matrix) !Matrix {
                assert(a.cols == b.rows);
                var res = try init(a.allocptr, a.rows, b.cols);
                for (0..a.rows) |i| {
                    for (0..a.entAt(i)) |j| {
                        const r = a.colAt(i, j);
                        for (0..b.entAt(r)) |k| {
                            const c = b.colAt(r, k);
                            const l = res.indAt(i, c);
                            try res.setAt(i, l.ind, l.ex, c, res.valAt(i, l.ind, l.ex).add(a.valAt(i, j, true).mul(b.valAt(r, k, true))));
                        }
                    }
                }
                return res;
            }

            /// a <- a * b
            /// return a
            /// O(n*m)
            pub fn divS(a: Matrix, b: Scalar) Matrix {
                for (0..a.rows) |i| {
                    const e = a.entAt(i);
                    for (1..e + 1) |j_| {
                        const j = e - j_;
                        a.setAt(i, j, true, 0, a.valAt(i, j, true).div(b)) catch unreachable;
                    }
                }
                return a;
            }
        };
    };
}

test "creation" {
    const allocator = std.testing.allocator;
    const F = @import("../scalar.zig").Float(f32);
    const M = MatrixType(F){ .allocptr = &allocator };

    const n = 5;
    const m = 8;
    var a = try M.zero(n, m);
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

test "removing entries" {
    const allocator = std.testing.allocator;
    const F = @import("../scalar.zig").Float(f32);
    const M = MatrixType(F){ .allocptr = &allocator };

    var a = try M.zero(3, 3);
    defer a.deinit();
    // 2 1 3
    // 6 0 5
    // 8 9 7

    try a.set(0, 1, F.from(1));
    try a.set(0, 0, F.from(2));
    try a.set(0, 2, F.from(3));
    try a.set(1, 1, F.zero);
    try a.set(1, 2, F.from(5));
    try a.set(1, 0, F.from(6));
    try a.set(2, 2, F.from(7));
    try a.set(2, 0, F.from(8));
    try a.set(2, 1, F.from(9));

    try testing.expect(a.at(0, 0).cmp(.equal, F.from(2)));
    try testing.expect(a.at(0, 1).cmp(.equal, F.from(1)));
    try testing.expect(a.at(0, 2).cmp(.equal, F.from(3)));
    try testing.expect(a.at(1, 0).cmp(.equal, F.from(6)));
    try testing.expect(a.at(1, 1).cmp(.equal, F.from(0)));
    try testing.expect(a.at(1, 2).cmp(.equal, F.from(5)));
    try testing.expect(a.at(2, 0).cmp(.equal, F.from(8)));
    try testing.expect(a.at(2, 1).cmp(.equal, F.from(9)));
    try testing.expect(a.at(2, 2).cmp(.equal, F.from(7)));

    try testing.expectEqual(@as(usize, 3), a.val[0].len);
    try testing.expectEqual(@as(usize, 2), a.val[1].len);

    try a.set(2, 2, F.from(0));
    try a.set(2, 0, F.from(0));
    try a.set(2, 1, F.from(0));

    try testing.expectEqual(@as(usize, 0), a.val[2].len);
}

test "mulpiplication with scalar" {
    const allocator = std.testing.allocator;
    const F = @import("../scalar.zig").Float(f32);
    const M = MatrixType(F){ .allocptr = &allocator };

    const n = 3;
    const a_ = F.from(-3.14);
    var a = (try M.eye(n)).mulS(a_);
    defer a.deinit();
    for (0..n) |i| {
        for (0..n) |j| {
            if (i == j) {
                try testing.expect(a.at(i, j).cmp(.equal, a_));
            } else {
                try testing.expect(a.at(i, j).cmp(.equal, F.zero));
            }
        }
    }
    _ = a.divS(a_);
    for (0..n) |i| {
        for (0..n) |j| {
            if (i == j) {
                try testing.expect(a.at(i, j).cmp(.equal, a_.div(a_)));
            } else {
                try testing.expect(a.at(i, j).cmp(.equal, F.zero));
            }
        }
    }
}

test "mulpiplication with vector" {
    const allocator = std.testing.allocator;
    const F = @import("../scalar.zig").Float(f32);
    const V = @import("vector.zig").VectorType(F){ .allocptr = &allocator };
    const M = MatrixType(F){ .allocptr = &allocator };

    const n = 3;
    var a = try M.eye(n);
    defer a.deinit();
    try a.set(0, 2, F.from(-1));
    try a.set(1, 0, F.from(-1));
    try a.set(2, 1, F.from(-1));

    var v = try V.rep(F.zero, n);
    defer v.deinit();
    v.set(1, F.from(1));
    v.set(2, F.from(2));

    const av = try a.mulV(v);
    defer av.deinit();
    try testing.expect(av.at(0).cmp(.equal, F.from(-2)));
    try testing.expect(av.at(1).cmp(.equal, F.from(1)));
    try testing.expect(av.at(2).cmp(.equal, F.from(1)));
}

test "transpose" {
    const allocator = std.testing.allocator;
    const F = @import("../scalar.zig").Float(f32);
    const M = MatrixType(F){ .allocptr = &allocator };

    var a = try M.zero(3, 3);
    defer a.deinit();
    // 2 1 3
    // 6 0 5
    // 8 9 7

    try a.set(0, 1, F.from(1));
    try a.set(0, 0, F.from(2));
    try a.set(0, 2, F.from(3));
    try a.set(1, 1, F.zero);
    try a.set(1, 2, F.from(5));
    try a.set(1, 0, F.from(6));
    try a.set(2, 2, F.from(7));
    try a.set(2, 0, F.from(8));
    try a.set(2, 1, F.from(9));

    var b = try a.transpose();
    defer b.deinit();

    try testing.expect(b.at(0, 0).cmp(.equal, F.from(2)));
    try testing.expect(b.at(1, 0).cmp(.equal, F.from(1)));
    try testing.expect(b.at(2, 0).cmp(.equal, F.from(3)));
    try testing.expect(b.at(0, 1).cmp(.equal, F.from(6)));
    try testing.expect(b.at(1, 1).cmp(.equal, F.from(0)));
    try testing.expect(b.at(2, 1).cmp(.equal, F.from(5)));
    try testing.expect(b.at(0, 2).cmp(.equal, F.from(8)));
    try testing.expect(b.at(1, 2).cmp(.equal, F.from(9)));
    try testing.expect(b.at(2, 2).cmp(.equal, F.from(7)));
}

test "matrix multiplication" {
    const allocator = std.testing.allocator;
    const F = @import("../scalar.zig").Float(f32);
    const M = MatrixType(F){ .allocptr = &allocator };

    // 2 1 3   1 0 1   2 1 6
    // 6 0 5 * 0 1 1 = 6 0 11
    // 8 9 7   0 0 1   8 9 24

    var a = try M.zero(3, 3);
    defer a.deinit();

    try a.set(0, 1, F.from(1));
    try a.set(0, 0, F.from(2));
    try a.set(0, 2, F.from(3));
    try a.set(1, 1, F.zero);
    try a.set(1, 2, F.from(5));
    try a.set(1, 0, F.from(6));
    try a.set(2, 2, F.from(7));
    try a.set(2, 0, F.from(8));
    try a.set(2, 1, F.from(9));

    var b = try M.eye(3);
    defer b.deinit();

    try b.set(0, 2, F.from(1));
    try b.set(1, 2, F.from(1));

    var c = try a.mul(b);
    defer c.deinit();

    try testing.expect(c.at(0, 0).cmp(.equal, F.from(2)));
    try testing.expect(c.at(0, 1).cmp(.equal, F.from(1)));
    try testing.expect(c.at(0, 2).cmp(.equal, F.from(6)));
    try testing.expect(c.at(1, 0).cmp(.equal, F.from(6)));
    try testing.expect(c.at(1, 1).cmp(.equal, F.from(0)));
    try testing.expect(c.at(1, 2).cmp(.equal, F.from(11)));
    try testing.expect(c.at(2, 0).cmp(.equal, F.from(8)));
    try testing.expect(c.at(2, 1).cmp(.equal, F.from(9)));
    try testing.expect(c.at(2, 2).cmp(.equal, F.from(24)));
}
