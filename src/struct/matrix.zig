const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;
const Allocator = std.mem.Allocator;

//TODO:
// + add
// + sub
// + LU
// + solve

/// return struct that can create new matrices
/// operations by this struct result in a newly allocated matrix
pub fn MatrixType(comptime Scalar: type) type {
    return struct {
        const Matrix = @This();
        const Vector = @import("vector.zig").VectorType(Scalar);
        const Row = std.MultiArrayList(struct { col: usize, val: Scalar });
        val: [*]Row,
        rows: usize,
        cols: usize,
        allocator: Allocator,


        ///allocate empty matrix
        pub fn zero(rows: usize, cols: usize, allocator: Allocator) !Matrix {
            var res = Matrix{
                .val = (try allocator.alloc(Matrix.Row, rows)).ptr,
                .rows = rows,
                .cols = cols,
                .allocator = allocator,
            };
            for (0..res.rows) |i| {
                res.val[i] = Matrix.Row{};
            }
            return res;
        }

        /// allocate the identity matrix (square)
        pub fn eye(size: usize, allocator: Allocator) !Matrix {
            var res = try zero(size, size, allocator);
            for (0..res.rows) |i| {
                try res.set(i, i, Scalar.eye);
            }
            return res;
        }

        //TODO: create matrix with 2 vectors
        //TODO: convert vector to matrix
        //TODO: convert scalar to matrix

        /// deinitialize matrix
        pub fn deinit(a: Matrix) void {
            for (0..a.rows) |i| {
                a.val[i].deinit(a.allocator);
            }
            a.allocator.free(a.val[0..a.rows]);
        }

        /// b := a
        /// return b
        pub fn copy(a: Matrix, allocator: Allocator) !Matrix {
            var res = try zero(a.rows, a.cols, allocator);
            for (0..a.rows) |i| {
                res.val[i] = a.val[i].clone(allocator);
            }
            return res;
        }

        /// b := a^T
        /// return b
        /// O(m*n)
        pub fn transpose(a: Matrix, allocator: Allocator) !Matrix {
            var res = try zero(a.cols, a.rows, allocator);
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
                    try a.val[row].insert(a.allocator, i, .{ .col = col, .val = b });
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
        pub fn mulV(a: Matrix, b: Vector, allocator: Allocator) !Vector {
            assert(a.cols == b.len);
            var res = try Vector.rep(Scalar.zero, a.rows, allocator);
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
        pub fn mul(a: Matrix, b: Matrix, allocator: Allocator) !Matrix {
            assert(a.cols == b.rows);
            var res = try zero(a.rows, b.cols, allocator);
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

        /// a <- a / b
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
}

test "creation" {
    const ally = std.testing.allocator;
    const F = @import("../scalar.zig").Float(f32);
    const M = MatrixType(F);

    const n = 5;
    const m = 8;
    var a = try M.zero(n, m, ally);
    defer a.deinit();
    for (0..n) |i| {
        for (0..m) |j| {
            try testing.expect(a.at(i, j).cmp(.equal, F.zero));
        }
    }

    var b = try M.eye(n, ally);
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
    const ally = std.testing.allocator;
    const F = @import("../scalar.zig").Float(f32);
    const M = MatrixType(F);

    var a = try M.zero(3, 3, ally);
    defer a.deinit();
    // 2 1 3
    // 6 0 5
    // 8 9 7

    try a.set(0, 1, F.from(1,1));
    try a.set(0, 0, F.from(2,1));
    try a.set(0, 2, F.from(3,1));
    try a.set(1, 1, F.zero);
    try a.set(1, 2, F.from(5,1));
    try a.set(1, 0, F.from(6,1));
    try a.set(2, 2, F.from(7,1));
    try a.set(2, 0, F.from(8,1));
    try a.set(2, 1, F.from(9,1));

    try testing.expect(a.at(0, 0).cmp(.equal, F.from(2,1)));
    try testing.expect(a.at(0, 1).cmp(.equal, F.from(1,1)));
    try testing.expect(a.at(0, 2).cmp(.equal, F.from(3,1)));
    try testing.expect(a.at(1, 0).cmp(.equal, F.from(6,1)));
    try testing.expect(a.at(1, 1).cmp(.equal, F.from(0,1)));
    try testing.expect(a.at(1, 2).cmp(.equal, F.from(5,1)));
    try testing.expect(a.at(2, 0).cmp(.equal, F.from(8,1)));
    try testing.expect(a.at(2, 1).cmp(.equal, F.from(9,1)));
    try testing.expect(a.at(2, 2).cmp(.equal, F.from(7,1)));

    try testing.expectEqual(@as(usize, 3), a.val[0].len);
    try testing.expectEqual(@as(usize, 2), a.val[1].len);

    try a.set(2, 2, F.zero);
    try a.set(2, 0, F.zero);
    try a.set(2, 1, F.zero);

    try testing.expectEqual(@as(usize, 0), a.val[2].len);
}

test "mulpiplication with scalar" {
    const ally = std.testing.allocator;
    const F = @import("../scalar.zig").Float(f32);
    const M = MatrixType(F);

    const n = 3;
    const a_ = F.from(-314,100);
    var a = (try M.eye(n, ally)).mulS(a_);
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
    const ally= std.testing.allocator;
    const F = @import("../scalar.zig").Float(f32);
    const V = @import("vector.zig").VectorType(F);
    const M = MatrixType(F);

    const n = 3;
    var a = try M.eye(n, ally);
    defer a.deinit();
    try a.set(0, 2, F.from(-1, 1));
    try a.set(1, 0, F.from(-1, 1));
    try a.set(2, 1, F.from(-1, 1));

    var v = try V.rep(F.zero, n, ally);
    defer v.deinit(ally);
    v.set(1, F.from(1,1));
    v.set(2, F.from(2,1));

    const av = try a.mulV(v, ally);
    defer av.deinit(ally);
    try testing.expect(av.at(0).cmp(.equal, F.from(-2,1)));
    try testing.expect(av.at(1).cmp(.equal, F.from(1,1)));
    try testing.expect(av.at(2).cmp(.equal, F.from(1,1)));
}

test "transpose" {
    const ally = std.testing.allocator;
    const F = @import("../scalar.zig").Float(f32);
    const M = MatrixType(F);

    var a = try M.zero(3, 3, ally);
    defer a.deinit();
    // 2 1 3
    // 6 0 5
    // 8 9 7

    try a.set(0, 1, F.from(1,1));
    try a.set(0, 0, F.from(2,1));
    try a.set(0, 2, F.from(3,1));
    try a.set(1, 1, F.zero);
    try a.set(1, 2, F.from(5,1));
    try a.set(1, 0, F.from(6,1));
    try a.set(2, 2, F.from(7,1));
    try a.set(2, 0, F.from(8,1));
    try a.set(2, 1, F.from(9,1));

    var b = try a.transpose(ally);
    defer b.deinit();

    try testing.expect(b.at(0, 0).cmp(.equal, F.from(2,1)));
    try testing.expect(b.at(1, 0).cmp(.equal, F.from(1,1)));
    try testing.expect(b.at(2, 0).cmp(.equal, F.from(3,1)));
    try testing.expect(b.at(0, 1).cmp(.equal, F.from(6,1)));
    try testing.expect(b.at(1, 1).cmp(.equal, F.from(0,1)));
    try testing.expect(b.at(2, 1).cmp(.equal, F.from(5,1)));
    try testing.expect(b.at(0, 2).cmp(.equal, F.from(8,1)));
    try testing.expect(b.at(1, 2).cmp(.equal, F.from(9,1)));
    try testing.expect(b.at(2, 2).cmp(.equal, F.from(7,1)));
}

test "matrix multiplication" {
    const ally = std.testing.allocator;
    const F = @import("../scalar.zig").Float(f32);
    const M = MatrixType(F);

    // 2 1 3   1 0 1   2 1 6
    // 6 0 5 * 0 1 1 = 6 0 11
    // 8 9 7   0 0 1   8 9 24

    var a = try M.zero(3, 3, ally);
    defer a.deinit();

    try a.set(0, 1, F.from(1,1));
    try a.set(0, 0, F.from(2,1));
    try a.set(0, 2, F.from(3,1));
    try a.set(1, 1, F.zero);
    try a.set(1, 2, F.from(5,1));
    try a.set(1, 0, F.from(6,1));
    try a.set(2, 2, F.from(7,1));
    try a.set(2, 0, F.from(8,1));
    try a.set(2, 1, F.from(9,1));

    var b = try M.eye(3, ally);
    defer b.deinit();

    try b.set(0, 2, F.from(1,1));
    try b.set(1, 2, F.from(1,1));

    var c = try a.mul(b, ally);
    defer c.deinit();

    try testing.expect(c.at(0, 0).cmp(.equal, F.from(2,1)));
    try testing.expect(c.at(0, 1).cmp(.equal, F.from(1,1)));
    try testing.expect(c.at(0, 2).cmp(.equal, F.from(6,1)));
    try testing.expect(c.at(1, 0).cmp(.equal, F.from(6,1)));
    try testing.expect(c.at(1, 1).cmp(.equal, F.from(0,1)));
    try testing.expect(c.at(1, 2).cmp(.equal, F.from(11,1)));
    try testing.expect(c.at(2, 0).cmp(.equal, F.from(8,1)));
    try testing.expect(c.at(2, 1).cmp(.equal, F.from(9,1)));
    try testing.expect(c.at(2, 2).cmp(.equal, F.from(24,1)));
}
