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

/// return struct that can allocate and free matrices
/// operations by this struct result in a newly allocated matrix
pub fn MatrixType(comptime Scalar: type) type {
    return struct {
        const Self = @This();
        const VectorType = @import("vector.zig").VectorType(Scalar);
        const Vector = VectorType.Vector;
        allocptr: *const std.mem.Allocator,

        fn init(allocptr: *const std.mem.Allocator, rows: usize, cols: usize) !Matrix {
            return Matrix{
                .val = (try allocptr.alloc(Matrix.Row, rows)).ptr,
                .rows = rows,
                .cols = cols,
                .allocptr = allocptr,
            };
        }

        /// create an empty matrix with n rows and m columns
        /// do not forget to call deinit() on result
        pub fn zero(self: Self, rows: usize, cols: usize) !Matrix {
            var res = try init(self.allocptr, rows, cols);
            for (0..res.rows) |i| {
                res.val[i] = Matrix.Row{};
            }
            return res;
        }

        /// create the identity matrix of size n
        /// do not forget to call deinit() on result
        pub fn eye(self: Self, size: usize) !Matrix {
            var res = try self.zero(size, size);
            for (0..res.rows) |i| {
                try res.set(i, i, Scalar.eye);
            }
            return res;
        }

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

            /// allocates a copy of matrix a
            pub fn copy(a: Matrix) !Matrix {
                var res = try init(a.allocptr, a.rows, a.cols);
                for (0..a.rows) |i| {
                    res.val[i] = a.val[i].clone(a.allocptr.*);
                }
                return res;
            }

            /// return index of col in row
            /// performs binary search
            /// O(log(m)), O(1) if dense
            fn indAt(a: Matrix, row: usize, col: usize) struct { index: usize, exists: bool } {
                const r = a.val[row].items(.col);
                var min = r.len -| (a.cols - col);
                var max = @min(r.len, col);
                while (min != max) {
                    const pivot = @divFloor(min + max, 2);
                    if (col <= r[pivot]) {
                        max = pivot;
                    } else {
                        min = pivot + 1;
                    }
                }
                return (.{ .index = min, .exists = (min < r.len and col == r[min]) });
            }

            /// return number of entries in the row
            fn entAt(a: Matrix, row: usize) usize {
                return a.val[row].len;
            }

            /// return element at row and columnIndex
            /// O(1)
            fn valAt(a: Matrix, row: usize, colIndex: usize) Scalar {
                assert(row < a.rows);
                assert(colIndex < a.val[row].len);
                return a.val[row].items(.val)[colIndex];
            }

            /// return column at row and columnIndex
            /// O(1)
            fn colAt(a: Matrix, row: usize, colIndex: usize) usize {
                assert(row < a.rows);
                assert(colIndex < a.val[row].len);
                return a.val[row].items(.col)[colIndex];
            }

            /// set element at row and columnIndex
            /// O(1)
            fn setAt(a: Matrix, row: usize, colIndex: usize, b: Scalar) void {
                assert(row < a.rows);
                assert(colIndex < a.val[row].len);
                a.val[row].items(.val)[colIndex] = b;
            }

            /// return element at row and column
            /// O(log(m))
            pub fn at(a: Matrix, row: usize, col: usize) Scalar {
                assert(row < a.rows);
                assert(col < a.cols);
                const i = a.indAt(row, col);
                if (i.exists) {
                    return a.val[row].items(.val)[i.index];
                } else {
                    return Scalar.zero;
                }
            }

            // set element at row i and column j to b
            // O(m)
            pub fn set(a: Matrix, row: usize, col: usize, b: Scalar) !void {
                assert(row < a.rows);
                assert(col < a.cols);
                const i = a.indAt(row, col);
                if (b.cmp(.equal, Scalar.zero)) {
                    if (i.exists) {
                        a.val[row].orderedRemove(i.index);
                    } // else do nothing
                } else {
                    if (i.exists) {
                        a.setAt(row, i.index, b);
                    } else {
                        try a.val[row].insert(a.allocptr.*, i.index, .{ .col = col, .val = b });
                    }
                }
                return;
            }

            /// a <- a * b
            /// return a
            /// O(n*m)
            pub fn mulS(a: Matrix, b: Scalar) Matrix {
                for (0..a.rows) |i| {
                    for (0..a.entAt(i)) |j| {
                        a.setAt(i, j, a.valAt(i, j).mul(b));
                    }
                }
                return a;
            }

            /// return a * b
            /// O(n*m)
            pub fn mulV(a: Matrix, b: Vector) !Vector {
                assert(a.cols == b.len);
                var res = try (VectorType{ .allocptr = a.allocptr }).rep(Scalar.zero, a.rows);
                for (0..a.rows) |i| {
                    for (0..a.entAt(i)) |j| {
                        res.set(i, res.at(i).add(a.valAt(i, j).mul(b.at(a.colAt(i, j)))));
                    }
                }
                return res;
            }

            /// a <- a * b
            /// return a
            pub fn divS(a: Matrix, b: Scalar) Matrix {
                for (0..a.rows) |i| {
                    for (0..a.entAt(i)) |j| {
                        a.setAt(i, j, a.valAt(i, j).div(b));
                    }
                }
                return a;
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
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();
    defer {
        const deinit_status = gpa.deinit();
        //fail test; can't try in defer as defer is executed after we return
        if (deinit_status == .leak) testing.expect(false) catch @panic("TEST FAIL");
    }

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
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();
    defer {
        const deinit_status = gpa.deinit();
        //fail test; can't try in defer as defer is executed after we return
        if (deinit_status == .leak) testing.expect(false) catch @panic("TEST FAIL");
    }

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
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();
    defer {
        const deinit_status = gpa.deinit();
        //fail test; can't try in defer as defer is executed after we return
        if (deinit_status == .leak) testing.expect(false) catch @panic("TEST FAIL");
    }

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
