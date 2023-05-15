const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;

/// return struct that can allocate and free matrices
/// operations by this struct result in a newly allocated matrix
/// operation by the matrix struct are inplace
pub fn MatrixType(comptime Scalar: type) type {
    return struct {
        const Self = @This();
        const MultiArrayList = std.MultiArrayList(struct { col: usize, val: Scalar });
        allocator: std.mem.Allocator,

        fn init(self: Self, rows: usize, cols: usize) !Matrix {
            return Matrix{
                .rows = try self.allocator.alloc(MultiArrayList, rows),
                .cols = cols,
                .allocator = self.allocator,
            };
        }

        /// deinitialize matrix
        pub fn deinit(self: Self, a: Matrix) void {
            for (0..a.rows.len) |i| {
                a.rows[i].deinit(self.allocator);
            }
            self.allocator.free(a.rows);
        }

        /// create an empty matrix with n rows and m columns
        /// do not forget to call deinit() on result
        pub fn zero(self: Self, rows: usize, cols: usize) !Matrix {
            var res = try self.init(rows, cols);
            for (0..rows) |i| {
                res.rows[i] = MultiArrayList{};
            }
            return res;
        }

        /// create the identity matrix of size n
        /// do not forget to call deinit() on result
        pub fn eye(self: Self, size: usize) !Matrix {
            var res = try self.zero(size, size);
            for (0..size) |i| {
                try res.set(i, i, Scalar.eye);
            }
            return res;
        }

        /// allocates a copy of matrix a
        pub fn copy(self: Self, a: Matrix) !Matrix {
            var res = try self.init(a.rows.len, a.cols);
            for (0..a.rows.len) |i| {
                res.rows[i] = a.rows[i].clone(self.allocator);
            }
            return res;
        }

        /// sparse matrix with runtime size
        /// all functions in this struct are inplace
        /// they might require memory allocation since this is a sparse datastructure
        const Matrix = struct {
            rows: []MultiArrayList,
            cols: usize,
            allocator: std.mem.Allocator,

            /// return index of col in row
            /// performs binary search
            /// O(log(m))
            fn findIndex(a: Matrix, row: usize, col: usize) struct { index: usize, exists: bool } {
                const r = a.rows[row].items(.col);
                var max = r.len;
                if (max == 0) {
                    return .{ .index = 0, .exists = false };
                }
                var min: usize = 0;
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

            /// return element at row and column
            /// O(log(m))
            pub fn at(a: Matrix, row: usize, col: usize) Scalar {
                assert(row < a.rows.len);
                assert(col < a.cols);
                const i = a.findIndex(row, col);
                if (i.exists) {
                    return a.rows[row].items(.val)[i.index];
                } else {
                    return Scalar.zero;
                }
            }

            // set element at row i and column j to b
            // O(m)
            pub fn set(a: Matrix, row: usize, col: usize, b: Scalar) !void {
                assert(row < a.rows.len);
                assert(col < a.cols);
                const i = a.findIndex(row, col);
                if (b.cmp(.equal, Scalar.zero)) {
                    if (i.exists) {
                        a.rows[row].orderedRemove(i.index);
                    } // else do nothing
                } else {
                    if (i.exists) {
                        a.rows[row].items(.val)[i.index] = b;
                    } else {
                        try a.rows[row].insert(a.allocator, i.index, .{ .col = col, .val = b });
                    }
                }
                return;
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
    var a = try M.zero(n, m);
    defer M.deinit(a);
    for (0..n) |i| {
        for (0..m) |j| {
            try testing.expect(a.at(i, j).cmp(.equal, F.zero));
        }
    }

    var b = try M.eye(n);
    defer M.deinit(b);
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
    const M = MatrixType(F){ .allocator = allocator };

    var a = try M.zero(3, 3);
    defer M.deinit(a);
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

    try testing.expectEqual(@as(usize, 3), a.rows[0].len);
    try testing.expectEqual(@as(usize, 2), a.rows[1].len);

    try a.set(2, 2, F.from(0));
    try a.set(2, 0, F.from(0));
    try a.set(2, 1, F.from(0));

    try testing.expectEqual(@as(usize, 0), a.rows[2].len);
}
