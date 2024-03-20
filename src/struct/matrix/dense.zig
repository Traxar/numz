const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;
const Allocator = std.mem.Allocator;
const PermutationQueue = @import("../utils/permutationQueue.zig").PermutationQueueType;

/// struct for matrix operations
/// Element must be a field
pub fn MatrixType(comptime Element: type) type {
    const Vector = @import("../vector.zig").VectorType(Element);
    const SIMDElement = Element.SIMDType(null);
    return struct {
        const Matrix = @This();

        rows: usize,
        cols: usize,
        val: []SIMDElement,

        /// deinitialize matrix
        pub fn deinit(a: Matrix, allocator: Allocator) void {
            allocator.free(a.val);
        }

        ///allocate empty matrix
        pub fn init(rows: usize, cols: usize, allocator: Allocator) !Matrix {
            const n_SIMD = 1 + @divFloor(cols - 1, SIMDElement.SIMDsize); //divCeil
            const res = Matrix{
                .val = try allocator.alloc(SIMDElement, n_SIMD * rows),
                .rows = rows,
                .cols = cols,
            };
            for (0..rows) |i| {
                res.rowCast(i).SIMDsetTail(SIMDElement.zero);
            }
            return res;
        }

        ///allocate empty matrix with same dimensions as a
        pub fn like(a: Matrix, allocator: Allocator) !Matrix {
            return init(a.rows, a.cols, allocator);
        }

        pub fn fill(a: Matrix, b: Element) void {
            for (0..a.rows) |i| {
                a.rowCast(i).fill(b);
            }
        }

        /// res <- a
        pub fn copy(a: Matrix, res: Matrix) void {
            assert(a.rows == res.rows);
            assert(a.cols == res.cols);
            if (a.val.ptr != res.val.ptr) {
                @memcpy(res.val, a.val);
            }
        }

        fn rowCast(a: Matrix, row: usize) Vector {
            assert(row < a.rows);
            const n_SIMD = @divExact(a.val.len, a.rows);
            const start = n_SIMD * row;
            return Vector{ .val = a.val[start .. start + n_SIMD], .len = a.cols };
        }

        /// return element at row and column
        pub fn at(a: Matrix, row: usize, col: usize) Element {
            assert(col < a.cols); //matrix dimensions
            return a.rowCast(row).at(col);
        }

        /// sets element at row and column
        pub fn set(a: Matrix, row: usize, col: usize, b: Element) void {
            assert(col < a.cols); //matrix dimensions
            a.rowCast(row).set(col, b);
        }

        /// res <- A^T
        /// O(m*n)
        pub fn transpose(a: Matrix, res: Matrix) void {
            assert(a.cols == res.rows);
            assert(a.rows == res.cols);
            assert(a.val.ptr != res.val.ptr); //no inplace
            for (0..res.rows) |i| {
                for (0..res.cols) |j| {
                    res.set(i, j, a.at(j, i));
                }
            }
        }

        /// res <- b * a
        pub fn mulE(a: Matrix, b: Element, res: Matrix) void {
            assert(res.rows == a.rows);
            assert(res.cols == a.cols);
            if (b.cmp(.eq, Element.zero)) {
                const a_ = Vector{ .val = a.val, .len = a.val.len * SIMDElement.SIMDsize };
                a_.fill(Element.zero);
            } else {
                for (0..res.rows) |i| {
                    const a_ = a.rowCast(i);
                    const res_ = res.rowCast(i);
                    a_.mulE(b, res_);
                }
            }
        }

        /// res <- a / b
        pub fn divE(a: Matrix, b: Element, res: Matrix) void {
            assert(res.rows == a.rows);
            assert(res.cols == a.cols);
            assert(b.cmp(.neq, Element.zero));
            for (0..res.rows) |i| {
                const a_ = a.rowCast(i);
                const res_ = res.rowCast(i);
                a_.divE(b, res_);
            }
        }

        /// res <- a * b
        pub fn mulV(a: Matrix, b: Vector, res: Vector) void {
            assert(a.cols == b.len);
            assert(a.rows == res.len);
            assert(b.val.ptr != res.val.ptr);
            for (0..a.rows) |i| {
                res.set(i, a.rowCast(i).dot(b));
            }
        }

        /// res <- a^T * b
        pub fn mulTV(a: Matrix, b: Vector, res: Vector) void {
            assert(a.cols == b.len);
            assert(a.rows == res.len);
            assert(b.val.ptr != res.val.ptr);
            a.rowCast(0).mulE(b.at(0), res);
            for (1..a.rows) |i| {
                res.mulEAdd(b.at(i), a.rowCast(i), res);
            }
        }

        /// res <- a * b
        pub fn mul(a: Matrix, b: Matrix, res: Matrix) void {
            assert(a.cols == b.rows);
            assert(res.rows == a.rows);
            assert(res.cols == b.cols);
            assert(a.val.ptr != res.val.ptr);
            assert(b.val.ptr != res.val.ptr);
            for (0..res.rows) |i| {
                b.mulTV(a.rowCast(i), res.rowCast(i));
            }
        }

        /// res <- a + b
        pub fn add(a: Matrix, b: Matrix, res: Matrix) void {
            assert(a.rows == b.rows and a.rows == res.rows);
            assert(a.cols == b.cols and a.cols == res.cols);
            for (0..res.rows) |i| {
                a.rowCast(i).add(b.rowCast(i), res.rowCast(i));
            }
        }

        /// res <- a - b
        pub fn sub(a: Matrix, b: Matrix, res: Matrix) void {
            assert(a.rows == b.rows and a.rows == res.rows);
            assert(a.cols == b.cols and a.cols == res.cols);
            for (0..res.rows) |i| {
                a.rowCast(i).sub(b.rowCast(i), res.rowCast(i));
            }
        }

        /// frobenius norm = sqrt of sum of all squared entries
        pub fn frobenius(a: Matrix) Element {
            var sum = Element.zero;
            for (0..a.rows) |i| {
                const row = a.rowCast(i);
                sum = sum.add(row.dot(row));
            }
            return sum.sqrt();
        }

        pub const LQ = struct {
            l: Matrix,
            q: Matrix,

            pub fn init(rows: usize, cols: usize, allocator: Allocator) !LQ {
                var res: LQ = undefined;
                res.l = try Matrix.init(rows, cols, allocator);
                errdefer res.l.deinit(allocator);
                res.q = try Matrix.init(cols, cols, allocator);
                errdefer res.q.deinit(allocator);
                return res;
            }

            pub fn deinit(a: LQ, allocator: Allocator) void {
                a.l.deinit(allocator);
                a.q.deinit(allocator);
            }

            pub fn det(a: LQ) Element {
                var res = Element.eye;
                for (0..a.l.rows) |i| {
                    res = res.mul(a.l.at(i, i));
                }
                return res;
            }

            ///modified Gram-Schmidt Solver
            pub fn fromMGS(lq: LQ, a: Matrix, eps: Element) !void {
                a.copy(lq.q);
                lq.l.fill(Element.zero);
                for (0..lq.q.rows) |k| {
                    const row_k = lq.q.rowCast(k);
                    const norm = row_k.norm();
                    if (norm.cmp(.lt, eps)) return error.NearlySingular;
                    lq.l.set(k, k, norm);
                    row_k.divE(norm, row_k);
                    for (k + 1..lq.q.rows) |j| {
                        const row_j = lq.q.rowCast(j);
                        const dot = row_k.dot(row_j);
                        lq.l.set(j, k, dot);
                        row_j.mulEAdd(dot.neg(), row_k, row_j);
                    }
                }
            }
        };
    };
}

test "matrix creation" {
    const n = 5;
    const m = 8;

    const ally = testing.allocator;
    const F = @import("../field.zig").Float(f32);
    const M = MatrixType(F);

    const a = try M.init(n, m, ally); //const only refers to the dimensions
    defer a.deinit(ally);
    for (0..n) |i| {
        for (0..m) |j| {
            a.set(i, j, F.from(@intCast(i * m + j), 1));
        }
    }
    for (0..n) |i| {
        for (0..m) |j| {
            try testing.expectEqual(F.from(@intCast(i * m + j), 1), a.at(i, j));
        }
    }
}

test "matrix transpose" {
    const ally = std.testing.allocator;
    const F = @import("../field.zig").Float(f32);
    const M = MatrixType(F);

    const a = try M.init(3, 4, ally);
    defer a.deinit(ally);
    var c: isize = 1;
    for (0..3) |i| {
        for (0..4) |j| {
            if (j >= i) {
                a.set(i, j, F.from(c, 1));
                c += 1;
            }
        }
    }

    const b = try M.init(4, 3, ally);
    defer b.deinit(ally);
    a.transpose(b);
    c = 1;
    for (0..3) |i| {
        for (0..4) |j| {
            if (j >= i) {
                try testing.expectEqual(F.from(c, 1), b.at(j, i));
                c += 1;
            }
        }
    }
}

test "matrix addition" {
    const ally = std.testing.allocator;
    const F = @import("../field.zig").Float(f32);
    const M = MatrixType(F);

    // 1 2 3   1 0 0   0 2 3
    // 0 0 4 - 2 0 0 =-2 0 4
    // 0 0 0   3 4 0  -3-4 0

    const n = 3;
    const a = try M.init(n, n, ally);
    defer a.deinit(ally);
    a.fill(F.zero);
    a.set(0, 0, F.from(1, 1));
    a.set(0, 1, F.from(2, 1));
    a.set(0, 2, F.from(3, 1));
    a.set(1, 2, F.from(4, 1));

    const b = try M.init(n, n, ally);
    defer b.deinit(ally);
    a.transpose(b);
    b.mulE(F.eye.neg(), b);

    var c = try M.init(n, n, ally);
    defer c.deinit(ally);
    a.add(b, c);

    try testing.expectEqual(F.from(0, 1), c.at(0, 0));
    try testing.expectEqual(F.from(2, 1), c.at(0, 1));
    try testing.expectEqual(F.from(3, 1), c.at(0, 2));
    try testing.expectEqual(F.from(-2, 1), c.at(1, 0));
    try testing.expectEqual(F.from(0, 1), c.at(1, 1));
    try testing.expectEqual(F.from(4, 1), c.at(1, 2));
    try testing.expectEqual(F.from(-3, 1), c.at(2, 0));
    try testing.expectEqual(F.from(-4, 1), c.at(2, 1));
    try testing.expectEqual(F.from(0, 1), c.at(2, 2));
}

test "matrix multiplication" {
    const ally = std.testing.allocator;
    const F = @import("../field.zig").Float(f32);
    const M = MatrixType(F);

    // 1 2 3   1 0-1   1 5 0
    // 0 4 5 * 0 1-1 = 0 9 1
    // 0 0 6   0 1 1   0 6 6

    const n = 3;
    const a = try M.init(n, n, ally);
    defer a.deinit(ally);
    a.fill(F.zero);
    a.set(0, 0, F.from(1, 1));
    a.set(0, 1, F.from(2, 1));
    a.set(0, 2, F.from(3, 1));
    a.set(1, 1, F.from(4, 1));
    a.set(1, 2, F.from(5, 1));
    a.set(2, 2, F.from(6, 1));

    const b = try M.init(n, n, ally);
    defer b.deinit(ally);
    b.fill(F.zero);
    b.set(0, 0, F.from(1, 1));
    b.set(0, 2, F.from(-1, 1));
    b.set(1, 1, F.from(1, 1));
    b.set(1, 2, F.from(-1, 1));
    b.set(2, 1, F.from(1, 1));
    b.set(2, 2, F.from(1, 1));

    var c = try M.init(n, n, ally);
    defer c.deinit(ally);
    a.mul(b, c);

    try testing.expectEqual(F.from(1, 1), c.at(0, 0));
    try testing.expectEqual(F.from(5, 1), c.at(0, 1));
    try testing.expectEqual(F.from(0, 1), c.at(0, 2));
    try testing.expectEqual(F.from(0, 1), c.at(1, 0));
    try testing.expectEqual(F.from(9, 1), c.at(1, 1));
    try testing.expectEqual(F.from(1, 1), c.at(1, 2));
    try testing.expectEqual(F.from(0, 1), c.at(2, 0));
    try testing.expectEqual(F.from(6, 1), c.at(2, 1));
    try testing.expectEqual(F.from(6, 1), c.at(2, 2));
}

test "matrix multiplication with element" {
    const ally = std.testing.allocator;
    const F = @import("../field.zig").Float(f32);
    const M = MatrixType(F);

    const n = 3;
    const a = try M.init(n, n, ally);
    defer a.deinit(ally);
    a.fill(F.zero);
    for (0..n) |i| {
        a.set(i, i, F.eye);
    }
    const b = F.from(-314, 100);

    a.mulE(b, a);
    for (0..n) |i| {
        for (0..n) |j| {
            try testing.expect(a.at(i, j).cmp(.eq, if (i == j) b else F.zero));
        }
    }

    a.divE(b, a);
    for (0..n) |i| {
        for (0..n) |j| {
            try testing.expect(a.at(i, j).cmp(.eq, if (i == j) F.eye else F.zero));
        }
    }
}

test "matrix mulpiplication with vector" {
    const ally = std.testing.allocator;
    const F = @import("../field.zig").Float(f32);
    const V = @import("../vector.zig").VectorType(F);
    const M = MatrixType(F);

    // 1 2 3   -2   -1
    // 0 4 5 * -1 =  1
    // 0 0 6    1    6

    const n = 3;
    const a = try M.init(n, n, ally);
    defer a.deinit(ally);
    a.fill(F.zero);
    a.set(0, 0, F.from(1, 1));
    a.set(0, 1, F.from(2, 1));
    a.set(0, 2, F.from(3, 1));
    a.set(1, 1, F.from(4, 1));
    a.set(1, 2, F.from(5, 1));
    a.set(2, 2, F.from(6, 1));

    const b = try V.init(n, ally);
    defer b.deinit(ally);
    b.set(0, F.from(-2, 1));
    b.set(1, F.from(-1, 1));
    b.set(2, F.from(1, 1));

    const c = try V.init(n, ally);
    defer c.deinit(ally);
    a.mulV(b, c);

    try testing.expect(c.at(0).cmp(.eq, F.from(-1, 1)));
    try testing.expect(c.at(1).cmp(.eq, F.from(1, 1)));
    try testing.expect(c.at(2).cmp(.eq, F.from(6, 1)));
}

test "matrix norm" {
    const ally = std.testing.allocator;
    const F = @import("../field.zig").Float(f32);
    const M = MatrixType(F);

    // 1  2  3
    // 2  1  2
    // 0  1  1

    const n = 3;
    const a = try M.init(n, n, ally);
    defer a.deinit(ally);
    a.fill(F.zero);
    a.set(0, 0, F.from(1, 1));
    a.set(0, 1, F.from(2, 1));
    a.set(0, 2, F.from(3, 1));
    a.set(1, 0, F.from(2, 1));
    a.set(1, 1, F.from(1, 1));
    a.set(1, 2, F.from(2, 1));
    a.set(2, 1, F.from(1, 1));
    a.set(2, 2, F.from(1, 1));

    try testing.expectEqual(F.from(5, 1), a.frobenius());
}

test "matrix LQ" {
    const ally = std.testing.allocator;
    const F = @import("../field.zig").Float(f32);
    const M = MatrixType(F);

    // 1  2  3
    // 2  2  2
    // 0  1  1

    const n = 3;
    const a = try M.init(n, n, ally);
    defer a.deinit(ally);
    a.fill(F.zero);
    a.set(0, 0, F.from(1, 1));
    a.set(0, 1, F.from(2, 1));
    a.set(0, 2, F.from(3, 1));
    a.set(1, 0, F.from(2, 1));
    a.set(1, 1, F.from(2, 1));
    a.set(1, 2, F.from(2, 1));
    a.set(2, 1, F.from(1, 1));
    a.set(2, 2, F.from(1, 1));

    const eps = F.from(1, 1E5);

    const lq = try M.LQ.init(n, n, ally);
    defer lq.deinit(ally);
    try lq.fromMGS(a, eps);

    const b = try a.like(ally);
    defer b.deinit(ally);
    lq.l.mul(lq.q, b);

    for (0..n) |i| {
        for (0..n) |j| {
            try testing.expectEqual(a.at(i, j), b.at(i, j));
        }
    }

    const c = try lq.q.like(ally);
    defer c.deinit(ally);
    lq.q.transpose(c);

    const d = try lq.q.like(ally);
    defer d.deinit(ally);
    lq.q.mul(c, d);

    for (0..n) |i| {
        d.set(i, i, d.at(i, i).sub(F.eye));
    }

    try testing.expect(d.frobenius().cmp(.lt, eps));
}
