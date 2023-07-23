const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;
const Allocator = std.mem.Allocator;
const PPQ = @import("utils/permutationPriorityQueue.zig").PermutationPriorityQueue;
const Permutation = @import("permutation.zig").Permutation;

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
                res.val[i] = try a.val[i].clone(allocator);
            }
            return res;
        }

        /// b := a^T
        /// return b
        /// O(m*n)
        pub fn transpose(a: Matrix, allocator: Allocator) !Matrix {
            var res = try zero(a.cols, a.rows, allocator);
            var nz_col = try allocator.alloc(usize, a.cols);
            defer allocator.free(nz_col);
            for (0..a.cols) |j| {
                nz_col[j] = 0;
            }
            for (0..a.rows) |i| {
                for (0..a.entAt(i)) |j| {
                    nz_col[a.colAt(i, j)] += 1;
                }
            }
            for (0..a.cols) |j| {
                try res.val[j].ensureTotalCapacity(allocator, nz_col[j]);
            }
            for (0..a.rows) |i| {
                for (0..a.entAt(i)) |j| {
                    res.val[a.colAt(i, j)].appendAssumeCapacity(.{ .col = i, .val = a.valAt(i, j, true) });
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
        fn valAt(a: Matrix, row: usize, ind: usize, exists: bool) Scalar {
            if (!exists) {
                return Scalar.zero;
            } else {
                return a.val[row].items(.val)[ind];
            }
        }

        /// return column at row and columnIndex i
        /// O(1)
        fn colAt(a: Matrix, row: usize, ind: usize) usize {
            return a.val[row].items(.col)[ind];
        }

        /// set element at row and columnIndex i
        /// O(1)
        fn setAt(a: Matrix, row: usize, ind: usize, exists: bool, col: usize, b: Scalar) !void {
            if (b.cmp(.eq, Scalar.zero)) {
                if (exists) {
                    a.val[row].orderedRemove(ind);
                } // else do nothing
            } else {
                if (exists) {
                    a.val[row].items(.val)[ind] = b;
                } else {
                    try a.val[row].insert(a.allocator, ind, .{ .col = col, .val = b });
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

        /// c := a + b
        /// return c
        /// O(n*m)
        pub fn add(a: Matrix, b: Matrix, allocator: Allocator) !Matrix {
            assert(a.rows == b.rows);
            assert(a.cols == b.cols);
            var res = try zero(a.rows, a.cols, allocator);
            for (0..a.rows) |i| {
                try res.val[i].ensureTotalCapacity(allocator, a.entAt(i) + b.entAt(i));
                var j_a: usize = 0;
                var j_b: usize = 0;
                while (true) {
                    const col_a = if (j_a < a.entAt(i)) a.colAt(i, j_a) else a.cols;
                    const col_b = if (j_b < b.entAt(i)) b.colAt(i, j_b) else b.cols;
                    if (col_a == col_b) {
                        if (col_a == a.cols) break;
                        const val = a.valAt(i, j_a, true).add(b.valAt(i, j_b, true));
                        if (val.cmp(.neq, Scalar.zero)) {
                            res.val[i].appendAssumeCapacity(.{
                                .col = col_a,
                                .val = val,
                            });
                        }
                        j_a += 1;
                        j_b += 1;
                    } else if (col_a < col_b) {
                        res.val[i].appendAssumeCapacity(.{ .col = col_a, .val = a.valAt(i, j_a, true) });
                        j_a += 1;
                    } else { //col_b < col_a
                        res.val[i].appendAssumeCapacity(.{ .col = col_b, .val = b.valAt(i, j_b, true) });
                        j_b += 1;
                    }
                }
                res.val[i].shrinkAndFree(allocator, res.entAt(i));
            }
            return res;
        }

        /// c := a + b
        /// return c
        /// O(n*m)
        pub fn sub(a: Matrix, b: Matrix, allocator: Allocator) !Matrix {
            assert(a.rows == b.rows);
            assert(a.cols == b.cols);
            var res = try zero(a.rows, a.cols, allocator);
            for (0..a.rows) |i| {
                try res.val[i].ensureTotalCapacity(allocator, a.entAt(i) + b.entAt(i));
                var j_a: usize = 0;
                var j_b: usize = 0;
                while (true) {
                    const col_a = if (j_a < a.entAt(i)) a.colAt(i, j_a) else a.cols;
                    const col_b = if (j_b < b.entAt(i)) b.colAt(i, j_b) else b.cols;
                    if (col_a == col_b) {
                        if (col_a == a.cols) break;
                        const val = a.valAt(i, j_a, true).sub(b.valAt(i, j_b, true));
                        if (val.cmp(.neq, Scalar.zero)) {
                            res.val[i].appendAssumeCapacity(.{
                                .col = col_a,
                                .val = val,
                            });
                        }
                        j_a += 1;
                        j_b += 1;
                    } else if (col_a < col_b) {
                        res.val[i].appendAssumeCapacity(.{ .col = col_a, .val = a.valAt(i, j_a, true) });
                        j_a += 1;
                    } else { //col_b < col_a
                        res.val[i].appendAssumeCapacity(.{ .col = col_b, .val = b.valAt(i, j_b, true).neg() });
                        j_b += 1;
                    }
                }
                res.val[i].shrinkAndFree(allocator, res.entAt(i));
            }
            return res;
        }

        fn elimAt(a: Matrix, row_src: usize, row_trg: usize, ind_pvt: usize, nz_col: []IndexSet) !Scalar {
            //allocate result
            var res = Row{};
            try res.ensureTotalCapacity(a.allocator, a.entAt(row_src) + a.entAt(row_trg) - 1);

            //get factor for elimination
            const col_pvt = a.colAt(row_src, ind_pvt);
            const factor = a.at(row_trg, col_pvt).div(a.valAt(row_src, ind_pvt, true)).neg();

            //main loop
            var i_src: usize = 0;
            var i_trg: usize = 0;
            while (true) {
                const col_src = if (i_src < a.entAt(row_src)) a.colAt(row_src, i_src) else a.cols;
                const col_trg = if (i_trg < a.entAt(row_trg)) a.colAt(row_trg, i_trg) else a.cols;
                if (col_src == col_trg) {
                    if (col_src == a.cols) break; //end
                    if (col_src != col_pvt) { //skip at pivot column
                        const val = a.valAt(row_trg, i_trg, true).add(factor.mul(a.valAt(row_src, i_src, true)));
                        if (val.cmp(.neq, Scalar.zero)) {
                            res.appendAssumeCapacity(.{
                                .col = col_src,
                                .val = val,
                            });
                        } else {
                            _ = nz_col[col_trg].swapRemove(row_trg);
                        }
                    }
                    i_src += 1;
                    i_trg += 1;
                } else if (col_src < col_trg) {
                    res.appendAssumeCapacity(.{ .col = col_src, .val = factor.mul(a.valAt(row_src, i_src, true)) });
                    try nz_col[col_trg].put(a.allocator, row_trg, undefined);
                    i_src += 1;
                } else { //col_b < col_a
                    res.appendAssumeCapacity(.{ .col = col_trg, .val = a.valAt(row_trg, i_trg, true) });
                    i_trg += 1;
                }
            }
            res.shrinkAndFree(a.allocator, res.len);

            //replace row
            a.val[row_trg].deinit(a.allocator);
            a.val[row_trg] = res;
            return factor.neg();
        }

        const LU = struct {
            //PAQ=LU
            p: [*]usize,
            q: [*]usize,
            lt: Matrix,
            u: Matrix,

            pub fn deinit(lu: LU, allocator: Allocator) void {
                const n = lu.u.rows;
                allocator.free(lu.p[0..n]);
                allocator.free(lu.q[0..n]);
                lu.lt.deinit();
                lu.u.deinit();
            }

            pub fn det(lu: LU) Scalar {
                var res = Scalar.eye;
                const n = lu.u.rows;
                for (0..n) |i| {
                    res = res.mul(lu.u.at(lu.p[i], lu.q[i]));
                }
                return res;
            }

            pub fn solve(lu: LU, b: Vector, allocator: Allocator) !Vector {
                const n = b.len;

                //apply l-1
                var c = try b.copy(allocator);
                defer c.deinit(allocator);
                for (0..n) |i| {
                    const row = lu.p[i];
                    for (0..lu.lt.entAt(row)) |j| {
                        const col = lu.lt.colAt(row, j);
                        c.set(col, c.at(col).sub(lu.lt.valAt(row, j, true).mul(c.at(row))));
                    }
                }

                var x = try Vector.init(n, allocator);
                //apply u-1
                for (0..n) |i| {
                    const row = lu.p[n - 1 - i];
                    const col_diag = lu.q[n - 1 - i];
                    var diag: Scalar = undefined;
                    var sum = c.at(row);
                    for (0..lu.u.entAt(row)) |j| {
                        const col = lu.u.colAt(row, j);
                        if (col == col_diag) {
                            diag = lu.u.valAt(row, j, true);
                        } else {
                            sum = sum.sub(lu.u.valAt(row, j, true).mul(x.at(col)));
                        }
                    }
                    x.set(col_diag, sum.div(diag));
                }
                return x;
            }
        };

        fn debugprint(self: Matrix) void {
            std.debug.print("\n", .{});
            for (0..self.rows) |i| {
                for (0..self.cols) |j| {
                    std.debug.print("{}, ", .{self.at(i, j).f});
                }
                std.debug.print("\n", .{});
            }
        }

        const IndexSet = std.AutoArrayHashMapUnmanaged(usize, void);

        const SortContext = struct {
            const Self = @This();
            val: IndexSet.DataList,
            pub fn lessThan(ctx: Self, a_index: usize, b_index: usize) bool {
                return ctx.val.items(.key)[a_index] < ctx.val.items(.key)[b_index];
            }
        };

        const RowPriority = struct {
            const Self = @This();
            nonzeros: usize,
            norm1: Scalar,

            pub fn from(row: Row) Self {
                var row_priority = RowPriority{
                    .nonzeros = row.len,
                    .norm1 = Scalar.zero,
                };
                for (0..row.len) |j| {
                    row_priority.norm1 = row_priority.norm1.add(row.items(.val)[j]);
                }
                return row_priority;
            }

            pub fn before(self: Self, other: Self) bool {
                if (self.nonzeros < other.nonzeros) return true;
                if (self.nonzeros == other.nonzeros and self.norm1.cmp(.gt, other.norm1)) return true;
                return false;
            }
        };

        /// O(n*m*(m+log(n)))
        pub fn decompLU(self: Matrix, allocator: Allocator) !LU {
            assert(self.rows == self.cols);
            const n = self.rows;

            //allocating result structs
            var u = try self.copy(allocator);
            errdefer u.deinit();
            var lt = try Matrix.zero(n, n, allocator); //l transpose for faster fill
            errdefer lt.deinit();
            var q = try allocator.alloc(usize, n);
            errdefer allocator.free(q);

            // count nonzeros for each column O(n*m)
            var nz_col = try allocator.alloc(IndexSet, n);
            defer {
                for (0..n) |i| {
                    nz_col[i].deinit(allocator);
                }
                allocator.free(nz_col);
            }
            for (0..n) |i| {
                nz_col[i] = IndexSet{};
            }
            for (0..n) |i| { //n
                for (0..u.entAt(i)) |j| { //m
                    try nz_col[u.colAt(i, j)].put(allocator, i, undefined); //O(1)
                }
            }

            //initialize the priority queue for row pivoting
            var row_queue = try PPQ(RowPriority, RowPriority.before).init(n, allocator);
            errdefer row_queue.deinit(allocator);
            for (0..n) |i| { //n
                row_queue.set(i, RowPriority.from(u.val[i])); //O(log(n)+m)
            }

            //main loop
            for (0..n) |i| { //n
                //get pivot row from queue
                const row_pvt = row_queue.remove(); //O(log(n))
                if (u.entAt(row_pvt) == 0) return error.NonInvertible;

                //find pivot column O(m)
                var ind_pvt: usize = undefined;
                {
                    var nz_pvt: usize = n;
                    var abs_pvt = Scalar.zero;
                    for (0..u.entAt(row_pvt)) |j| { //m
                        const col_j = nz_col[u.colAt(row_pvt, j)].count();
                        const abs_j = u.valAt(row_pvt, j, true).abs();
                        if (col_j < nz_pvt or (col_j == nz_pvt and abs_j.cmp(.gt, abs_pvt))) {
                            nz_pvt = col_j;
                            abs_pvt = abs_j;
                            ind_pvt = j;
                        }
                    }
                }
                const col_pvt = u.colAt(row_pvt, ind_pvt);
                q[i] = col_pvt;

                //remove pivot row form nonzeros
                for (0..u.entAt(row_pvt)) |j| {
                    _ = nz_col[u.colAt(row_pvt, j)].swapRemove(row_pvt);
                }

                //elimination
                try lt.val[row_pvt].ensureTotalCapacity(allocator, nz_col[col_pvt].count());
                nz_col[col_pvt].sort(SortContext{ .val = nz_col[col_pvt].entries }); //O(m*log(m)) or O(m^2) ?
                var elim_iter = nz_col[col_pvt].iterator();
                while (elim_iter.next()) |entry| { //m
                    const row_trg = entry.key_ptr.*;
                    const factor = try u.elimAt(row_pvt, row_trg, ind_pvt, nz_col); //O(m)
                    lt.val[row_pvt].appendAssumeCapacity(.{ .col = row_trg, .val = factor });
                    row_queue.set(row_trg, RowPriority.from(u.val[row_trg])); //update priority //O(log(n))
                }
            }
            const p = row_queue.deinitToPermutation(allocator);
            allocator.free(p.inv[0..n]);

            return LU{
                .p = p.val,
                .q = q.ptr,
                .lt = lt,
                .u = u,
            };
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
            try testing.expect(a.at(i, j).cmp(.eq, F.zero));
        }
    }

    var b = try M.eye(n, ally);
    defer b.deinit();
    for (0..n) |i| {
        for (0..n) |j| {
            if (i == j) {
                try testing.expect(b.at(i, j).cmp(.eq, F.eye));
            } else {
                try testing.expect(b.at(i, j).cmp(.eq, F.zero));
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

    try a.set(0, 1, F.from(1, 1));
    try a.set(0, 0, F.from(2, 1));
    try a.set(0, 2, F.from(3, 1));
    try a.set(1, 1, F.zero);
    try a.set(1, 2, F.from(5, 1));
    try a.set(1, 0, F.from(6, 1));
    try a.set(2, 2, F.from(7, 1));
    try a.set(2, 0, F.from(8, 1));
    try a.set(2, 1, F.from(9, 1));

    try testing.expect(a.at(0, 0).cmp(.eq, F.from(2, 1)));
    try testing.expect(a.at(0, 1).cmp(.eq, F.from(1, 1)));
    try testing.expect(a.at(0, 2).cmp(.eq, F.from(3, 1)));
    try testing.expect(a.at(1, 0).cmp(.eq, F.from(6, 1)));
    try testing.expect(a.at(1, 1).cmp(.eq, F.from(0, 1)));
    try testing.expect(a.at(1, 2).cmp(.eq, F.from(5, 1)));
    try testing.expect(a.at(2, 0).cmp(.eq, F.from(8, 1)));
    try testing.expect(a.at(2, 1).cmp(.eq, F.from(9, 1)));
    try testing.expect(a.at(2, 2).cmp(.eq, F.from(7, 1)));

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
    const a_ = F.from(-314, 100);
    var a = (try M.eye(n, ally)).mulS(a_);
    defer a.deinit();
    for (0..n) |i| {
        for (0..n) |j| {
            if (i == j) {
                try testing.expect(a.at(i, j).cmp(.eq, a_));
            } else {
                try testing.expect(a.at(i, j).cmp(.eq, F.zero));
            }
        }
    }
    _ = a.divS(a_);
    for (0..n) |i| {
        for (0..n) |j| {
            if (i == j) {
                try testing.expect(a.at(i, j).cmp(.eq, a_.div(a_)));
            } else {
                try testing.expect(a.at(i, j).cmp(.eq, F.zero));
            }
        }
    }
}

test "mulpiplication with vector" {
    const ally = std.testing.allocator;
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
    v.set(1, F.from(1, 1));
    v.set(2, F.from(2, 1));

    const av = try a.mulV(v, ally);
    defer av.deinit(ally);
    try testing.expect(av.at(0).cmp(.eq, F.from(-2, 1)));
    try testing.expect(av.at(1).cmp(.eq, F.from(1, 1)));
    try testing.expect(av.at(2).cmp(.eq, F.from(1, 1)));
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

    try a.set(0, 1, F.from(1, 1));
    try a.set(0, 0, F.from(2, 1));
    try a.set(0, 2, F.from(3, 1));
    try a.set(1, 1, F.zero);
    try a.set(1, 2, F.from(5, 1));
    try a.set(1, 0, F.from(6, 1));
    try a.set(2, 2, F.from(7, 1));
    try a.set(2, 0, F.from(8, 1));
    try a.set(2, 1, F.from(9, 1));

    var b = try a.transpose(ally);
    defer b.deinit();

    try testing.expect(b.at(0, 0).cmp(.eq, F.from(2, 1)));
    try testing.expect(b.at(1, 0).cmp(.eq, F.from(1, 1)));
    try testing.expect(b.at(2, 0).cmp(.eq, F.from(3, 1)));
    try testing.expect(b.at(0, 1).cmp(.eq, F.from(6, 1)));
    try testing.expect(b.at(1, 1).cmp(.eq, F.from(0, 1)));
    try testing.expect(b.at(2, 1).cmp(.eq, F.from(5, 1)));
    try testing.expect(b.at(0, 2).cmp(.eq, F.from(8, 1)));
    try testing.expect(b.at(1, 2).cmp(.eq, F.from(9, 1)));
    try testing.expect(b.at(2, 2).cmp(.eq, F.from(7, 1)));
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

    try a.set(0, 1, F.from(1, 1));
    try a.set(0, 0, F.from(2, 1));
    try a.set(0, 2, F.from(3, 1));
    try a.set(1, 1, F.zero);
    try a.set(1, 2, F.from(5, 1));
    try a.set(1, 0, F.from(6, 1));
    try a.set(2, 2, F.from(7, 1));
    try a.set(2, 0, F.from(8, 1));
    try a.set(2, 1, F.from(9, 1));

    var b = try M.eye(3, ally);
    defer b.deinit();

    try b.set(0, 2, F.from(1, 1));
    try b.set(1, 2, F.from(1, 1));

    var c = try a.mul(b, ally);
    defer c.deinit();

    try testing.expect(c.at(0, 0).cmp(.eq, F.from(2, 1)));
    try testing.expect(c.at(0, 1).cmp(.eq, F.from(1, 1)));
    try testing.expect(c.at(0, 2).cmp(.eq, F.from(6, 1)));
    try testing.expect(c.at(1, 0).cmp(.eq, F.from(6, 1)));
    try testing.expect(c.at(1, 1).cmp(.eq, F.from(0, 1)));
    try testing.expect(c.at(1, 2).cmp(.eq, F.from(11, 1)));
    try testing.expect(c.at(2, 0).cmp(.eq, F.from(8, 1)));
    try testing.expect(c.at(2, 1).cmp(.eq, F.from(9, 1)));
    try testing.expect(c.at(2, 2).cmp(.eq, F.from(24, 1)));
}

test "matrix addition and subtraction" {
    const ally = std.testing.allocator;
    const F = @import("../scalar.zig").Float(f32);
    const M = MatrixType(F);

    // 2 0 1   1 0 1   3 0 2
    // 0 0 5 + 0 1 1 = 0 1 6
    // 8 9-1   0 0 1   8 9 0

    var a = try M.zero(3, 3, ally);
    defer a.deinit();

    try a.set(0, 0, F.from(2, 1));
    try a.set(0, 1, F.zero);
    try a.set(0, 2, F.from(1, 1));
    try a.set(1, 0, F.zero);
    try a.set(1, 1, F.zero);
    try a.set(1, 2, F.from(5, 1));
    try a.set(2, 0, F.from(8, 1));
    try a.set(2, 1, F.from(9, 1));
    try a.set(2, 2, F.from(-1, 1));

    var b = try M.eye(3, ally);
    defer b.deinit();

    try b.set(0, 2, F.from(1, 1));
    try b.set(1, 2, F.from(1, 1));

    var c = try a.add(b, ally);
    defer c.deinit();

    try testing.expect(c.at(0, 0).cmp(.eq, F.from(3, 1)));
    try testing.expect(c.at(0, 1).cmp(.eq, F.zero));
    try testing.expect(c.at(0, 2).cmp(.eq, F.from(2, 1)));
    try testing.expect(c.at(1, 0).cmp(.eq, F.zero));
    try testing.expect(c.at(1, 1).cmp(.eq, F.from(1, 1)));
    try testing.expect(c.at(1, 2).cmp(.eq, F.from(6, 1)));
    try testing.expect(c.at(2, 0).cmp(.eq, F.from(8, 1)));
    try testing.expect(c.at(2, 1).cmp(.eq, F.from(9, 1)));
    try testing.expect(c.at(2, 2).cmp(.eq, F.zero));

    // 2 0 1   1 0 1   1 0 0
    // 0 0 5 - 0 1 1 = 0-1 4
    // 8 9-1   0 0 1   8 9-2

    var d = try a.sub(b, ally);
    defer d.deinit();

    try testing.expect(d.at(0, 0).cmp(.eq, F.from(1, 1)));
    try testing.expect(d.at(0, 1).cmp(.eq, F.zero));
    try testing.expect(d.at(0, 2).cmp(.eq, F.zero));
    try testing.expect(d.at(1, 0).cmp(.eq, F.zero));
    try testing.expect(d.at(1, 1).cmp(.eq, F.from(-1, 1)));
    try testing.expect(d.at(1, 2).cmp(.eq, F.from(4, 1)));
    try testing.expect(d.at(2, 0).cmp(.eq, F.from(8, 1)));
    try testing.expect(d.at(2, 1).cmp(.eq, F.from(9, 1)));
    try testing.expect(d.at(2, 2).cmp(.eq, F.from(-2, 1)));
}

test "LU solve" {
    const ally = std.testing.allocator;
    const F = @import("../scalar.zig").Float(f64);
    const M = MatrixType(F);

    var a = try M.zero(3, 3, ally);
    defer a.deinit();
    // 2  1 3
    // 6  0 5
    // 8 10 7
    // det = 40+180-100-42 = 78

    try a.set(0, 0, F.from(2, 1));
    try a.set(0, 1, F.from(1, 1));
    try a.set(0, 2, F.from(3, 1));
    try a.set(1, 0, F.from(6, 1));
    try a.set(1, 1, F.zero);
    try a.set(1, 2, F.from(5, 1));
    try a.set(2, 0, F.from(8, 1));
    try a.set(2, 1, F.from(10, 1));
    try a.set(2, 2, F.from(7, 1));

    var b = try M.Vector.rep(F.eye, 3, ally);
    defer b.deinit(ally);
    b.set(1, F.from(2, 1));
    b.set(2, F.from(3, 1));

    //decompose
    const lu = try a.decompLU(ally);
    defer lu.deinit(ally);

    //solve
    var x = try lu.solve(b, ally);
    defer x.deinit(ally);

    //check solve
    var b_ = try a.mulV(x, ally);
    defer b_.deinit(ally);

    try testing.expect(b.at(0).cmp(.eq, b_.at(0)));
    try testing.expect(b.at(1).cmp(.eq, b_.at(1)));
    try testing.expect(b.at(2).cmp(.eq, b_.at(2)));

    //check det
    try testing.expect(F.from(78, 1).cmp(.eq, lu.det()));

    //reconstuct a from LU
    var l = try lu.lt.transpose(ally);
    defer l.deinit();
    for (0..3) |i| {
        try l.set(i, i, F.eye);
    }

    const a_ = try l.mul(lu.u, ally);
    defer a_.deinit();

    try testing.expect(a.at(0, 0).cmp(.eq, a_.at(0, 0)));
    try testing.expect(a.at(0, 1).cmp(.eq, a_.at(0, 1)));
    try testing.expect(a.at(0, 2).cmp(.eq, a_.at(0, 2)));
    try testing.expect(a.at(1, 0).cmp(.eq, a_.at(1, 0)));
    try testing.expect(a.at(1, 1).cmp(.eq, a_.at(1, 1)));
    try testing.expect(a.at(1, 2).cmp(.eq, a_.at(1, 2)));
    try testing.expect(a.at(2, 0).cmp(.eq, a_.at(2, 0)));
    try testing.expect(a.at(2, 1).cmp(.eq, a_.at(2, 1)));
    try testing.expect(a.at(2, 2).cmp(.eq, a_.at(2, 2)));
}
