const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;
const Allocator = std.mem.Allocator;
const PermutationQueue = @import("../utils/permutationQueue.zig").PermutationQueueType;

/// struct for matrix operations
/// Element must be a field
/// Index must be an unsigned integer <= usize
pub fn MatrixType(comptime Element: type, comptime Index: type) type {
    assert(@typeInfo(Index) == .Int);
    assert(@typeInfo(Index).Int.signedness == .unsigned);
    assert(@typeInfo(Index).Int.bits <= 64);
    const Vector = @import("../vector.zig").VectorType(Element);
    return struct {
        const Matrix = @This();

        const Entry = struct { col: Index, val: Element };
        const EntryList = std.MultiArrayList(Entry);
        const EntrySlice = EntryList.Slice;

        rows: Index,
        cols: Index,
        rptr: [*]usize,
        val: *EntrySlice,
        allocator: Allocator,

        /// deinitialize matrix
        pub fn deinit(a: Matrix) void {
            a.val.deinit(a.allocator);
            a.allocator.destroy(a.val);
            a.allocator.free(a.rptr[0 .. a.rows + 1]);
        }

        ///allocate empty matrix
        pub fn init(rows: Index, cols: Index, allocator: Allocator) !Matrix {
            const r = rows + 1;
            const val = try allocator.create(EntrySlice);
            errdefer allocator.destroy(val);
            const rptr = try allocator.alloc(usize, r);
            errdefer allocator.free(rptr);
            const res = Matrix{
                .val = val,
                .rptr = rptr.ptr,
                .rows = rows,
                .cols = cols,
                .allocator = allocator,
            };
            res.val.* = (EntryList{}).slice();
            @memset(res.rptr[0..r], @as(Index, 0));
            return res;
        }

        ///allocate empty matrix with same dimensions as a
        pub fn like(a: Matrix, allocator: Allocator) !Matrix {
            return init(a.rows, a.cols, allocator);
        }

        /// ensures matrix has given capacity
        /// invalidates values !!
        fn ensureTotalCapacity(a: Matrix, new_capacity: usize) !void {
            if (a.val.capacity < new_capacity) {
                a.val.deinit(a.allocator);
                var new: EntryList = .{};
                try new.setCapacity(a.allocator, new_capacity);
                a.val.* = new.slice();
            }
        }

        /// res <- a
        pub fn copy(a: Matrix, res: Matrix) !void {
            assert(a.rows == res.rows);
            assert(a.cols == res.cols);
            if (a.val != res.val) {
                try res.ensureTotalCapacity(a.val.len);
                res.val.len = a.val.len;
                @memcpy(res.val.items(.col), a.val.items(.col));
                @memcpy(res.val.items(.val), a.val.items(.val));
                const r = res.rows + 1;
                @memcpy(res.rptr[1..r], a.rptr[1..r]);
            }
        }

        const Builder = struct {
            res: Matrix,
            row: Index = 0,
            col: Index = 0,

            /// sets value of an element
            /// only to be used in order
            pub fn set(a: *Builder, row: Index, col: Index, b: Element) void {
                assert(row < a.res.rows and col < a.res.cols); //matrix dimensions
                assert(row > a.row or (row == a.row and col >= a.col)); //input order
                if (b.cmp(.neq, Element.zero)) {
                    @memset(a.res.rptr[a.row + 1 .. row + 1], a.res.val.len);

                    a.res.appendAssumeCapacity(.{ .col = col, .val = b });
                    a.row = row;
                    a.col = col + 1;
                }
            }

            /// finalize the matrix being built
            pub fn fin(a: *Builder) void {
                @memset(a.res.rptr[a.row + 1 .. a.res.rows + 1], a.res.val.len);
            }
        };

        /// clears matrix, ensures capacity and returns a Builder
        pub fn build(res: Matrix, nonzeros: usize) !Builder {
            try res.ensureTotalCapacity(nonzeros);
            res.val.len = 0;
            return .{ .res = res };
        }

        const I = struct { ind: usize, ex: bool };

        /// return index of row, col in matrix
        /// performs binary search
        /// O(log(m)), O(1) if dense
        fn indAt(a: Matrix, row: Index, col: Index) I {
            const r = a.val.items(.col)[a.rptr[row]..a.rptr[row + 1]];
            if (r.len == 0) {
                return .{ .ind = a.rptr[row], .ex = false };
            } else {
                var min: Index = @truncate(r.len -| (a.cols - col));
                var max: Index = @truncate(@min(r.len - 1, col));
                while (min < max) {
                    const pivot = @divFloor(min + max, 2);
                    if (col <= r[pivot]) {
                        max = pivot;
                    } else {
                        min = pivot + 1;
                    }
                }
                return .{ .ind = a.rptr[row] + min, .ex = col == r[min] };
            }
        }

        /// return element at row and column
        /// O(log(m))
        pub fn at(a: Matrix, row: Index, col: Index) Element {
            assert(row < a.rows and col < a.cols); //matrix dimensions
            const i = a.indAt(row, col);
            if (i.ex) {
                return a.val.items(.val)[i.ind];
            } else {
                return Element.zero;
            }
        }

        /// res <- A^T
        /// O(m*n)
        pub fn transpose(a: Matrix, res: Matrix) !void {
            assert(a.cols == res.rows);
            assert(a.rows == res.cols);
            assert(a.val != res.val); //no inplace
            try res.ensureTotalCapacity(a.val.len);
            res.val.len = a.val.len;
            // prepare rowpointers
            @memset(res.rptr[1 .. res.rows + 1], 0);
            for (a.val.items(.col)) |c| {
                if (c + 2 <= res.rows) {
                    res.rptr[c + 2] += 1;
                }
            }
            for (2..res.rows) |i| {
                res.rptr[i + 1] += res.rptr[i];
            }
            // write values + fix rowpointers
            for (0..a.rows) |r| {
                for (a.rptr[r]..a.rptr[r + 1]) |i| {
                    const c = a.val.items(.col)[i];
                    const v = a.val.items(.val)[i];
                    res.val.set(res.rptr[c + 1], .{ .col = @truncate(r), .val = v });
                    res.rptr[c + 1] += 1;
                }
            }
        }

        /// slow way to set element
        /// if possible use Builder instead !!
        pub fn edit(a: Matrix, row: Index, col: Index, b: Element) !void {
            const i = a.indAt(row, col);
            if (b.cmp(.eq, Element.zero)) {
                if (i.ex) {
                    a.val.orderedRemove(i.ind);
                    for (row + 1..a.rows + 1) |j| {
                        a.rptr[j] -= 1;
                    }
                }
            } else {
                if (i.ex) {
                    a.val.items(.col)[i.ind] = b;
                } else {
                    try a.val.insert(a.allocator, i.ind, .{ .col = col, .val = b });
                    for (row + 1..a.rows + 1) |j| {
                        a.rptr[j] += 1;
                    }
                }
            }
        }

        /// res <- b * a
        pub fn mulE(a: Matrix, b: Element, res: Matrix) !void {
            assert(res.rows == a.rows);
            assert(res.cols == a.cols);
            if (b.cmp(.eq, Element.zero)) {
                res.val.len = 0;
                @memset(res.rptr[1 .. res.rows + 1], @as(Index, 0));
            } else {
                if (res.val != a.val) {
                    try res.ensureTotalCapacity(a.val.len);
                    res.val.len = a.val.len;
                    @memcpy(res.rptr[1 .. res.rows + 1], a.rptr[1 .. a.rows + 1]);
                }
                for (0..res.val.len) |i| {
                    res.val.items(.val)[i] = b.mul(a.val.items(.val)[i]);
                }
            }
        }

        /// res <- a / b
        pub fn divE(a: Matrix, b: Element, res: Matrix) !void {
            assert(res.rows == a.rows);
            assert(res.cols == a.cols);
            assert(b.cmp(.neq, Element.zero));
            if (res.val != a.val) {
                try res.ensureTotalCapacity(a.val.len);
                res.val.len = a.val.len;
                @memcpy(res.rptr[1 .. res.rows + 1], a.rptr[1 .. a.rows + 1]);
            }
            for (0..res.val.len) |i| {
                res.val.items(.val)[i] = a.val.items(.val)[i].div(b);
            }
        }

        /// res <- a * b
        pub fn mulV(a: Matrix, b: Vector, res: Vector) void {
            assert(a.cols == b.len);
            assert(a.rows == res.len);
            assert(b.val.ptr != res.val.ptr);
            for (0..a.rows) |r| {
                var sum = Element.zero;
                for (a.rptr[r]..a.rptr[r + 1]) |i| {
                    sum = sum.add(a.val.items(.val)[i].mul(b.at(a.val.items(.col)[i])));
                }
                res.set(r, sum);
            }
        }

        inline fn appendAssumeCapacity(a: Matrix, b: Entry) void {
            assert(a.val.capacity > a.val.len);
            a.val.len += 1;
            a.val.set(a.val.len - 1, b);
        }

        /// count nonzeros of a * b
        fn count_mul(a: Matrix, b: Matrix, buf: []usize) !usize {
            var n: usize = 0;
            for (0..a.rows) |i| {
                for (a.rptr[i]..a.rptr[i + 1]) |j| {
                    const r = a.val.items(.col)[j];
                    buf[r] = b.rptr[r];
                }
                while (true) {
                    var min = b.cols;
                    for (a.rptr[i]..a.rptr[i + 1]) |j| {
                        const r = a.val.items(.col)[j];
                        if (buf[r] < b.rptr[r + 1]) {
                            const c = b.val.items(.col)[buf[r]];
                            if (c < min) {
                                min = c;
                            }
                        }
                    }
                    if (min == b.cols) {
                        break;
                    }
                    for (a.rptr[i]..a.rptr[i + 1]) |j| {
                        const r = a.val.items(.col)[j];
                        if (buf[r] < b.rptr[r + 1]) {
                            const c = b.val.items(.col)[buf[r]];
                            if (c == min) {
                                buf[r] += 1;
                            }
                        }
                    }
                    n += 1;
                }
            }
            return n;
        }

        /// res <- a * b
        pub fn mul(a: Matrix, b: Matrix, res: Matrix) !void {
            assert(a.cols == b.rows);
            assert(res.rows == a.rows);
            assert(res.cols == b.cols);
            assert(a.val != res.val);
            assert(b.val != res.val);
            const buf = try res.allocator.alloc(usize, b.rows);
            defer res.allocator.free(buf);
            try res.ensureTotalCapacity(try a.count_mul(b, buf));
            res.val.len = 0;
            for (0..a.rows) |i| {
                for (a.rptr[i]..a.rptr[i + 1]) |j| {
                    const r = a.val.items(.col)[j];
                    buf[r] = b.rptr[r];
                }
                while (true) {
                    var min = b.cols;
                    for (a.rptr[i]..a.rptr[i + 1]) |j| {
                        const r = a.val.items(.col)[j];
                        if (buf[r] < b.rptr[r + 1]) {
                            const c = b.val.items(.col)[buf[r]];
                            if (c < min) {
                                min = c;
                            }
                        }
                    }
                    if (min == b.cols) {
                        res.rptr[i + 1] = res.val.len;
                        break;
                    }
                    var sum = Element.zero;
                    for (a.rptr[i]..a.rptr[i + 1]) |j| {
                        const r = a.val.items(.col)[j];
                        if (buf[r] < b.rptr[r + 1]) {
                            const c = b.val.items(.col)[buf[r]];
                            const v = b.val.items(.val)[buf[r]];
                            if (c == min) {
                                buf[r] += 1;
                                sum = sum.add(a.val.items(.val)[j].mul(v));
                            }
                        }
                    }
                    if (sum.cmp(.neq, Element.zero)) {
                        res.appendAssumeCapacity(.{ .col = min, .val = sum });
                    }
                }
            }
        }

        /// count nonzeros of a + b
        fn count_add(a: Matrix, b: Matrix) usize {
            var n: usize = 0;
            for (0..a.rows) |r| {
                var i_a = a.rptr[r];
                var i_b = b.rptr[r];
                while (true) {
                    const c_a = if (i_a < a.rptr[r + 1]) a.val.items(.col)[i_a] else a.cols;
                    const c_b = if (i_b < b.rptr[r + 1]) b.val.items(.col)[i_b] else b.cols;
                    if (c_a == c_b) {
                        if (c_a == a.cols) {
                            break;
                        }
                        i_a += 1;
                        i_b += 1;
                    } else if (c_a < c_b) {
                        i_a += 1;
                    } else { //col_b < col_a
                        i_b += 1;
                    }
                    n += 1;
                }
            }
            return n;
        }

        /// res <- a + b
        pub fn add(a: Matrix, b: Matrix, res: Matrix) !void {
            assert(a.rows == b.rows and a.rows == res.rows);
            assert(a.cols == b.cols and a.cols == res.cols);

            try res.ensureTotalCapacity(a.count_add(b));
            res.val.len = 0;

            for (0..a.rows) |r| {
                var i_a = a.rptr[r];
                var i_b = b.rptr[r];
                while (true) {
                    const c_a = if (i_a < a.rptr[r + 1]) a.val.items(.col)[i_a] else a.cols;
                    const c_b = if (i_b < b.rptr[r + 1]) b.val.items(.col)[i_b] else a.cols;
                    if (c_a == c_b) {
                        if (c_a == a.cols) {
                            res.rptr[r + 1] = res.val.len;
                            break;
                        }
                        const sum = a.val.items(.val)[i_a].add(b.val.items(.val)[i_b]);
                        if (sum.cmp(.neq, Element.zero)) {
                            res.appendAssumeCapacity(.{ .col = c_a, .val = sum });
                        }
                        i_a += 1;
                        i_b += 1;
                    } else if (c_a < c_b) {
                        res.appendAssumeCapacity(.{ .col = c_a, .val = a.val.items(.val)[i_a] });
                        i_a += 1;
                    } else { //col_b < col_a
                        res.appendAssumeCapacity(.{ .col = c_b, .val = b.val.items(.val)[i_b] });
                        i_b += 1;
                    }
                }
            }
        }

        /// res <- a - b
        pub fn sub(a: Matrix, b: Matrix, res: Matrix) !void {
            assert(a.rows == b.rows and a.rows == res.rows);
            assert(a.cols == b.cols and a.cols == res.cols);

            try res.ensureTotalCapacity(a.count_add(b));
            res.val.len = 0;

            for (0..a.rows) |r| {
                var i_a = a.rptr[r];
                var i_b = b.rptr[r];
                while (true) {
                    const c_a = if (i_a < a.rptr[r + 1]) a.val.items(.col)[i_a] else a.cols;
                    const c_b = if (i_b < b.rptr[r + 1]) b.val.items(.col)[i_b] else a.cols;
                    if (c_a == c_b) {
                        if (c_a == a.cols) {
                            res.rptr[r + 1] = res.val.len;
                            break;
                        }
                        const sum = a.val.items(.val)[i_a].sub(b.val.items(.val)[i_b]);
                        if (sum.cmp(.neq, Element.zero)) {
                            res.appendAssumeCapacity(.{ .col = c_a, .val = sum });
                        }
                        i_a += 1;
                        i_b += 1;
                    } else if (c_a < c_b) {
                        res.appendAssumeCapacity(.{ .col = c_a, .val = a.val.items(.val)[i_a] });
                        i_a += 1;
                    } else { //col_b < col_a
                        res.appendAssumeCapacity(.{ .col = c_b, .val = b.val.items(.val)[i_b].neg() });
                        i_b += 1;
                    }
                }
            }
        }

        /// frobenius norm = sqrt of sum of all squared entries
        pub fn frobenius(a: Matrix) Element {
            var sum = Element.zero;
            for (a.val.items(.val)) |a_| {
                sum = sum.add(a_.mul(a_));
            }
            return sum.sqrt();
        }

        pub const LU = struct {
            p: [*]Index,
            q: [*]Index,
            lT: Matrix,
            u: Matrix,
            buf: Vector,

            pub fn init(n: Index, allocator: Allocator) !LU {
                var res: LU = undefined;
                res.u = try Matrix.init(n, n, allocator);
                errdefer res.u.deinit();
                res.lT = try Matrix.init(n, n, allocator);
                errdefer res.lT.deinit();
                res.buf = try Vector.init(n, allocator);
                errdefer res.buf.deinit(allocator);
                res.p = (try allocator.alloc(Index, n)).ptr;
                errdefer allocator.free(res.p[0..n]);
                res.q = (try allocator.alloc(Index, n)).ptr;
                errdefer allocator.free(res.q[0..n]);
                return res;
            }

            pub fn deinit(a: LU) void {
                a.u.allocator.free(a.q[0..a.u.rows]);
                a.u.allocator.free(a.p[0..a.u.rows]);
                a.buf.deinit(a.u.allocator);
                a.lT.deinit();
                a.u.deinit();
            }

            pub fn from(a: LU, b: Matrix) !void {
                const s = try Solver.init(b.rows, a.u.allocator);
                defer s.deinit();
                try s.run(b, a);
            }

            pub fn solve(a: LU, b: Vector, res: Vector) void {
                const n = a.buf.len;
                assert(n == b.len and n == res.len);
                b.copy(a.buf);
                //apply lT-1
                for (0..n) |i| {
                    const row = a.p[i];
                    for (a.lT.rptr[row]..a.lT.rptr[row + 1]) |j| {
                        const lT = a.lT.val.get(j);
                        a.buf.set(lT.col, a.buf.at(lT.col).sub(lT.val.mul(a.buf.at(row))));
                    }
                }
                //apply u-1
                for (0..n) |i| {
                    const row = a.p[n - 1 - i];
                    const col_diag = a.q[n - 1 - i];
                    var elem_diag: Element = undefined;
                    var sum = a.buf.at(row);
                    for (a.u.rptr[row]..a.u.rptr[row + 1]) |j| {
                        const u = a.u.val.get(j);
                        if (u.col == col_diag) {
                            elem_diag = u.val;
                        } else {
                            sum = sum.sub(u.val.mul(res.at(u.col)));
                        }
                    }
                    res.set(col_diag, sum.div(elem_diag));
                }
            }

            pub fn det(a: LU) Element {
                var res = Element.eye;
                for (0..a.u.rows) |i| {
                    res = res.mul(a.u.at(a.p[i], a.q[i]));
                }
                return res;
            }

            pub const Solver = struct {
                const NonZeroList = std.ArrayListUnmanaged(Index);
                const PQ = PermutationQueue(Priority, Priority.before, Index);
                u: []EntrySlice,
                lT: []EntrySlice,
                queue: PQ,
                buf: *EntrySlice,
                allocator: Allocator,

                const Priority = struct {
                    nzrows: NonZeroList = .{},
                    n: Index = undefined,
                    norm1: Element = undefined,

                    fn before(self: Priority, other: Priority) bool {
                        if (self.n < other.n) return true;
                        if (self.n == other.n and self.norm1.cmp(.gt, other.norm1)) return true;
                        return false;
                    }
                };

                pub fn init(n: Index, allocator: Allocator) !Solver {
                    var res: Solver = undefined;
                    res.u = try allocator.alloc(EntrySlice, n);
                    errdefer allocator.free(res.u);
                    @memset(res.u, (EntryList{}).slice());
                    res.lT = try allocator.alloc(EntrySlice, n);
                    errdefer allocator.free(res.lT);
                    @memset(res.lT, (EntryList{}).slice());
                    res.queue = try PQ.init(n, allocator);
                    errdefer res.queue.deinit(allocator);
                    @memset(res.queue.priorities[0..n], Priority{});
                    res.buf = try allocator.create(EntrySlice);
                    errdefer allocator.destroy(res.buf);
                    res.buf.* = (EntryList{}).slice();
                    res.allocator = allocator;
                    return res;
                }

                pub fn deinit(a: Solver) void {
                    for (0..a.u.len) |i| {
                        a.u[i].deinit(a.allocator);
                    }
                    a.allocator.free(a.u);
                    for (0..a.lT.len) |i| {
                        a.lT[i].deinit(a.allocator);
                    }
                    a.allocator.free(a.lT);
                    for (0..a.queue.items.len) |i| {
                        a.queue.priorities[i].nzrows.deinit(a.allocator);
                    }
                    a.queue.deinit(a.allocator);
                    a.buf.deinit(a.allocator);
                    a.allocator.destroy(a.buf);
                }

                /// ensures Buffer has given capacity
                /// invalidates values !!
                fn ensureEntrySliceCapacity(a: *EntrySlice, new_capacity: usize, allocator: Allocator) !void {
                    if (a.capacity < new_capacity) {
                        a.deinit(allocator);
                        var new: EntryList = .{};
                        try new.setCapacity(allocator, new_capacity);
                        a.* = new.slice();
                    }
                }

                fn set(a: Solver, b: Matrix) !void {
                    assert(a.u.len == b.rows);
                    assert(a.queue.items.len == b.cols);
                    for (0..b.rows) |i| {
                        try Solver.ensureEntrySliceCapacity(&a.u[i], b.rptr[i + 1] - b.rptr[i], a.allocator);
                        a.u[i].len = b.rptr[i + 1] - b.rptr[i];
                        a.lT[i].len = 0;
                        const cols = b.val.items(.col);
                        @memcpy(a.u[i].items(.col), cols[b.rptr[i]..b.rptr[i + 1]]);
                        @memcpy(a.u[i].items(.val), b.val.items(.val)[b.rptr[i]..b.rptr[i + 1]]);
                        for (b.rptr[i]..b.rptr[i + 1]) |j| {
                            a.queue.priorities[cols[j]].nzrows.items.len += 1;
                        }
                    }
                    for (0..b.cols) |i| {
                        const n = a.queue.priorities[i].nzrows.items.len;
                        a.queue.priorities[i].nzrows.items.len = 0;
                        try a.queue.priorities[i].nzrows.ensureTotalCapacity(a.allocator, n);
                    }
                    var i: Index = 0;
                    while (i < b.rows) : (i += 1) {
                        const cols = b.val.items(.col);
                        for (b.rptr[i]..b.rptr[i + 1]) |j| {
                            a.queue.priorities[cols[j]].nzrows.appendAssumeCapacity(i);
                        }
                    }
                    i = 0;

                    while (i < b.cols) : (i += 1) {
                        a.updateQueue(i);
                    }
                }

                /// return index of row, col in matrix
                /// performs binary search
                /// O(log(m)), O(1) if dense
                /// TODO: refactor with Matrix.indAt
                fn indAt(a: Solver, row: Index, col: Index) Index {
                    const r = a.u[row].items(.col);
                    if (r.len == 0) {
                        return 0;
                    } else {
                        var min: Index = @truncate(r.len -| (a.queue.items.len - col));
                        var max: Index = @truncate(@min(r.len - 1, col));
                        while (min < max) {
                            const pivot = @divFloor(min + max, 2);
                            if (col <= r[pivot]) {
                                max = pivot;
                            } else {
                                min = pivot + 1;
                            }
                        }
                        return min;
                    }
                }

                fn updateQueue(a: Solver, col: Index) void {
                    var norm1 = Element.zero;
                    const prio = &a.queue.priorities[col];
                    var i: usize = 0;
                    while (i < prio.nzrows.items.len) {
                        const row = prio.nzrows.items[i];
                        const ind = a.indAt(row, col);
                        if (col == a.u[row].items(.col)[ind]) {
                            norm1 = norm1.add(a.u[row].items(.val)[ind].abs());
                            i += 1;
                        } else {
                            _ = prio.nzrows.swapRemove(i);
                        }
                    }
                    prio.n = @intCast(prio.nzrows.items.len);
                    prio.norm1 = norm1;
                    a.queue.set(col, prio.*);
                }

                fn countElim(a: Solver, src_row: Index, trg_row: Index) usize {
                    var n: Index = 0;
                    var i_src: Index = 0;
                    var i_trg: Index = 0;
                    while (true) {
                        const c_src = if (i_src < a.u[src_row].len) a.u[src_row].items(.col)[i_src] else a.queue.items.len;
                        const c_trg = if (i_trg < a.u[trg_row].len) a.u[trg_row].items(.col)[i_trg] else a.queue.items.len;
                        if (c_src == c_trg) {
                            if (c_src == a.queue.items.len) {
                                break;
                            }
                            i_src += 1;
                            i_trg += 1;
                        } else if (c_src < c_trg) {
                            i_src += 1;
                        } else { //col_b < col_a
                            i_trg += 1;
                        }
                        n += 1;
                    }
                    return n - 1;
                }

                fn elim(a: Solver, ind_src: Index, row_src: Index, row_trg: Index) !Element {
                    try ensureEntrySliceCapacity(a.buf, a.countElim(row_src, row_trg), a.allocator);
                    a.buf.len = 0;

                    const vals_src = a.u[row_src].items(.val);
                    const vals_trg = a.u[row_trg].items(.val);
                    const cols_src = a.u[row_src].items(.col);
                    const cols_trg = a.u[row_trg].items(.col);
                    const col = cols_src[ind_src];

                    const factor = blk: {
                        const val_src = vals_src[ind_src];
                        const val_trg = vals_trg[a.indAt(row_trg, col)];
                        break :blk val_trg.div(val_src);
                    };

                    var i_src: Index = 0;
                    var i_trg: Index = 0;
                    while (true) {
                        const col_src = if (i_src < a.u[row_src].len) cols_src[i_src] else a.queue.items.len;
                        const col_trg = if (i_trg < a.u[row_trg].len) cols_trg[i_trg] else a.queue.items.len;
                        if (col_src == col_trg) {
                            if (col_src == a.queue.items.len) break; //end
                            if (col_src != col) { //skip at pivot column
                                const val = vals_trg[i_trg].sub(factor.mul(vals_src[i_src]));
                                if (val.cmp(.neq, Element.zero)) {
                                    a.buf.len += 1;
                                    a.buf.set(a.buf.len - 1, .{ .col = col_src, .val = val });
                                }
                            }
                            i_src += 1;
                            i_trg += 1;
                        } else if (col_src < col_trg) {
                            a.buf.len += 1;
                            a.buf.set(a.buf.len - 1, .{ .col = col_src, .val = factor.mul(vals_src[i_src]).neg() });
                            try a.queue.priorities[col_src].nzrows.append(a.allocator, row_trg);
                            i_src += 1;
                        } else { //col_b < col_a
                            a.buf.len += 1;
                            a.buf.set(a.buf.len - 1, .{ .col = col_src, .val = vals_trg[i_trg] });
                            i_trg += 1;
                        }
                    }

                    // swap buf and trg
                    const helper = a.buf.*;
                    a.buf.* = a.u[row_trg];
                    a.u[row_trg] = helper;

                    return factor;
                }

                fn choosePvtRow(a: Solver, col: Index) struct { row: Index, ind: Index } {
                    const nzrows = a.queue.priorities[col].nzrows;
                    var row: Index = nzrows.items[0];
                    var ind: Index = a.indAt(row, col);
                    var min: usize = a.u[0].len;
                    var max = a.u[0].items(.val)[ind].abs();
                    for (1..nzrows.items.len) |j| {
                        const r = nzrows.items[j];
                        if (a.u[r].len > min) continue;
                        const i = a.indAt(r, col);
                        const abs = a.u[r].items(.val)[i].abs();
                        if (a.u[r].len < min or abs.cmp(.gt, max)) {
                            min = a.u[r].len;
                            max = abs;
                            ind = i;
                            row = r;
                        }
                    }
                    return .{ .row = row, .ind = ind };
                }

                fn removeRowFromComb(a: Solver, row: Index) void {
                    for (a.u[row].items(.col)) |c| {
                        for (a.queue.priorities[c].nzrows.items, 0..) |r, i| {
                            if (r == row) {
                                _ = a.queue.priorities[c].nzrows.swapRemove(i);
                                break;
                            }
                        }
                    }
                }

                fn elimCol(a: Solver, col: Index) !Index {
                    const pvt = a.choosePvtRow(col);
                    a.removeRowFromComb(pvt.row);
                    const n = a.queue.priorities[col].nzrows.items.len;
                    try ensureEntrySliceCapacity(&a.lT[pvt.row], n, a.allocator);
                    for (0..n) |_| {
                        var r: Index = 0;
                        const nz = a.queue.priorities[col].nzrows.items;
                        var i: Index = 1;
                        while (i < nz.len) : (i += 1) {
                            if (nz[i] < nz[r]) {
                                r = i;
                            }
                        }
                        const trg_row = nz[r];
                        _ = a.queue.priorities[col].nzrows.swapRemove(r);

                        const factor = try a.elim(pvt.ind, pvt.row, trg_row);
                        a.lT[pvt.row].len += 1;
                        a.lT[pvt.row].set(a.lT[pvt.row].len - 1, .{ .val = factor, .col = trg_row });
                    }
                    for (0..a.u[pvt.row].len) |i| {
                        if (i != pvt.ind) {
                            a.updateQueue(a.u[pvt.row].items(.col)[i]);
                        }
                    }
                    return pvt.row;
                }

                fn build(rows: []EntrySlice, res: Matrix) !void {
                    assert(rows.len == res.rows);
                    var count: usize = 0;
                    for (rows) |row| {
                        count += row.len;
                    }
                    try res.ensureTotalCapacity(count);
                    res.val.len = count;
                    for (rows, 0..) |row, i| {
                        res.rptr[i + 1] = res.rptr[i] + row.len;
                        @memcpy(res.val.items(.val)[res.rptr[i]..res.rptr[i + 1]].ptr, row.items(.val));
                        @memcpy(res.val.items(.col)[res.rptr[i]..res.rptr[i + 1]].ptr, row.items(.col));
                    }
                }

                pub fn run(a: Solver, b: Matrix, lu: LU) !void {
                    try a.set(b);
                    assert(lu.u.rows == a.u.len);
                    for (0..a.queue.items.len) |i| {
                        const col = a.queue.remove();
                        lu.q[i] = col;
                        if (a.queue.priorities[col].n == 0) return error.NotInvertible;
                        const row = try a.elimCol(col);
                        lu.p[i] = row;
                    }
                    try Solver.build(a.u, lu.u);
                    try Solver.build(a.lT, lu.lT);
                }
            };
        };
    };
}

test "matrix builder" {
    const n = 5;
    const m = 8;

    const ally = testing.allocator;
    const F = @import("../field.zig").Float(f32);
    const M = MatrixType(F, u8);

    const a = try M.init(n, m, ally); //const only refers to the dimensions
    defer a.deinit();
    try testing.expect(a.val.len == 0);
    var b = try a.build(n * m);
    try testing.expect(a.val.capacity == n * m);
    for (0..n) |i| {
        for (0..m) |j| {
            b.set(@intCast(i), @intCast(j), F.from(@intCast(i * m + j), 1));
        }
    }
    b.fin();
    try testing.expect(a.val.len == n * m - 1);
    for (0..n) |i| {
        for (0..m) |j| {
            try testing.expectEqual(F.from(@intCast(i * m + j), 1), a.at(@intCast(i), @intCast(j)));
        }
    }

    b = try a.build(n);
    try testing.expect(a.val.capacity == n * m);
    for (0..@min(n, m)) |i| {
        b.set(@intCast(i), @intCast(i), F.from(@intCast(i), 1));
    }
    b.fin();
    try testing.expect(a.val.len == n - 1);
    for (0..n) |i| {
        for (0..m) |j| {
            if (i == j) {
                try testing.expectEqual(F.from(@intCast(i), 1), a.at(@intCast(i), @intCast(j)));
            } else {
                try testing.expectEqual(F.zero, a.at(@intCast(i), @intCast(j)));
            }
        }
    }
}

test "matrix transpose" {
    const ally = std.testing.allocator;
    const F = @import("../field.zig").Float(f32);
    const M = MatrixType(F, u8);

    const a = try M.init(3, 4, ally);
    defer a.deinit();
    var a_ = try a.build(9);
    var c: isize = 1;
    for (0..3) |i| {
        for (0..4) |j| {
            if (j >= i) {
                a_.set(@intCast(i), @intCast(j), F.from(c, 1));
                c += 1;
            }
        }
    }
    a_.fin();

    const b = try M.init(4, 3, ally);
    defer b.deinit();
    try a.transpose(b);
    c = 1;
    for (0..3) |i| {
        for (0..4) |j| {
            if (j >= i) {
                try testing.expectEqual(F.from(c, 1), b.at(@intCast(j), @intCast(i)));
                c += 1;
            }
        }
    }
}

test "matrix addition" {
    const ally = std.testing.allocator;
    const F = @import("../field.zig").Float(f32);
    const M = MatrixType(F, u8);

    // 1 2 3   1 0 0   0 2 3
    // 0 0 4 - 2 0 0 =-2 0 4
    // 0 0 0   3 4 0  -3-4 0

    const n = 3;
    const a = try M.init(n, n, ally);
    defer a.deinit();
    var a_ = try a.build(4);
    a_.set(0, 0, F.from(1, 1));
    a_.set(0, 1, F.from(2, 1));
    a_.set(0, 2, F.from(3, 1));
    a_.set(1, 2, F.from(4, 1));
    a_.fin();

    const b = try M.init(n, n, ally);
    defer b.deinit();
    try a.transpose(b);
    try b.mulE(F.from(-1, 1), b);

    var c = try M.init(n, n, ally);
    defer c.deinit();
    try a.add(b, c);

    try testing.expect(c.val.capacity == 7);
    try testing.expect(c.val.len == 6);
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
    const M = MatrixType(F, u8);

    // 1 2 3   1 0-1   1 5 0
    // 0 4 5 * 0 1-1 = 0 9 1
    // 0 0 6   0 1 1   0 6 6

    const n = 3;
    const a = try M.init(n, n, ally);
    defer a.deinit();
    var a_ = try a.build(6);
    a_.set(0, 0, F.from(1, 1));
    a_.set(0, 1, F.from(2, 1));
    a_.set(0, 2, F.from(3, 1));
    a_.set(1, 1, F.from(4, 1));
    a_.set(1, 2, F.from(5, 1));
    a_.set(2, 2, F.from(6, 1));
    a_.fin();

    const b = try M.init(n, n, ally);
    defer b.deinit();
    var b_ = try b.build(6);
    b_.set(0, 0, F.from(1, 1));
    b_.set(0, 2, F.from(-1, 1));
    b_.set(1, 1, F.from(1, 1));
    b_.set(1, 2, F.from(-1, 1));
    b_.set(2, 1, F.from(1, 1));
    b_.set(2, 2, F.from(1, 1));
    b_.fin();

    var c = try M.init(n, n, ally);
    defer c.deinit();
    try a.mul(b, c);

    try testing.expect(c.val.capacity == 7);
    try testing.expect(c.val.len == 6);
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
    const M = MatrixType(F, u8);

    const n = 3;
    const a = try M.init(n, n, ally);
    defer a.deinit();
    var a_ = try a.build(n);
    for (0..n) |i| {
        a_.set(@intCast(i), @intCast(i), F.eye);
    }
    a_.fin();
    const b = F.from(-314, 100);

    a.mulE(b, a) catch unreachable; //since inplace
    for (0..n) |i| {
        for (0..n) |j| {
            try testing.expect(a.at(@intCast(i), @intCast(j)).cmp(.eq, if (i == j) b else F.zero));
        }
    }

    a.divE(b, a) catch unreachable; //since inplace
    for (0..n) |i| {
        for (0..n) |j| {
            try testing.expect(a.at(@intCast(i), @intCast(j)).cmp(.eq, if (i == j) F.eye else F.zero));
        }
    }
}

test "matrix mulpiplication with vector" {
    const ally = std.testing.allocator;
    const F = @import("../field.zig").Float(f32);
    const V = @import("../vector.zig").VectorType(F);
    const M = MatrixType(F, u8);

    // 1 2 3   -2   -1
    // 0 4 5 * -1 =  1
    // 0 0 6    1    6

    const n = 3;
    const a = try M.init(n, n, ally);
    defer a.deinit();
    var a_ = try a.build(6);
    a_.set(0, 0, F.from(1, 1));
    a_.set(0, 1, F.from(2, 1));
    a_.set(0, 2, F.from(3, 1));
    a_.set(1, 1, F.from(4, 1));
    a_.set(1, 2, F.from(5, 1));
    a_.set(2, 2, F.from(6, 1));
    a_.fin();

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
    const M = MatrixType(F, u8);

    // 1  2  3
    // 2  1  2
    // 0  1  1

    const n = 3;
    const a = try M.init(n, n, ally);
    defer a.deinit();
    var a_ = try a.build(8);
    a_.set(0, 0, F.from(1, 1));
    a_.set(0, 1, F.from(2, 1));
    a_.set(0, 2, F.from(3, 1));
    a_.set(1, 0, F.from(2, 1));
    a_.set(1, 1, F.from(1, 1));
    a_.set(1, 2, F.from(2, 1));
    a_.set(2, 1, F.from(1, 1));
    a_.set(2, 2, F.from(1, 1));
    a_.fin();

    try testing.expectEqual(F.from(5, 1), a.frobenius());
}

test "matrix solve" {
    const ally = std.testing.allocator;
    const F = @import("../field.zig").Float(f32);
    const M = MatrixType(F, u8);
    const V = @import("../vector.zig").VectorType(F);

    // 1  2  3  -0.5   1
    // 2  2  2 * 1.5 = 1
    // 0  1  1  -0.5   1

    const n = 3;
    const a = try M.init(n, n, ally);
    defer a.deinit();
    var a_ = try a.build(8);
    a_.set(0, 0, F.from(1, 1));
    a_.set(0, 1, F.from(2, 1));
    a_.set(0, 2, F.from(3, 1));
    a_.set(1, 0, F.from(2, 1));
    a_.set(1, 1, F.from(2, 1));
    a_.set(1, 2, F.from(2, 1));
    a_.set(2, 1, F.from(1, 1));
    a_.set(2, 2, F.from(1, 1));
    a_.fin();

    const b = try V.init(n, ally);
    defer b.deinit(ally);
    b.fill(F.eye);

    const lu = try M.LU.init(n, ally);
    defer lu.deinit();
    try lu.from(a);

    const x = try b.like(ally);
    defer x.deinit(ally);
    lu.solve(b, x);

    try testing.expectEqual(F.from(-1, 2), x.at(0));
    try testing.expectEqual(F.from(3, 2), x.at(1));
    try testing.expectEqual(F.from(-1, 2), x.at(2));

    try testing.expectEqual(F.from(2, 1), lu.det());
}
