const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;
const Allocator = std.mem.Allocator;
const PPQ = @import("utils/permutationPriorityQueue.zig").PermutationPriorityQueue;

/// struct for matrix operations
/// Element must be a field
/// Index must be an unsigned integer <= usize
pub fn MatrixType(comptime Element: type, comptime Index: type) type {
    //TODO: assertions for Index type
    return struct {
        const Matrix = @This();
        const Vector = @import("vector.zig").VectorType(Element);

        const Entry = struct { col: Index, val: Element };
        const MultiArrayList = std.MultiArrayList(Entry);
        const MultiSlice = MultiArrayList.Slice;
        const max_index = ~@as(Index, 0);

        rows: Index,
        cols: Index,
        rptr: [*]usize,
        val: *MultiSlice,
        //TODO: check usize/Index
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
            const val = try allocator.create(MultiSlice);
            errdefer allocator.destroy(val);
            const rptr = try allocator.alloc(usize, r);
            const res = Matrix{
                .val = val,
                .rptr = rptr.ptr,
                .rows = rows,
                .cols = cols,
                .allocator = allocator,
            };
            res.val.* = (MultiArrayList{}).slice();
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
                var new: MultiArrayList = .{};
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
                    res.val.items(.val)[i] = a.val.items(.val)[i].div(b);
                }
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
            a.val.len += 1;
            a.val.set(a.val.len - 1, b);
        }

        /// count nonzeros of a * b
        fn count_mul(a: Matrix, b: Matrix, buf: []Index) !usize {
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
            const buf = try res.allocator.alloc(Index, b.rows);
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
                    const c_a = if (i_a < a.rptr[r + 1]) a.val.items(.col)[i_a] else max_index;
                    const c_b = if (i_b < b.rptr[r + 1]) b.val.items(.col)[i_b] else max_index;
                    if (c_a == c_b) {
                        if (c_a == max_index) {
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
                    const c_a = if (i_a < a.rptr[r + 1]) a.val.items(.col)[i_a] else max_index;
                    const c_b = if (i_b < b.rptr[r + 1]) b.val.items(.col)[i_b] else max_index;
                    if (c_a == c_b) {
                        if (c_a == max_index) {
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

        const Solver = struct {
            const IndexSet = std.AutoArrayHashMapUnmanaged(Index, void);
            min_ptr: *Index,
            max_ptr: *Index,
            val: []MultiSlice,
            nzs: []IndexSet,
            buf: MultiSlice,
            allocator: Allocator,

            pub fn from(a: Matrix, allocator: Allocator) !Solver {
                const val = try allocator.alloc(MultiSlice, a.rows);
                errdefer a.allocator.free(val);
                const nzs = try allocator.alloc(IndexSet, a.cols);
                errdefer a.allocator.free(nzs);
                var res = .{
                    .min_ptr = &a.val.items(.col)[0],
                    .max_ptr = &a.val.items(.col)[a.val.len - 1],
                    .val = val,
                    .nzs = nzs,
                    .buf = (MultiArrayList{}).slice(),
                    .allocator = allocator,
                };
                const col_ptr: [*]Index = @ptrCast(@alignCast(a.val.ptrs[0]));
                const val_ptr: [*]Index = @ptrCast(@alignCast(a.val.ptrs[1]));
                @memset(res.nzs, IndexSet{});
                for (0..a.rows) |i| {
                    res.val[i].len = a.rptr[i + 1] - a.rptr[i];
                    res.val[i].capacity = res.val[i].len;
                    res.val[i].ptrs[0] = @ptrCast(@alignCast(&col_ptr[a.rptr[i]]));
                    res.val[i].ptrs[1] = @ptrCast(@alignCast(&val_ptr[a.rptr[i]]));
                    for (a.rptr[i]..a.rptr[i + 1]) |j| {
                        res.nzs[a.val.items(.col)[j]].entries.len += 1;
                    }
                }
                for (0..a.cols) |i| {
                    const cap = res.nzs[i].entries.len;
                    res.nzs[i].entries.len = 0;
                    try res.nzs[i].ensureTotalCapacity(a.allocator, cap);
                    errdefer res.nzs[i].deinit(a.allocator);
                }
                for (0..a.rows) |i| {
                    for (a.rptr[i]..a.rptr[i + 1]) |j| {
                        res.nzs[a.val.items(.col)[j]].putAssumeCapacity(i, undefined);
                    }
                }
                return res;
            }

            inline fn readOnly(a: *Solver, slice: MultiSlice) bool {
                if (slice.capacity == 0) return true; //since deinitializing at cap 0 causes segfaults
                const ptr: *Index = &slice.items(.col)[0];
                return @intFromPtr(a.min_ptr) <= @intFromPtr(ptr) and @intFromPtr(ptr) <= @intFromPtr(a.max_ptr);
            }

            fn ensureBufferCapacity(a: *Solver, new_capacity: usize) !void {
                const r = a.readOnly(a.buf);
                if (r or a.buf.capacity < new_capacity) {
                    if (!r) a.buf.deinit(a.allocator);
                    var new: MultiArrayList = .{};
                    try new.setCapacity(a.allocator, new_capacity);
                    a.buf = new.slice();
                }
            }

            pub fn deinit(a: *Solver) void {
                for (0..a.val.len) |i| {
                    if (!a.readOnly(a.val[i])) {
                        a.val[i].deinit(a.allocator);
                    }
                }
                a.allocator.free(a.val);
                for (0..a.nzs.len) |i| {
                    a.nzs[i].deinit(a.allocator);
                }
                a.allocator.free(a.nzs);
                if (!a.readOnly(a.buf)) {
                    a.buf.deinit(a.allocator);
                }
                a.* = undefined;
            }

            /// return index of row, col in matrix
            /// performs binary search
            /// O(log(m)), O(1) if dense
            /// TODO: refactor with Matrix.indAt
            fn indAt(a: Solver, row: Index, col: Index) Index {
                const r = a.val[row].items(.col);
                if (r.len == 0) {
                    return 0;
                } else {
                    var min: Index = @truncate(r.len -| (a.nzs.len - col));
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

            /// count nonzeros of eliminated row
            /// TODO: refactor with Matrix
            fn count_elim(a: Solver, src_row: Index, trg_row: Index) usize {
                var n: Index = 0;
                var i_src: Index = 0;
                var i_trg: Index = 0;
                while (true) {
                    const c_src = if (i_src < a.val[src_row].len) a.val[src_row].items(.col)[i_src] else max_index;
                    const c_trg = if (i_trg < a.val[trg_row].len) a.val[trg_row].items(.col)[i_trg] else max_index;
                    if (c_src == c_trg) {
                        if (c_src == max_index) {
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

            /// eliminate entry in traget row and col using the source row
            /// TODO: refactor with Matrix
            fn elim(a: *Solver, row_src: Index, row_trg: Index, ind_src: Index) !void {
                try a.ensureBufferCapacity(a.count_elim(row_src, row_trg));
                a.buf.len = 0;

                const vals_src = a.val[row_src].items(.val);
                const vals_trg = a.val[row_trg].items(.val);
                const cols_src = a.val[row_src].items(.col);
                const cols_trg = a.val[row_trg].items(.col);

                const col_pvt = cols_src[ind_src];
                const ind_trg = a.indAt(row_trg, col_pvt);
                const factor = vals_trg[ind_trg].div(vals_src[ind_src]).neg();

                // TODO: buf <- trg + factor * src
                var i_src: Index = 0;
                var i_trg: Index = 0;
                while (true) {
                    const col_src = if (i_src < a.val[row_src].len) cols_src[i_src] else max_index;
                    const col_trg = if (i_trg < a.val[row_trg].len) cols_trg[i_trg] else max_index;
                    if (col_src == col_trg) {
                        if (col_src == max_index) break; //end
                        if (col_src != col_pvt) { //skip at pivot column
                            const val = vals_trg[i_trg].add(factor.mul(vals_src[i_src]));
                            if (val.cmp(.neq, Element.zero)) {
                                a.buf.len += 1;
                                a.buf.set(a.buf.len - 1, .{ .col = col_src, .val = val });
                            } else {
                                _ = a.nzs[col_trg].swapRemove(row_trg);
                            }
                        }
                        i_src += 1;
                        i_trg += 1;
                    } else if (col_src < col_trg) {
                        a.buf.len += 1;
                        a.buf.set(a.buf.len - 1, .{ .col = col_src, .val = factor.mul(vals_src[i_src]) });
                        try a.nzs[col_trg].put(a.allocator, row_trg, undefined);
                        i_src += 1;
                    } else { //col_b < col_a
                        a.buf.len += 1;
                        a.buf.set(a.buf.len - 1, .{ .col = col_src, .val = vals_trg[i_trg] });
                        i_trg += 1;
                    }
                }

                // swap buf and trg
                const helper = a.buf;
                a.buf = a.val[row_trg];
                a.val[row_trg] = helper;
            }
        };

        // const LU = struct {
        //     //PAQ=LU
        //     p: [*]Index,
        //     q: [*]Index,
        //     lt: Matrix,
        //     u: Matrix,
        //     buffer: Vector,

        //     const RowPriority = struct {
        //         const Self = @This();
        //         nonzeros: Index,
        //         norm1: Element,

        //         pub fn from(row: Row) Self {
        //             var row_priority = RowPriority{
        //                 .nonzeros = row.len,
        //                 .norm1 = Element.zero,
        //             };
        //             for (0..row.len) |j| {
        //                 row_priority.norm1 = row_priority.norm1.add(row.items(.val)[j].abs());
        //             }
        //             return row_priority;
        //         }

        //         pub fn before(self: Self, other: Self) bool {
        //             if (self.nonzeros < other.nonzeros) return true;
        //             if (self.nonzeros == other.nonzeros and self.norm1.cmp(.gt, other.norm1)) return true;
        //             return false;
        //         }
        //     };

        //     pub fn from(a: Matrix, allocator: Allocator) !LU {
        //         var u = try a.like(allocator);
        //         errdefer u.deinit();
        //         try a.copy(u);
        //         return try LU.decompose(u, allocator);
        //     }

        //     /// O(n*m*(m+log(n)))
        //     pub fn decompose(u: Matrix, allocator: Allocator) !LU {
        //         assert(u.rows == u.cols);
        //         const n = u.rows;

        //         //allocating result structs
        //         var lt = try Matrix.init(n, n, allocator); //l^T for faster fill
        //         errdefer lt.deinit();
        //         var q = try allocator.alloc(Index, n);
        //         errdefer allocator.free(q);

        //         // count nonzeros for each column O(n*m)
        //         var nz_col = try allocator.alloc(IndexSet, n);
        //         defer {
        //             for (0..n) |i| {
        //                 nz_col[i].deinit(allocator);
        //             }
        //             allocator.free(nz_col);
        //         }
        //         for (0..n) |i| {
        //             nz_col[i] = IndexSet{};
        //         }
        //         for (0..n) |i| { //n
        //             for (0..u.lenAt(i)) |j| { //m
        //                 try nz_col[u.colAt(i, j)].put(allocator, i, undefined); //O(1)
        //             }
        //         }

        //         //initialize the priority queue for row pivoting
        //         var row_queue = try PPQ(RowPriority, RowPriority.before, Index).init(n, allocator); //O(n)
        //         errdefer row_queue.deinit(allocator);
        //         for (0..n) |i| { //n
        //             row_queue.set(i, RowPriority.from(u.val[i])); //O(log(n)+m)
        //         }

        //         //main loop
        //         for (0..n) |i| { //n
        //             //get pivot row from queue
        //             const row_pvt = row_queue.remove(); //O(log(n))
        //             if (u.lenAt(row_pvt) == 0) return error.NonInvertible;

        //             //find pivot column O(m)
        //             var ind_pvt: usize = undefined;
        //             {
        //                 var nz_pvt: usize = n;
        //                 var abs_pvt = Element.zero;
        //                 for (0..u.lenAt(row_pvt)) |j| { //m
        //                     const nz_j = nz_col[u.colAt(row_pvt, j)].count();
        //                     const abs_j = u.valAt(row_pvt, j, true).abs();
        //                     if (nz_j < nz_pvt or (nz_j == nz_pvt and abs_j.cmp(.gt, abs_pvt))) {
        //                         nz_pvt = nz_j;
        //                         abs_pvt = abs_j;
        //                         ind_pvt = j;
        //                     }
        //                 }
        //             }
        //             const col_pvt = u.colAt(row_pvt, ind_pvt);
        //             q[i] = col_pvt;

        //             //remove pivot row from nonzeros
        //             for (0..u.lenAt(row_pvt)) |j| {
        //                 _ = nz_col[u.colAt(row_pvt, j)].swapRemove(row_pvt);
        //             }

        //             //elimination
        //             try lt.val[row_pvt].ensureTotalCapacity(allocator, nz_col[col_pvt].count());
        //             var elim_iter = nz_col[col_pvt].iterator();
        //             while (elim_iter.next()) |entry| { //m
        //                 const row_trg = entry.key_ptr.*;
        //                 const factor = try u.elimAt(row_pvt, row_trg, ind_pvt, nz_col); //O(m)
        //                 lt.val[row_pvt].appendAssumeCapacity(.{ .col = row_trg, .val = factor });
        //                 row_queue.set(row_trg, RowPriority.from(u.val[row_trg])); //update priority //O(log(n))
        //             }
        //         }
        //         const p = row_queue.deinitToPermutation(allocator);
        //         allocator.free(p.inv[0..n]);

        //         return LU{ .p = p.val, .q = q.ptr, .lt = lt, .u = u, .buffer = try Vector.init(n, allocator) };
        //     }

        //     pub fn deinit(lu: LU, allocator: Allocator) void {
        //         const n = lu.u.rows;
        //         allocator.free(lu.p[0..n]);
        //         allocator.free(lu.q[0..n]);
        //         lu.lt.deinit();
        //         lu.u.deinit();
        //         lu.buffer.deinit(allocator);
        //     }

        //     pub fn det(lu: LU) Element {
        //         var res = Element.eye;
        //         const n = lu.u.rows;
        //         for (0..n) |i| {
        //             res = res.mul(lu.u.at(lu.p[i], lu.q[i]));
        //         }
        //         return res;
        //     }

        //     pub fn solve(lu: LU, b: Vector, x: Vector) !void {
        //         const n = b.len;
        //         assert(x.len == n);
        //         assert(lu.u.rows == n);

        //         //apply l-1
        //         b.copy(lu.buffer);
        //         for (0..n) |i| {
        //             const row = lu.p[i];
        //             for (0..lu.lt.lenAt(row)) |j| {
        //                 const col = lu.lt.colAt(row, j);
        //                 lu.buffer.set(col, lu.buffer.at(col).sub(lu.lt.valAt(row, j, true).mul(lu.buffer.at(row))));
        //             }
        //         }

        //         //apply u-1
        //         for (0..n) |i| {
        //             const row = lu.p[n - 1 - i];
        //             const col_diag = lu.q[n - 1 - i];
        //             var elem_diag: Element = undefined;
        //             var sum = lu.buffer.at(row);
        //             for (0..lu.u.lenAt(row)) |j| {
        //                 const col = lu.u.colAt(row, j);
        //                 if (col == col_diag) {
        //                     elem_diag = lu.u.valAt(row, j, true);
        //                 } else {
        //                     sum = sum.sub(lu.u.valAt(row, j, true).mul(x.at(col)));
        //                 }
        //             }
        //             x.set(col_diag, sum.div(elem_diag));
        //         }
        //     }
        // };

        // fn debugprint(self: Matrix) void {
        //     std.debug.print("\n", .{});
        //     for (0..self.rows) |i| {
        //         for (0..self.cols) |j| {
        //             std.debug.print("{}, ", .{self.at(i, j).f});
        //         }
        //         std.debug.print("\n", .{});
        //     }
        // }

        // /// O(n*m*(m+log(n)))
        // pub fn decomp(self: Matrix, allocator: Allocator) !LU {
        //     return LU.from(self, allocator);
        // }
    };
}

test "matrix builder" {
    const n = 5;
    const m = 8;

    const ally = testing.allocator;
    const F = @import("field.zig").Float(f32);
    const M = MatrixType(F, usize);

    const a = try M.init(n, m, ally); //const only refers to the dimensions
    defer a.deinit();
    try testing.expect(a.val.len == 0);
    var b = try a.build(n * m);
    try testing.expect(a.val.capacity == n * m);
    for (0..n) |i| {
        for (0..m) |j| {
            b.set(i, j, F.from(@intCast(i * m + j), 1));
        }
    }
    b.fin();
    try testing.expect(a.val.len == n * m - 1);
    for (0..n) |i| {
        for (0..m) |j| {
            try testing.expectEqual(F.from(@intCast(i * m + j), 1), a.at(i, j));
        }
    }

    b = try a.build(n);
    try testing.expect(a.val.capacity == n * m);
    for (0..@min(n, m)) |i| {
        b.set(i, i, F.from(@intCast(i), 1));
    }
    b.fin();
    try testing.expect(a.val.len == n - 1);
    for (0..n) |i| {
        for (0..m) |j| {
            if (i == j) {
                try testing.expectEqual(F.from(@intCast(i), 1), a.at(i, j));
            } else {
                try testing.expectEqual(F.zero, a.at(i, j));
            }
        }
    }
}

test "matrix transpose" {
    const ally = std.testing.allocator;
    const F = @import("field.zig").Float(f32);
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

test "matrix multiplication" {
    const ally = std.testing.allocator;
    const F = @import("field.zig").Float(f32);
    const M = MatrixType(F, usize);

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

test "matrix addition" {
    const ally = std.testing.allocator;
    const F = @import("field.zig").Float(f32);
    const M = MatrixType(F, usize);

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

test "matrix fragmentation" {
    const ally = std.testing.allocator;
    const F = @import("field.zig").Float(f32);
    const M = MatrixType(F, usize);

    // 1 2 3
    // 0 0 4
    // 0 0 0

    const n = 3;
    const a = try M.init(n, n, ally);
    defer a.deinit();
    var a_ = try a.build(4);
    a_.set(0, 0, F.from(1, 1));
    a_.set(0, 1, F.from(2, 1));
    a_.set(0, 2, F.from(3, 1));
    a_.set(1, 2, F.from(4, 1));
    a_.fin();

    var s = try M.Solver.from(a, ally);
    defer s.deinit();

    try testing.expect(s.count_elim(0, 1) == 2);

    try s.elim(1, 0, 0);

    try testing.expect(s.val[0].len == 2);
}

test "matrix multiplication with element" {
    const ally = std.testing.allocator;
    const F = @import("field.zig").Float(f32);
    const M = MatrixType(F, usize);

    const n = 3;
    const a = try M.init(n, n, ally);
    defer a.deinit();
    var a_ = try a.build(n);
    for (0..n) |i| {
        a_.set(i, i, F.eye);
    }
    a_.fin();
    const b = F.from(-314, 100);

    a.mulE(b, a) catch unreachable; //since inplace
    for (0..n) |i| {
        for (0..n) |j| {
            try testing.expect(a.at(i, j).cmp(.eq, if (i == j) b else F.zero));
        }
    }

    a.divE(b, a) catch unreachable; //since inplace
    for (0..n) |i| {
        for (0..n) |j| {
            try testing.expect(a.at(i, j).cmp(.eq, if (i == j) F.eye else F.zero));
        }
    }
}

test "matrix mulpiplication with vector" {
    const ally = std.testing.allocator;
    const F = @import("field.zig").Float(f32);
    const V = @import("vector.zig").VectorType(F);
    const M = MatrixType(F, usize);

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

// test "matrix LU solve" {
//     const ally = std.testing.allocator;
//     const F = @import("field.zig").Float(f32);
//     const M = MatrixType(F, usize);

//     const n = 3;
//     var a = try M.init(n, n, ally);
//     defer a.deinit();
//     // 2  1 3
//     // 6  0 5
//     // 8 10 7
//     // det = 40+180-100-42 = 78

//     try a.set(0, 0, F.from(2, 1));
//     try a.set(0, 1, F.from(1, 1));
//     try a.set(0, 2, F.from(3, 1));
//     try a.set(1, 0, F.from(6, 1));
//     try a.set(1, 1, F.zero);
//     try a.set(1, 2, F.from(5, 1));
//     try a.set(2, 0, F.from(8, 1));
//     try a.set(2, 1, F.from(10, 1));
//     try a.set(2, 2, F.from(7, 1));

//     var b = try M.Vector.init(n, ally);
//     defer b.deinit(ally);
//     b.set(0, F.from(1, 1));
//     b.set(1, F.from(2, 1));
//     b.set(2, F.from(3, 1));

//     //decompose
//     const lu = try a.decomp(ally);
//     defer lu.deinit(ally);

//     //solve
//     var x = try M.Vector.init(n, ally);
//     defer x.deinit(ally);
//     try lu.solve(b, x);

//     //check solve
//     var b_ = try M.Vector.init(n, ally);
//     defer b_.deinit(ally);
//     a.mulV(x, b_);

//     try testing.expect(b.at(0).cmp(.eq, b_.at(0)));
//     try testing.expect(b.at(1).cmp(.eq, b_.at(1)));
//     try testing.expect(b.at(2).cmp(.eq, b_.at(2)));

//     //check det
//     try testing.expect(F.from(78, 1).cmp(.eq, lu.det()));

//     //reconstuct a from LU
//     const l = try M.init(n, n, ally);
//     defer l.deinit();
//     try lu.lt.transpose(l);
//     for (0..3) |i| {
//         try l.set(i, i, F.eye);
//     }

//     var a_ = try M.init(n, n, ally);
//     defer a_.deinit();
//     try l.mul(lu.u, a_);

//     try testing.expect(a.at(0, 0).cmp(.eq, a_.at(0, 0)));
//     try testing.expect(a.at(0, 1).cmp(.eq, a_.at(0, 1)));
//     try testing.expect(a.at(0, 2).cmp(.eq, a_.at(0, 2)));
//     try testing.expect(a.at(1, 0).cmp(.eq, a_.at(1, 0)));
//     try testing.expect(a.at(1, 1).cmp(.eq, a_.at(1, 1)));
//     try testing.expect(a.at(1, 2).cmp(.eq, a_.at(1, 2)));
//     try testing.expect(a.at(2, 0).cmp(.eq, a_.at(2, 0)));
//     try testing.expect(a.at(2, 1).cmp(.eq, a_.at(2, 1)));
//     try testing.expect(a.at(2, 2).cmp(.eq, a_.at(2, 2)));
// }
