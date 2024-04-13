const std = @import("std");
const assert = std.debug.assert;
const panic = std.debug.panic;

const testing = std.testing;
const Allocator = std.mem.Allocator;

const DenseUnmanagedType = @import("base/dense.zig").DenseUnmanagedType;
const Utils = @import("utils.zig");
const Majority = Utils.Major;

/// struct for dense matrix computations
pub fn DenseType(comptime Element: type, comptime majority: Majority) type {
    assert(Element.simd_size == 1);
    const DenseUnmanaged = DenseUnmanagedType(Element);
    const SimdElement = DenseUnmanaged.SimdElement;
    const simd_size = SimdElement.simd_size;
    return struct {
        const Dense = @This();
        const DenseT = DenseType(Element, majority.other());
        const major = majority;
        //const Sparse = @import("sparse.zig").SparseType(Element, major);
        val: DenseUnmanaged,
        rows: usize,
        cols: usize,

        fn n(a: Dense) usize {
            return switch (majority) {
                .row => a.rows,
                .col => a.cols,
            };
        }

        fn m(a: Dense) usize {
            return switch (majority) {
                .row => a.cols,
                .col => a.rows,
            };
        }

        pub fn init(rows: usize, cols: usize, allocator: Allocator) !Dense {
            var res = Dense{
                .val = undefined,
                .rows = rows,
                .cols = cols,
            };
            res.val = try DenseUnmanaged.init(res.n(), res.m(), allocator);
            return res;
        }

        pub fn deinit(a: Dense, allocator: Allocator) void {
            a.val.deinit(a.n(), allocator);
        }

        fn indAt(a: Dense, i: usize, j: usize) DenseUnmanaged.Index {
            assert(i <= a.rows and j <= a.cols);
            return switch (majority) {
                .row => a.val.indAt(i, j),
                .col => a.val.indAt(j, i),
            };
        }

        pub fn at(a: Dense, i: usize, j: usize) Element {
            assert(i < a.rows and j < a.cols);
            const ind = a.indAt(i, j);
            return a.val.at(.j_one, ind);
        }

        pub fn set(a: Dense, i: usize, j: usize, b: Element) void {
            assert(i < a.rows and j < a.cols);
            const ind = a.indAt(i, j);
            a.val.set(.j_one, ind, b);
        }

        /// casts the matrix to its transpose but in the other majority
        pub fn t(a: Dense) DenseT {
            return .{
                .val = a.val,
                .rows = a.cols,
                .cols = a.rows,
            };
        }

        /// res <- a^T
        pub fn transpose(res: Dense, a: Dense) void {
            assert(res.rows == a.cols and res.cols == a.rows);
            const inplace = res.val.val == a.val.val;
            var iter = res.indAt(res.rows, res.cols);
            var iter_ = a.indAt(a.rows, a.cols);
            const iter_a = if (inplace) &iter else &iter_;
            while (iter.prev(.i_one, res.val)) {
                assert(iter_a.prev(.j_one, a.val));
                var ind = iter;
                var ind_a = iter_a.*;
                while (ind.prev(.j_one, res.val)) {
                    assert(ind_a.prev(.i_one, a.val));
                    if (inplace) {
                        const h = res.val.at(.j_one, ind);
                        res.val.set(.j_one, ind, a.val.at(.j_one, ind_a));
                        a.val.set(.j_one, ind_a, h);
                    } else {
                        res.val.set(.j_one, ind, a.val.at(.j_one, ind_a));
                    }
                }
            }
        }

        inline fn argsMajority(comptime fn_type: type, comptime args_type: type) ?Majority {
            comptime {
                const args_info = @typeInfo(args_type).Struct;
                const fn_info = @typeInfo(fn_type).Fn;
                if (fn_info.params.len != args_info.fields.len) @compileError("number of inputs do not match");
                var maj_: ?Majority = null;
                for (args_info.fields) |arg| {
                    switch (arg.type) {
                        Element => {},
                        Dense, DenseT => {
                            if (maj_) |maj| {
                                if (maj != arg.type.major) @compileError("matrices must have same majority");
                            } else maj_ = arg.type.major;
                        },
                        else => @compileError("unsupported type"),
                    }
                }
                return maj_;
            }
        }

        /// asserts that all argument dimesions are the same and returns them in the form of a matrix with no values
        fn argsDimensions(args: anytype) Dense {
            const args_info = @typeInfo(@TypeOf(args)).Struct;
            var r: ?usize = null;
            var c: ?usize = null;
            inline for (args_info.fields) |arg| {
                if (arg.type != Element) {
                    const a = @field(args, arg.name);
                    if (r == null and c == null) {
                        r = a.rows;
                        c = a.cols;
                    } else if (a.rows != r or a.cols != c) panic(
                        "expected argument dimensions: {any}x{any}, found: {}x{}",
                        .{ r, c, a.rows, a.cols },
                    );
                }
            }
            if (r == null or c == null) panic("expected atleast one matrix argument, found none", .{});
            var res = Dense{ .rows = r.?, .cols = c.?, .val = undefined };
            res.val.m_simd = 1 + @divFloor(res.m() - 1, simd_size);
            return res;
        }

        fn Args(comptime step: DenseUnmanaged.Index.Step, args_type: type) type {
            return Utils.allFieldsOfAToB(
                args_type,
                switch (step) {
                    .j_simd => SimdElement,
                    else => Element,
                },
            );
        }

        fn prepArgs(comptime step: DenseUnmanaged.Index.Step, args: anytype) Args(step, @TypeOf(args)) {
            const args_info = @typeInfo(@TypeOf(args)).Struct;
            var a: Args(step, @TypeOf(args)) = undefined;
            inline for (args_info.fields) |arg| {
                if (arg.type == Element) {
                    @field(a, arg.name) = switch (step) {
                        .j_simd => SimdElement.simdSplat(@field(args, arg.name)),
                        else => @field(args, arg.name),
                    };
                }
            }
            return a;
        }

        fn setArgs(comptime step: DenseUnmanaged.Index.Step, ind: DenseUnmanaged.Index, args: anytype, a: *Args(step, @TypeOf(args))) void {
            const args_info = @typeInfo(@TypeOf(args)).Struct;
            inline for (args_info.fields) |arg| {
                if (arg.type != Element) {
                    @field(a, arg.name) = @field(args, arg.name).val.at(step, ind);
                }
            }
        }

        // res <- op(args)
        pub fn ew(res: Dense, comptime op: Element.Operator, args: anytype) if (op.ErrorSet()) |e| e!void else void {
            if (argsMajority(@TypeOf(op.f()), @TypeOf(args))) |maj| {
                if (maj != majority) @compileError("argument majority must match result");
                const size = argsDimensions(args);
                if (size.rows != res.rows or size.cols != res.cols) panic(
                    "expected argument dimensions: {}x{}, found: {}x{}",
                    .{ res.rows, res.cols, size.rows, size.cols },
                );
                const simd_op: SimdElement.Operator = @enumFromInt(@intFromEnum(op));
                const ops = .{ op, simd_op };
                const steps = .{ .j_one, .j_simd };
                var args_ = .{ prepArgs(.j_one, args), prepArgs(.j_simd, args) };
                var iter = size.indAt(size.rows, size.cols);
                if (op.ErrorSet() == null) iter.expandToSimd();
                while (iter.prev(.i_one, size.val)) {
                    var ind = iter;
                    inline for (steps, &args_, ops) |s, *a, o| {
                        while (true) {
                            if (simd_size > 1 and s == .j_one) {
                                if (ind.sub == 0) break;
                                assert(ind.prev(s, size.val));
                            } else if (!ind.prev(s, size.val)) break;
                            setArgs(s, ind, args, a);
                            const r_ = @call(.always_inline, o.f(), a.*);
                            const r = if (o.ErrorSet()) |_| try r_ else r_;
                            res.val.set(s, ind, r);
                        }
                    }
                }
            } else {
                const r_ = @call(.always_inline, op.f(), args);
                const r = if (op.ErrorSet()) |_| try r_ else r_;
                res.val.fill(res.n(), r);
            }
        }

        /// <- op_red(op_ew(args))
        pub fn red(comptime op_red: Element.Operator, comptime op_ew: Element.Operator, args: anytype) if (op_ew.ErrorSet()) |e| e!Element else Element {
            if (op_red.ErrorSet()) |_| @compileError("error on reduction operator not supported");
            if (argsMajority(@TypeOf(op_ew.f()), @TypeOf(args))) |_| {
                const size = argsDimensions(args);
                const simd_op_red: SimdElement.Operator = @enumFromInt(@intFromEnum(op_red));
                const simd_op_ew: SimdElement.Operator = @enumFromInt(@intFromEnum(op_ew));
                const steps = .{ .j_one, .j_simd };
                const ops_red = .{ op_red, simd_op_red };
                const ops_ew = .{ op_ew, simd_op_ew };
                var args_ = .{ prepArgs(.j_one, args), prepArgs(.j_simd, args) };
                var sum_one: ?Element = null;
                var sum_simd: ?SimdElement = null;
                const sums = .{ &sum_one, &sum_simd };
                var iter = size.indAt(size.rows, size.cols);
                while (iter.prev(.i_one, size.val)) {
                    var ind = iter;
                    inline for (steps, &args_, ops_red, ops_ew, sums) |s, *a, red_, ew_, sum_| {
                        while (true) {
                            if (s == .j_one) {
                                if (ind.subIndex() == 0) break;
                                assert(ind.prev(s, size.val));
                            } else if (!ind.prev(s, size.val)) break;
                            setArgs(s, ind, args, a);
                            const r_ = @call(.always_inline, ew_.f(), a.*);
                            const r = if (ew_.ErrorSet()) |_| try r_ else r_;
                            sum_.* = if (sum_.*) |sum| @call(.always_inline, red_.f(), .{ sum, r }) else r;
                        }
                    }
                }
                if (sum_simd) |s_simd| {
                    const sum = s_simd.simdReduce(simd_op_red);
                    return if (sum_one) |s_ones| @call(.always_inline, op_red.f(), .{ sum, s_ones }) else sum;
                } else {
                    return sum_one.?;
                }
            } else {
                @compileError("expected atleast one matrix argument");
            }
        }

        // <- and(cp(args))
        pub fn cmp(comptime cp: Element.Comparator, args: anytype) bool {
            if (argsMajority(@TypeOf(cp.f()), @TypeOf(args))) |_| {
                const size = argsDimensions(args);
                const simd_cp: SimdElement.Comparator = @enumFromInt(@intFromEnum(cp));
                const steps = .{ .j_one, .j_simd };
                const cps = .{ cp, simd_cp };
                var args_ = .{ prepArgs(.j_one, args), prepArgs(.j_simd, args) };
                const Es = .{ Element, SimdElement };
                var iter = size.indAt(size.rows, size.cols);
                while (iter.prev(.i_one, size.val)) {
                    var ind = iter;
                    inline for (steps, &args_, cps, Es) |s, *a, c, E| {
                        while (true) {
                            if (s == .j_one) {
                                if (ind.subIndex() == 0) break;
                                assert(ind.prev(s, size.val));
                            } else if (!ind.prev(s, size.val)) break;
                            setArgs(s, ind, args, a);
                            const r = @call(.always_inline, c.f(), a.*);
                            if (!E.all(r)) return false;
                        }
                    }
                }
                return true;
            } else {
                return @call(.always_inline, cp.f(), args);
            }
        }

        /// res <- a * b
        pub fn mul(res: Dense, a: anytype, b: anytype) void {
            if (@TypeOf(res) == DenseType(Element, .col)) return res.t().mul(b.t(), a.t());
            if (a.cols != b.rows or a.rows != res.rows or b.cols != res.cols) panic("argument dimensions not fit for multiplication", .{});
            const a_is_row = @TypeOf(a) == Dense;
            const b_is_row = @TypeOf(b) == Dense;
            const i_step = .i_one;
            const i_step_a = if (a_is_row) .i_one else .j_one;
            const j_step = if (b_is_row) .j_simd else .j_one;
            const j_step_b = if (b_is_row) .j_simd else .i_one;
            const k_step_a = if (!a_is_row) .i_one else if (b_is_row) .j_one else .j_simd;
            const k_step_b = if (b_is_row) .i_one else if (a_is_row) .j_simd else .j_one;
            var iter = res.indAt(res.rows, res.cols);
            if (b_is_row) iter.expandToSimd();
            var iter_a = a.indAt(a.rows, a.cols);
            var init_b = b.indAt(b.rows, b.cols);
            if (b_is_row) init_b.expandToSimd();
            while (iter_a.prev(i_step_a, a.val)) { // i
                assert(iter.prev(i_step, res.val));
                var iter_b = init_b;
                var ind = iter;
                while (ind.prev(j_step, res.val)) { // j
                    assert(iter_b.prev(j_step_b, b.val));
                    var ind_a = iter_a;
                    var ind_b = iter_b;
                    var sum_ones = Element.zero;
                    var sum_simd = SimdElement.zero;
                    if (!b_is_row) {
                        while (true) { // k ones
                            if (a_is_row) {
                                if (simd_size > 1 and ind_a.sub > 0) {
                                    assert(ind_a.prev(.j_one, a.val));
                                } else break;
                            } else if (!ind_a.prev(k_step_a, a.val)) break;
                            assert(ind_b.prev(.j_one, b.val));
                            const a_ = a.val.at(.j_one, ind_a);
                            const b_ = b.val.at(.j_one, ind_b);
                            sum_ones = sum_ones.add(a_.mul(b_));
                        }
                    }
                    if (a_is_row or b_is_row) {
                        while (ind_a.prev(k_step_a, a.val)) { // k simd
                            assert(ind_b.prev(k_step_b, b.val));
                            const a__ = a.val.at(k_step_a, ind_a);
                            const a_ = if (@TypeOf(a__) == SimdElement) a__ else SimdElement.simdSplat(a__);
                            const b_ = b.val.at(.j_simd, ind_b);
                            sum_simd = sum_simd.add(a_.mul(b_));
                        }
                    }
                    const sum = if (b_is_row) sum_simd else if (!a_is_row) sum_ones else sum_ones.add(sum_simd.simdReduce(.add));
                    res.val.set(j_step, ind, sum);
                }
            }
        }

        // TODO: implement inverse
        /// res_i <- op_red(op_ew(args_ij))
        pub fn collapse(res: Dense, comptime op_red: Element.Operator, comptime op_ew: Element.Operator, args: anytype) if (op_ew.ErrorSet()) |e| e!void else void {
            if (op_red.ErrorSet()) |_| @compileError("error on reduction operator not supported");
            if (res.rows != 1) panic("expected result to have 1 row, found {}", .{res.rows});
            if (comptime argsMajority(@TypeOf(op_ew.f()), @TypeOf(args))) |maj| {
                const size = DenseType(Element, maj).argsDimensions(args);
                if (size.cols != res.cols) panic(
                    "expected arguments to have {} columns, found {}",
                    .{ res.cols, size.cols },
                );
                const simd_op_red: SimdElement.Operator = @enumFromInt(@intFromEnum(op_red));
                const simd_op_ew: SimdElement.Operator = @enumFromInt(@intFromEnum(op_ew));
                const ops_red = .{ op_red, simd_op_red };
                const ops_ew = .{ op_ew, simd_op_ew };
                const res_is_row = majority == .row;
                const arg_is_row = maj == .row;
                var args_ = .{ prepArgs(.j_one, args), prepArgs(.j_simd, args) };
                var ind = res.indAt(0, res.cols);
                var iter_arg = size.indAt(size.rows, size.cols);
                if (arg_is_row) {
                    const j_step = if (res_is_row) .{ .j_one, .j_simd } else .{ .i_one, .i_one };
                    const j_step_arg = .{ .j_one, .j_simd };
                    const i_step_arg = .i_one;
                    const Es = .{ Element, SimdElement };
                    inline for (j_step, j_step_arg, &args_, ops_red, ops_ew, Es) |s, s_a, *a, red_, ew_, E| {
                        while (ind.prev(s, res.val)) {
                            assert(iter_arg.prev(s_a, size.val));
                            var ind_arg = iter_arg;
                            var sum_: ?E = null;
                            while (ind_arg.prev(i_step_arg, size.val)) {
                                setArgs(s_a, ind_arg, args, a);
                                const r_ = @call(.always_inline, ew_.f(), a.*);
                                const r = if (ew_.ErrorSet()) |_| try r_ else r_;
                                sum_ = if (sum_) |sum| @call(.always_inline, red_.f(), .{ sum, r }) else r;
                            }
                            if (s == s_a) {
                                res.val.set(s, ind, sum_.?);
                            } else {
                                while (true) {
                                    const sub = ind.i - iter_arg.j;
                                    res.val.set(s, ind, sum_.?.simdAt(sub));
                                    if (sub == 0) break;
                                }
                            }
                        }
                    }
                } else {
                    const j_step = if (res_is_row) .j_one else .i_one;
                    const j_step_arg = .i_one;
                    const i_step_arg = .{ .j_one, .j_simd };
                    while (ind.prev(j_step, res.val)) {
                        assert(iter_arg.prev(j_step_arg, size.val));
                        var ind_arg = iter_arg;
                        var sum_one: ?Element = null;
                        var sum_simd: ?SimdElement = null;
                        const sums = .{ &sum_one, &sum_simd };
                        inline for (i_step_arg, &args_, ops_red, ops_ew, sums) |s, *a, red_, ew_, sum_| {
                            while (true) {
                                if (s == .j_one) {
                                    if (ind_arg.subIndex() == 0) break;
                                    assert(ind_arg.prev(s, size.val));
                                } else if (!ind_arg.prev(s, size.val)) break;
                                setArgs(s, ind_arg, args, a);
                                const r_ = @call(.always_inline, ew_.f(), a.*);
                                const r = if (ew_.ErrorSet()) |_| try r_ else r_;
                                sum_.* = if (sum_.*) |sum| @call(.always_inline, red_.f(), .{ sum, r }) else r;
                            }
                        }
                        res.val.set(j_step, ind, sum: {
                            if (sum_simd) |s_simd_| {
                                const s_simd = s_simd_.simdReduce(simd_op_red);
                                break :sum if (sum_one) |s_one| @call(.always_inline, op_red.f(), .{ s_simd, s_one }) else s_simd;
                            } else {
                                break :sum sum_one.?;
                            }
                        });
                    }
                }
            } else {
                @compileError("expected atleast one matrix argument");
            }
        }
    };
}

test "dense matrix elementwise" {
    const ally = std.testing.allocator;
    const F = @import("../field.zig").Float(f32);
    const C = @import("../field.zig").Complex(F);
    const M = DenseType(C, .col);

    const n = 3;
    const m = 5;
    const a = try M.init(n, m, ally);
    defer a.deinit(ally);
    a.ew(.id, .{C.i});
    a.set(2, 4, C.one.add(C.one));
    a.ew(.mul, .{ a, C.i });
    a.ew(.sqrt, .{a.t().t()});
    try testing.expectEqual(C.i, a.at(0, 0));
    try testing.expectEqual(C.one.add(C.i), a.at(2, 4));
    const sum = M.red(.add, .id, .{a});
    try testing.expectEqual(C.from(F.one, F.from(n * m, 1)), sum);
    try testing.expect(!M.cmp(.eq, .{ a, C.i }));
    try testing.expect(!M.cmp(.neq, .{ a, C.i }));
}

test "dense matrix multiplication" {
    const ally = std.testing.allocator;
    const F = @import("../field.zig").Float(f32);
    const majors = .{ .row, .col };
    inline for (majors) |maj_a| {
        const A = DenseType(F, maj_a);
        const a = try A.init(2, 3, ally);
        defer a.deinit(ally);
        for (0..2) |i| {
            for (0..3) |j| {
                a.set(i, j, F.from(@intCast(i * 3 + j), 1));
            }
        }
        inline for (majors) |maj_b| {
            const B = DenseType(F, maj_b);
            const b = try B.init(3, 1, ally);
            defer b.deinit(ally);
            for (0..3) |i| {
                b.set(i, 0, F.from(@intCast(i + 1), 1));
            }
            inline for (majors) |maj_c| {
                const C = DenseType(F, maj_c);
                const c = try C.init(2, 1, ally);
                defer c.deinit(ally);

                c.mul(a, b);

                try testing.expectEqual(F.from(8, 1), c.at(0, 0));
                try testing.expectEqual(F.from(26, 1), c.at(1, 0));
            }
        }
    }
}

test "dense matrix collapse" {
    const ally = std.testing.allocator;
    const F = @import("../field.zig").Float(f32);
    const majors = .{ .row, .col };
    inline for (majors) |maj_a| {
        const A = DenseType(F, maj_a);
        const a = try A.init(2, 3, ally);
        defer a.deinit(ally);
        for (0..2) |i| {
            for (0..3) |j| {
                a.set(i, j, F.from(@intCast(i * 3 + j), 1));
            }
        }
        inline for (majors) |maj_b| {
            const B = DenseType(F, maj_b);
            const b = try B.init(2, 1, ally);
            defer b.deinit(ally);

            try b.t().collapse(.add, .div, .{ a.t(), F.from(2, 1) });

            try testing.expectEqual(F.from(3, 2), b.at(0, 0));
            try testing.expectEqual(F.from(6, 1), b.at(1, 0));
        }
    }
}
