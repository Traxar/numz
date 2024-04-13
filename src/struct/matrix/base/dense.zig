const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;
const Allocator = std.mem.Allocator;

/// unmanaged base struct for dense matrix computations
pub fn DenseUnmanagedType(comptime Element: type) type {
    assert(Element.simd_size == 1);
    return struct {
        pub const SimdElement = Element.SimdType(null); // null gives prefered size by the Element itself
        const simd_size = SimdElement.simd_size;
        const simd = simd_size > 1;
        const DenseUnmanaged = @This();
        val: [*]SimdElement,
        m_simd: usize,

        /// initialize dense data for a nxm matrix with undefined values
        pub fn init(n: usize, m: usize, allocator: Allocator) !DenseUnmanaged {
            var res: DenseUnmanaged = undefined;
            res.m_simd = 1 + @divFloor(m - 1, simd_size); //divCeil
            res.val = (try allocator.alloc(SimdElement, res.m_simd * n)).ptr;
            return res;
        }

        pub fn deinit(a: DenseUnmanaged, n: usize, allocator: Allocator) void {
            allocator.free(a.val[0 .. n * a.m_simd]);
        }

        pub fn fill(a: DenseUnmanaged, n: usize, b: Element) void {
            @memset(a.val[0 .. n * a.m_simd], SimdElement.simdSplat(b));
        }

        pub fn at(a: DenseUnmanaged, comptime step: Index.Step, ind: Index) step.ElementType() {
            if (simd and !step.isSimd()) {
                return a.val[ind.simd].simdAt(ind.sub);
            } else {
                return a.val[ind.simd];
            }
        }

        pub fn set(a: DenseUnmanaged, comptime step: Index.Step, ind: Index, b: anytype) void {
            if (step.ElementType() != @TypeOf(b)) @compileError("step and type of new value do not match");
            if (simd and !step.isSimd()) {
                a.val[ind.simd].simdSet(ind.sub, b);
            } else {
                a.val[ind.simd] = b;
            }
        }

        pub fn indAt(a: DenseUnmanaged, i: usize, j: usize) Index {
            var res: Index = undefined;
            res.i = i;
            res.j = j;
            if (simd) {
                res.simd = @divFloor(j, simd_size);
                res.sub = j - res.simd * simd_size;
                res.simd += i * a.m_simd;
            } else {
                res.simd = i * a.m_simd + j;
            }
            return res;
        }

        pub fn majorCast(a: DenseUnmanaged, i: usize) DenseUnmanaged {
            return .{
                .val = &a.val[i * a.m_simd],
                .m_simd = a.m_simd,
            };
        }

        pub const Index = struct {
            i: usize,
            j: usize,
            simd: usize,
            sub: if (simd) usize else void = undefined,

            pub const Step = enum {
                i_one,
                j_one,
                j_simd,

                pub inline fn isSimd(comptime step: Step) bool {
                    return switch (step) {
                        .j_simd => true,
                        else => false,
                    };
                }

                pub fn ElementType(comptime step: Step) type {
                    return if (step.isSimd()) SimdElement else Element;
                }
            };

            pub inline fn subIndex(ind: Index) usize {
                return if (simd) ind.sub else 0;
            }

            pub fn expandToSimd(ind: *Index) void {
                if (simd and ind.sub != 0) {
                    ind.j += simd_size - ind.sub;
                    ind.sub = 0;
                    ind.simd += 1;
                }
            }

            pub fn prev(ind: *Index, comptime step: Step, a: DenseUnmanaged) bool {
                assert(step != .j_simd or !simd or ind.sub == 0); // only allow j_simd step if sub index is 0
                // check if prev exists
                switch (step) {
                    .i_one => if (ind.i == 0) return false,
                    else => if (ind.j == 0) return false,
                }
                // calc prev
                const step_eff: Step = if (!simd and step == .j_one) .j_simd else step;
                switch (step_eff) {
                    .i_one => {
                        ind.i -= 1;
                        ind.simd -= a.m_simd;
                    },
                    .j_one => {
                        ind.j -= 1;
                        if (ind.sub == 0) {
                            ind.sub = simd_size - 1;
                            ind.simd -= 1;
                        } else {
                            ind.sub -= 1;
                        }
                    },
                    .j_simd => {
                        ind.j -= simd_size;
                        ind.simd -= 1;
                    },
                }
                return true;
            }
        };
    };
}
