pub const Dense = @import("matrix/dense.zig").DenseType;
pub const Sparse = @import("matrix/sparse.zig").MatrixType;

test "matrices" {
    _ = Dense;
    _ = Sparse;
}
