pub const Permutation = @import("struct/permutation.zig").Permutation;
pub const Scalars = @import("struct/scalar.zig");
pub const Vector = @import("struct/vector.zig").VectorType;
pub const Matrix = @import("struct/matrix.zig").MatrixType;

test "structures" {
    _ = Scalars;
    _ = Vector;
    _ = Permutation;
    //_ = Matrix;
}
