pub const Permutation = @import("struct/permutation.zig").PermutationType;
pub const Vector = @import("struct/vector.zig").VectorType;
pub const Matrix = @import("struct/matrix.zig").MatrixType;

test "structures" {
    _ = Vector;
    _ = Permutation;
    _ = Matrix;
}
