pub const Permutation = @import("struct/permutation.zig").Permutation;
pub const BiPermutation = @import("struct/permutation.zig").BiPermutationType;
pub const Field = @import("struct/field.zig");
pub const Vector = @import("struct/vector.zig").VectorType;
//pub const Matrix = @import("struct/matrix.zig").MatrixType;

test "structures" {
    _ = Field;
    _ = Vector;
    _ = Permutation;
    _ = BiPermutation;
    //_ = Matrix;
}
