pub const Permutation = @import("struct/permutation.zig").PermutationType;
pub const Field = @import("struct/field.zig");
pub const Vector = @import("struct/vector.zig").VectorType;
pub const Matrix = @import("struct/matrix.zig");

test "structures" {
    _ = Field;
    _ = Vector;
    _ = Permutation;
    _ = Matrix;
}
