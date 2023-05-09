pub const Float = @import("scalar/float.zig").Float;
pub const Relation = @import("scalar/relation.zig").Relation;

test "scalars" {
    _ = Float;
    _ = Relation;
}
