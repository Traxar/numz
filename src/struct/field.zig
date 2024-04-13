pub const Float = @import("field/float.zig").FloatType;
pub const Complex = @import("field/complex.zig").ComplexType;

test "fields" {
    _ = Float;
    _ = Complex;
}
