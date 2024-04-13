const std = @import("std");
const testing = std.testing;

pub const Major = enum {
    row,
    col,

    pub fn other(a: Major) Major {
        return switch (a) {
            .row => .col,
            .col => .row,
        };
    }
};

pub fn allFieldsOfAToB(comptime A: type, comptime B: type) type {
    const info_A = @typeInfo(A).Struct;
    const info_B = @typeInfo(*B).Pointer;

    var info = info_A;
    var fields: [info_A.fields.len]std.builtin.Type.StructField = undefined;
    inline for (0..fields.len) |i| {
        fields[i] = .{
            .alignment = info_B.alignment,
            .default_value = null,
            .is_comptime = false,
            .name = info_A.fields[i].name,
            .type = B,
        };
    }
    info.fields = fields[0..];
    return @Type(std.builtin.Type{ .Struct = info });
}

test "allFieldsOfAToB" {
    const a = .{ 1, 2.0, true };
    const B = u32;
    const C = allFieldsOfAToB(@TypeOf(a), B);
    const info = @typeInfo(C).Struct;
    inline for (info.fields) |field| {
        try testing.expectEqual(B, field.type);
    }
    try testing.expectEqual(info.fields.len * @sizeOf(B), @sizeOf(C));
}
