const std = @import("std");
const testing = std.testing;

/// outcomes of comparions
pub const Relation = enum {
    less,
    lessEqual,
    equal,
    greaterEqual,
    greater,

    /// returns a => b
    pub fn imp(a: Relation, b: Relation) bool {
        return switch (b) {
            .lessEqual => switch (a) {
                .less, .lessEqual, .equal => true,
                else => false,
            },
            .greaterEqual => switch (a) {
                .greater, .greaterEqual, .equal => true,
                else => false,
            },
            else => a == b,
        };
    }
};

test "implications" {
    try testing.expect(Relation.less.imp(.less));
    try testing.expect(Relation.less.imp(.lessEqual));
    try testing.expect(!Relation.less.imp(.equal));
    try testing.expect(!Relation.less.imp(.greaterEqual));
    try testing.expect(!Relation.less.imp(.greater));

    try testing.expect(!Relation.lessEqual.imp(.less));
    try testing.expect(Relation.lessEqual.imp(.lessEqual));
    try testing.expect(!Relation.lessEqual.imp(.equal));
    try testing.expect(!Relation.lessEqual.imp(.greaterEqual));
    try testing.expect(!Relation.lessEqual.imp(.greater));

    try testing.expect(!Relation.equal.imp(.less));
    try testing.expect(Relation.equal.imp(.lessEqual));
    try testing.expect(Relation.equal.imp(.equal));
    try testing.expect(Relation.equal.imp(.greaterEqual));
    try testing.expect(!Relation.equal.imp(.greater));

    try testing.expect(!Relation.greaterEqual.imp(.less));
    try testing.expect(!Relation.greaterEqual.imp(.lessEqual));
    try testing.expect(!Relation.greaterEqual.imp(.equal));
    try testing.expect(Relation.greaterEqual.imp(.greaterEqual));
    try testing.expect(!Relation.greaterEqual.imp(.greater));

    try testing.expect(!Relation.greater.imp(.less));
    try testing.expect(!Relation.greater.imp(.lessEqual));
    try testing.expect(!Relation.greater.imp(.equal));
    try testing.expect(Relation.greater.imp(.greaterEqual));
    try testing.expect(Relation.greater.imp(.greater));
}
