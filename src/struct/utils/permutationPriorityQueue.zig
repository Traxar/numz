const std = @import("std");
const assert = std.debug.assert;
const testing = std.testing;
const Allocator = std.mem.Allocator;
const Permutation = @import("../../struct.zig").Permutation;

pub fn PermutationPriorityQueue(comptime Priority: type, comptime before: fn (a: Priority, b: Priority) bool) type {
    return struct {
        const Self = @This();

        items: Permutation,
        priorities: [*]Priority,
        len: usize, //current queue length

        ///Initialize and return priority queue
        pub fn init(n: usize, allocator: Allocator) !Self {
            return Self{
                .items = try Permutation.eye(n, allocator),
                .priorities = (try allocator.alloc(Priority, n)).ptr,
                .len = 0,
            };
        }

        pub fn deinit(self: Self, allocator: Allocator) void {
            allocator.free(self.priorities[0..self.items.len]);
            self.items.deinit(allocator);
        }

        pub fn capacity(self: Self) usize {
            return self.items.len;
        }

        //returns true if the the element is currently in the queue
        pub fn isQueued(self: Self, elem: usize) bool {
            return self.items.atInv(elem) >= self.capacity() - self.len;
        }

        pub fn set(self: *Self, elem: usize, priority: Priority) void {
            if (self.isQueued(elem)) { //update item
                unreachable; //TODO implement
            } else { //add item
                self.priorities[elem] = priority;
                const current_pos = self.items.atInv(elem);
                const queue_insert_pos = self.capacity() - 1 - self.len;
                _ = self.items.swap(current_pos, queue_insert_pos);
                self.len += 1;
                self.siftUp(queue_insert_pos);
            }
        }

        fn siftUp(self: *Self, start_index: usize) void {
            var child_index = start_index;
            const child = self.items.at(child_index);
            const child_priority = self.priorities[child];
            while (child_index < self.capacity() - 1) {
                var parent_index = self.capacity() - ((self.capacity() - child_index) >> 1);
                const parent = self.items.at(parent_index);
                const parent_priority = self.priorities[parent];
                if (!before(child_priority, parent_priority)) break;
                _ = self.items.swap(child_index, parent_index);
                child_index = parent_index;
            }
        }

        pub fn peek(self: *Self) usize {
            assert(self.len > 0);
            return self.items.at(self.capacity() - 1);
        }
    };
}

test "queueing elements" {
    const ally = testing.allocator;
    const PPQ = PermutationPriorityQueue(f32, struct {
        fn f(a: f32, b: f32) bool {
            return a < b;
        }
    }.f);
    var queue = try PPQ.init(6, ally);
    defer queue.deinit(ally);

    queue.set(1, 1);
    try testing.expectEqual(@as(usize, 1), queue.peek());
    queue.set(2, 2);
    try testing.expectEqual(@as(usize, 1), queue.peek());
    queue.set(0, 0);
    try testing.expectEqual(@as(usize, 0), queue.peek());
    queue.set(3, 3);
    try testing.expectEqual(@as(usize, 0), queue.peek());
    queue.set(4, -1);
    try testing.expectEqual(@as(usize, 4), queue.peek());
    try testing.expect(queue.isQueued(0));
    try testing.expect(queue.isQueued(1));
    try testing.expect(queue.isQueued(2));
    try testing.expect(queue.isQueued(3));
    try testing.expect(queue.isQueued(4));
    try testing.expect(!queue.isQueued(5));
}
