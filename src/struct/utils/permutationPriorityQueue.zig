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

        pub fn deinitToPermutation(self: Self, allocator: Allocator) Permutation {
            allocator.free(self.priorities[0..self.items.len]);
            return self.items;
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
                const old_priority = self.priorities[elem];
                const current_index = self.items.atInv(elem);
                self.priorities[elem] = priority;
                if (before(priority, old_priority)) {
                    self.siftUp(current_index);
                } else if (before(old_priority, priority)) {
                    self.siftDown(current_index);
                }
            } else { //add item
                self.priorities[elem] = priority;
                const current_index = self.items.atInv(elem);
                const insert_index = self.capacity() - 1 - self.len;
                self.items.swap(current_index, insert_index);
                self.len += 1;
                self.siftUp(insert_index);
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
                self.items.swap(child_index, parent_index);
                child_index = parent_index;
            }
        }

        fn siftDown(self: *Self, start_index: usize) void {
            var index = start_index;
            while (true) {
                var first_index = index;
                var first_priority = self.priorities[self.items.at(first_index)];
                //right
                const right_pos = (self.capacity() - index) << 1;
                if (right_pos <= self.len) {
                    const right_index = self.capacity() - right_pos;
                    const right_priority = self.priorities[self.items.at(right_index)];
                    if (before(right_priority, first_priority)) {
                        first_index = right_index;
                        first_priority = right_priority;
                    }
                } else return;
                //left
                const left_pos = right_pos + 1;
                if (left_pos <= self.len) {
                    const left_index = self.capacity() - left_pos;
                    const left_priority = self.priorities[self.items.at(left_index)];
                    if (before(left_priority, first_priority)) {
                        first_index = left_index;
                        first_priority = left_priority;
                    }
                }

                if (first_index == index) return;

                self.items.swap(first_index, index);
                index = first_index;
            }
        }

        pub fn peek(self: *Self) usize {
            assert(self.len > 0);
            return self.items.at(self.capacity() - 1);
        }

        pub fn remove(self: *Self) usize {
            assert(self.len > 0);
            const first_index = self.capacity() - 1;
            const first = self.items.at(first_index);
            self.items.swap(first_index, self.capacity() - self.len);
            self.len -= 1;
            self.siftDown(first_index);
            return first;
        }
    };
}

test "queueing elements" {
    const ally = testing.allocator;
    const P = struct {
        fn lt(a: f32, b: f32) bool {
            return a < b;
        }
    };
    const PPQ = PermutationPriorityQueue(f32, P.lt);
    var queue = try PPQ.init(6, ally);
    //defer queue.deinit(ally);

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

    queue.set(3, 1.5);
    queue.set(2, -0.5);
    queue.set(4, 4);

    try testing.expectEqual(@as(usize, 2), queue.remove());
    try testing.expectEqual(@as(usize, 0), queue.remove());
    try testing.expectEqual(@as(usize, 1), queue.remove());
    try testing.expectEqual(@as(usize, 3), queue.remove());
    try testing.expectEqual(@as(usize, 4), queue.remove());

    const order = queue.deinitToPermutation(ally);
    defer order.deinit(ally);
    try testing.expectEqual(@as(usize, 5), order.at(0));
    try testing.expectEqual(@as(usize, 2), order.at(1));
    try testing.expectEqual(@as(usize, 0), order.at(2));
    try testing.expectEqual(@as(usize, 1), order.at(3));
    try testing.expectEqual(@as(usize, 3), order.at(4));
    try testing.expectEqual(@as(usize, 4), order.at(5));
}
