# Num
This is a numerics library written in zig.

Planned Features:

    [x] type of scalars is easily interchangeable.
    [ ] provide truly deterministic decimals.
    [x] Vectors
    [ ] Permutations
    [ ] Sparse matrices: including LU solver
    [ ] Meshes
    [ ] PDE toolbox

Design Choices:

    * functions on structs (like vectors and matrices) only take an allocator as an argument if the result has to be freed.

