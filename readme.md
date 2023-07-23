# numz
This is a numerics library written in zig.

ToDo:

    * improve comments and doc comments
    * better test checks for vectors and matrices
    * 2nd permutation type without inverse
    * benchmarks 

Planned Features:

    [x] type of scalars is easily interchangeable.
    [ ] provide truly deterministic decimals.
    [x] Vectors
        [ ] SIMD support
    [x] Permutations
    [x] Sparse matrices
        [x] LU solver
    [ ] Meshes
    [ ] PDE toolbox
    [ ] visualization options, for a start a CLI print

Design Choices:

    * functions on structs (like vectors and matrices) only take an allocator as an argument if the result has to be freed.

