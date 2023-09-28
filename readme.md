# numz

This is a numerics library written in zig.

ToDo:

    * improve comments and doc comments
    * better test checks for vectors and matrices
    * benchmarks 

Planned Features:

    [x] type of scalars is easily interchangeable.
    [ ] provide truly deterministic decimals.
    [ ] support complex numbers
    [x] Vectors
        [x] SIMD support
    [x] Permutations
    [.] Sparse matrices (no SIMD support)
        [.] LU solver
    [ ] dense matrices (SIMD support)
        [ ] solver
    [ ] Meshes
    [ ] PDE toolbox
    [ ] visualization options, for a start a CLI print

Design Choices:

    * functions on structs (like vectors and matrices) only take an allocator as an argument if the result has to be freed. If a struct is supposed to support "inplace" operations (ex: sparse matrix set element) that need allocation, the struct has to store its allocator.
    * functions with allocated results (ex: vector addition) take the result location as an argument.
    * functions with stack results (ex: float addition) return their result.

