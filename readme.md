# numz

This is a numerics library written in `zig 0.12.0-dev.3336+dbb11915b`

ToDo:

    * inplace sparse matrix addition and elimination (internal optimization, hint: backwards)
    * optimize `rowCast` of dense matrix
    * QR, SVD
    * improve doc comments
    * better test checks for vectors and matrices
    * benchmarks 

Planned Features:

    [x] type of scalars is easily interchangeable.
    [ ] provide truly deterministic decimals.
    [x] support complex numbers
    [x] Vectors
        [x] SIMD support
    [x] Permutations
    [x] Sparse matrices (no SIMD support)
        [x] LU solver
    [x] dense matrices (SIMD support)
        [ ] QR decomp
            [x] MGS
            [ ] Householder
        [ ] LU solvers
    [ ] Meshes
    [ ] PDE toolbox
    [ ] visualization options, for a start a CLI print

Design Choices:

    * functions on structs (like vectors and matrices) only take an allocator as an argument if the result has to be freed. If a struct is supposed to support "inplace" operations (ex: sparse matrix set element) that need allocation, the struct has to store its allocator.
    * functions with allocated results (ex: vector addition) take the result location as an argument.
    * functions with stack results (ex: float addition) return their result.

