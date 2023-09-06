# GridapP4est

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gridap.github.io/GridapP4est.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gridap.github.io/GridapP4est.jl/dev)
[![Build Status](https://github.com/gridap/GridapP4est.jl/workflows/CI/badge.svg)](https://github.com/gridap/GridapP4est.jl/actions)
[![Coverage](https://codecov.io/gh/gridap/GridapP4est.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/gridap/GridapP4est.jl)

| ![](https://user-images.githubusercontent.com/38347633/134634010-2be9b499-201b-4166-80ac-e161f6adceb0.png)   |  ![](https://user-images.githubusercontent.com/38347633/134634023-83f37646-f6b9-435c-9f9f-291dea9f86c2.png) 
|:-------------:|:-------------:|

## Purpose 

The purpose of this package is to provide a `DistributedDiscreteModel` implementation (a parallel mesh data structure, see `GridapDistributed.jl` for more details) able to handle forests of quadtrees/octrees of the computational domain. To this end, it leverages the [`p4est` software library](https://p4est.github.io/) meshing engine under the hood (click [here](https://github.com/gridap/GridapP4est.jl/blob/main/test/PoissonUniformlyRefinedOctreeModelsTests.jl) for an example).

## Build 

Before using `GridapP4est.jl`, we have to build `MPI.jl` and 
`P4est_wrapper.jl`. We refer to the main `README.md` of the latter package (available [here](https://github.com/gridap/p4est_wrapper.jl)) for configuration instructions.
