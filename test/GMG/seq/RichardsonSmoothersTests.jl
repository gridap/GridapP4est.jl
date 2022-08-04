module RichardsonSmoothersTests
   using GridapP4est
   using Gridap
   using Test
   using LinearAlgebra
   smooth=RichardsonSmoother(JacobiLinearSolver(),10)
   A=rand(3,3)
   A=A'*A+10.0*I
   x=rand(3)
   b=rand(3)
   r=b-A*x
   r0=norm(r)
   solve!(x,smooth,A,r)
   @test norm(r)/norm(r0) < 1.0e-06
end
