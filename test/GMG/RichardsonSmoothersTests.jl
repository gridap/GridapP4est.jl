#module RichardsonSmoothersTests
   using GridapP4est
   using Gridap
   smooth=RichardsonSmoother(JacobiLinearSolver(),10)
   A=rand(3,3)
   A=A'*A
   x=rand(3)
   b=rand(3)
   r=b-A*x
   solve!(x,smooth,A,r)
#end
