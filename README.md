#A new function called updateG is implemented.

* updateG(double t, int j) locates at dynamic_system_solver.cc of MBSim kernel.
it calculates the projection direction matrix of the smooth contact force la.

* In the CVT simulation, the original updateG takes a very high percentage of the computational effort,
  this new updateG can reduced the total computational cost up to 40%, due to sparse matrix linear algebra
  techniques are applied.

* csparse library  is used for solving a L x = b, where L is a sparse lower triangular matrix.
