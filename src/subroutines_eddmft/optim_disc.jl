# Import the package and define the problem to optimize
using Optimization


rosenbrock(u, p) = (p[1] - u[1])^2 + p[2] * (u[2] - u[1]^2)^2
u0 = zeros(2)
p = [1.0, 100.0]

prob = OptimizationProblem(rosenbrock, u0, p)

# Import a solver package and solve the optimization problem
using OptimizationOptimJL
sol = solve(prob, NelderMead())

# Import a different solver package and solve the optimization problem a different way
using OptimizationBBO
prob = OptimizationProblem(rosenbrock, u0, p, lb = [-1.0, -1.0], ub = [1.0, 1.0])
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited())

