include("../src/read.jl")
include("../src/constraints.jl")
include("../src/model.jl")
include("../src/active.jl")
include("../src/partitions.jl")
include("../src/problem.jl")

jobs = build_jobs("/data/example.csv", [0, 0])
problem = Problem(jobs, 2, 10, 1e5)

# Box and 1-norm constraints on u.
base = ConstraintSet()
add_constraint!(base, GeneralConstraint(u -> u - 100))
add_constraint!(base, GeneralConstraint(u -> -u - 100))
add_constraint!(base, NormConstraint(1, 200))

# No other constraints in the single partition to start with.
partitions = [ConstraintSet()]

solved = solve_partition_model(problem, base, partitions)
println(getobjectivevalue(solved))

points = basic(solved, problem, base, partitions)
partitions = voronoi(points[1])

solved = solve_partition_model(problem, base, partitions)
println(getobjectivevalue(solved))

lwsolved= solve_lower_bound_model(problem,points[1])
println(getobjectivevalue(lwsolved))
