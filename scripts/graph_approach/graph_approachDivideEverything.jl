using CSV
using TickTock

include("../../src/read.jl")
include("../../src/constraints.jl")
include("../../src/model.jl")
include("../../src/active.jl")
include("../../src/partitions.jl")
include("../../src/problem.jl")
include("../../src/utils.jl")
include("../../src/uncertainty.jl")

"""
    Objective

On the smallest problem size, show how the number of partitions and variables grows on each
iteration when using the graph solving approach.
"""


function main()
    rm("graph_DivideEverything_Approach.csv",force=true)
    # Logging.
    out = stdout

    # Define run parameters.
    cost_backup = 1
    M = 1e5
    intervals = [0, 800]

    # Initialise termination criterion variables.
    max_iters = 20
    upper = Inf
    lower = -Inf
    tol = 0.05

    # Initialise array to store all points and the associated links matrices
    # used in computing the lower bounds.
    lb_points = []
    lb_link_matrices = []

    # Load jobs and create problem data.
    write_log(out, 0, "Loading problem.")
    jobs = build_jobs("/data/actual/runway_1_shift_57.csv", intervals)
    problemInterm = Problem(jobs, -1, 10, 1e5)
    num_worker = solve_deterministic_model(problemInterm)
    write_log(out, 0, "Nb workers need to solve the deterministic problem " * string(num_worker))
    problem = Problem(jobs, num_worker, cost_backup, M)

    # Standard constraints on u.
    write_log(out, 0, "Defining base uncertainty.")
    base = standard_uncertainty(problem)

    # No other constraints in the single partition to start with.
    partitions = [base]

    # Get the initial link matrix.
    write_log(out, 0, "Building initial link matrix.")
    link_matrices = [initial_link_matrix(problem, base)]

    # Create DataFrame to store results.
    results = DataFrame(iteration = [],
                        partitions = [],
                        upper_bound = [],
                        points = [],
                        lower_bound_maybe = [],
                        lower_bound = [],
                        gap = [])


    # Start the partition and bound procedure.
    write_log(out, 0, "Starting main loop.")
    tick()
    for i in 1:max_iters

        # Solve the graph model with partitioned uncertainty.
        write_log(out, i, "Solving graph model with $(length(link_matrices)) partition(s).")
        solved = solve_graph_model(problem, link_matrices, (NEVER, MAYBE),true)

        # Get upper bound from solution.
        upper = getobjectivevalue(solved)
        write_log(out, i, "Upper bound obtained: $upper")

        # For each current partition and corresponding set of active scenarios...
        write_log(out, i, "Building new partitions and link matrices.")
        partitions, link_matrices, lb_points, lb_link_matrices, Lb =
            build_new_partitions_link_matrices_relax_MAYBE(problem, solved, partitions, link_matrices, lb_points,lb_link_matrices)

        write_log(out, i, "$(length(partitions)) partitions built.")
        write_log(out, i, "$(length(lb_points)) points available to compute lower bound.")
        write_log(out,i, "Lower bound find with relax Maybe $(Lb)")

        # Solve deterministic model to obtain lower bound.
        write_log(out, i, "Solving lower bound model.")

        lower_solved = solve_graph_model(problem, lb_link_matrices, (NEVER, MAYBE),true)
        lower = max(Lb ,getobjectivevalue(lower_solved))
        write_log(out, i, "Lower bound obtained: $lower")

        # Check bound gap relative to tolerance, as the exit criterion.
        gap = (upper - lower) / abs(upper)
        write_log(out, i, "Relative bound gap at iteration $i is $gap")

        # Collect results.
        push!(results, [i,
                        length(partitions),
                        upper,
                        length(lb_points),
                        Lb,
                        lower,
                        gap])

        if gap < tol
            write_log(out, i, "Bound gap below tolerance. Exiting.")
            break
        end
        laptimer()
        # Record results in the script's directory.
        write_log(out, 0, "Writing results to file.")
        CSV.write("graph_DivideEverything_Approach.csv", results,append= true)
    end

    write_log(out, 0, "Procedure finished.")
end

main()
