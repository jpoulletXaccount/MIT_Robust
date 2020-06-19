using CSV
using Dates

include("../../src/read.jl")
include("../../src/constraints.jl")
include("../../src/model.jl")
include("../../src/active.jl")
include("../../src/partitions.jl")
include("../../src/problem.jl")
include("../../src/utils.jl")
include("../../src/uncertainty.jl")
include("../../src/deterministic.jl")
include("../../src/graph_model.jl")


"""
    Objective

"""

function main()

    # Get CL arg.
    possible = [50, 100, 200, 300, 400, 500]
    #num_jobs = possible[parse(Int, ARGS[1])]
    num_jobs = 100
    # Logging.
    out = stdout

    # Define run parameters.
    cost_backup = 1
    M = 1e5
    intervals = [0, 200]

    # Initialise termination criterion variables.
    max_iters = 8
    upper = Inf
    lower = -Inf
    tol = 0.05

    # Initialise array to store all points and the associated links matrices
    # used in computing the lower bounds.
    lb_points = []
    lb_link_matrices = []

    # Load jobs and create problem data.
    write_log(out, 0, "Loading problem.")
    jobs = build_jobs(joinpath(@__DIR__, "../../data/actual/spread_$num_jobs.csv"), intervals)

    # Get deterministic number of workers.
    first_problem = Problem(jobs, 0, 0, 0)
    num_workers = solve_deterministic_model(first_problem)
    write_log(out, 0, "# workers needed to solve the deterministic problem is $num_workers.")

    # Create actual problem.
    problem = Problem(jobs, num_workers, cost_backup, M)

    # Standard constraints on u.
    write_log(out, 0, "Defining base uncertainty.")
    base = standard_uncertainty(problem)

    # No other constraints in the single partition to start with.
    partitions = [base]

    # Get the initial link matrix.
    write_log(out, 0, "Building initial link matrix.")
    link_matrices = [initial_link_matrix(problem, base)]

    # Start the partition and bound procedure.
    write_log(out, 0, "Starting main loop.")
    for i in 1:max_iters

        # Record time for measuring iteration duration.
        start = now()

        # Solve the graph model with partitioned uncertainty.
        write_log(out, i, "Solving graph model with $(length(partitions)) partition(s).")
        solved = solve_graph_model(problem, link_matrices, (NEVER, MAYBE), true)

        # Get upper bound from solution.
        upper = getobjectivevalue(solved)
        write_log(out, i, "Upper bound obtained: $upper")

        # For each current partition and corresponding set of active scenarios...
        write_log(out, i, "Building new partitions and link matrices.")
        partitions, link_matrices, lb_points, lb_link_matrices, Lb =
            build_new_partitions_link_matrices_relax_MAYBE_worstPartition(problem, solved, partitions, link_matrices, lb_points, lb_link_matrices)

        write_log(out, i, "$(length(partitions)) partitions built.")
        write_log(out, i, "$(length(lb_points)) point(s) available to compute lower bound.")
        write_log(out,i, "Lower bound find with relax Maybe $(Lb)")

        # Solve deterministic model to obtain lower bound.
        write_log(out, i, "Solving lower bound model.")

        lower_solved = solve_graph_model(problem, lb_link_matrices, (NEVER, MAYBE), true)
        lower = max(Lb, getobjectivevalue(lower_solved))
        write_log(out, i, "Lower bound obtained: $lower")

        # Check bound gap relative to tolerance, as the exit criterion.
        gap = (upper - lower) / abs(upper)
        write_log(out, i, "Relative bound gap at iteration $i is $gap")

        # Create DataFrame to store results.
        results = DataFrame(iteration = i,
                            num_partitions = length(partitions),
                            upper_bound = upper,
                            num_points = length(lb_points),
                            lower_bound_maybe = Lb,
                            lower_bound = lower,
                            gap = gap,
                            duration = now() - start)

        # Record results in the script's directory.
        write_log(out, 0, "Writing results to file.")
        if i == 1
            CSV.write(joinpath(@__DIR__, "graph_worst_partition_approach_$num_jobs.csv"), results)
        else
            CSV.write(joinpath(@__DIR__, "graph_worst_partition_approach_$num_jobs.csv"), results, append=true)
        end

        if gap < tol
            write_log(out, i, "Bound gap below tolerance. Exiting.")
            break
        end
    end

    write_log(out, 0, "Procedure finished.")
end

main()
