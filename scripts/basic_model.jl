include("../src/read.jl")
include("../src/constraints.jl")
include("../src/model.jl")
include("../src/active.jl")
include("../src/partitions.jl")
include("../src/problem.jl")
include("../src/utils.jl")


function main()

    # Define run parameters.
    max_iters = 2
    out = stdout
    intervals = [0, 200]
    num_workers = 2
    cost_backup = 1
    M = 1e5

    # Load jobs and create problem data.
    write_log(out, 0, "Loading problem.")
    jobs = build_jobs("../data/actual/runway_1_shift_57.csv", intervals)
    problem = Problem(jobs, num_workers, cost_backup, M)

    # Box and norm constraints on u.
    write_log(out, 0, "Defining base uncertainty.")
    base = ConstraintSet()
    add_constraint!(base, GeneralConstraint(u -> u - 50))
    add_constraint!(base, GeneralConstraint(u -> -u - 50))
    add_constraint!(base, NormConstraint(1, 200))

    # No other constraints in the single partition to start with.
    partitions = [base]

    # Initialise termination criterion variables.
    upper = Inf
    lower = -Inf
    tol = 0.05

    # Start the partition and bound procedure.
    write_log(out, 0, "Starting main loop.")
    for i in 1:max_iters

        # Solve the robust model with partitioned uncertainty.
        write_log(out, i, "Solving partition model with $(length(partitions)) partition(s).")
        solved = solve_partition_model(problem, partitions)

        # Get upper bound from solution.
        upper = getobjectivevalue(solved)
        write_log(out, i, "Upper bound obtained: $upper")

        # For each current partition, get a set of active scenarios.
        write_log(out, i, "Computing active scenarios.")
        points = active_glouton(solved, problem, partitions)

        # Solve deterministic model to obtain lower bound.
        write_log(out, i, "Solving lower bound model with $(length(points)) points.")
        lower_solved = solve_lower_bound_model(problem, points[1])

        # Get lower bound from solution.
        lower = getobjectivevalue(lower_solved)
        write_log(out, i, "Upper bound obtained: $lower")

        # For each current partition and corresponding set of active scenarios...
        write_log(out, i, "Building new partitions.")
        new_partitions = ConstraintSet[]
        for (point_set, partition) in zip(points, partitions)

            # Get the new set of partition constraints, one for each subpartition to
            # be created within the current partition.
            subpartitions = voronoi(point_set)
            for subpartition in subpartitions

                # Add the new constraints to the old ones.
                new_partition = deepcopy(partition)
                for c in subpartition.constraints
                    add_constraint!(new_partition, c)
                end

                push!(new_partitions, new_partition)
            end
        end

        # Update the partitions for next iteration.
        partitions = new_partitions

        # Check bound gap relative to tolerance, as the exit criterion.
        gap = (upper - lower) / abs(lower)
        write_log(out, i, "Relative bound gap is $gap.")
        if gap < tol
            write_log(out, i, "Bound gap below tolerance. Exiting.")
            break
        end
    end

    write_log(out, 0, "Procedure finished")
end

main()
