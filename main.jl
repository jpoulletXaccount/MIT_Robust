include("src/read.jl")
include("src/constraints.jl")
include("src/model.jl")
include("src/active.jl")
include("src/partitions.jl")
include("src/problem.jl")
include("src/utils.jl")
include("src/deterministic.jl")

function main()

    MAX_ITER = 2
    out = stdout

    # Load jobs and create problem data.
    write_log(out, 0, "Loading problem.")
    jobs = build_jobs("/data/actual/runway_1_shift_57.csv", [0, 0])
    problemInterm = Problem(jobs, -1, 10, 1e5)
    num_worker = solve_deterministic_model(problemInterm)
    write_log(out, 0, "Nb workers need to solve the deterministic problem " * string(num_worker))
    problem = Problem(jobs, num_worker, 10, 1e5)

    # Box and norm constraints on u.
    base = ConstraintSet()
    add_constraint!(base, GeneralConstraint(u -> u - 50))
    add_constraint!(base, GeneralConstraint(u -> -u - 50))
    add_constraint!(base, NormConstraint(1, 200))

    # No other constraints in the single partition to start with.
    partitions = [base]

    # Initialise termination criterion variables.
    upperB = Inf
    lowerB = -Inf
    tol = 0.05

    # Start the partition and bound procedure.
    write_log(out, 0, "Starting main loop")
    for iter in 1:MAX_ITER

        # Check bound gap relative to tolerance, as the exit criterion.
        gap = (upperB - lowerB) / abs(lowerB)
        write_log(out, iter, "Relative bound gap at iteration $iter is $gap")
        if gap < tol
            write_log(out, iter, "Bound gap below tolerance at iteration $iter")
            break
        end

        # Solve the robust model with partitioned uncertainty.
        write_log(out, iter, "Solving partition model")
        solved = solve_partition_model(problem, partitions)
        upperB = getobjectivevalue(solved)
        println(upperB)

        # For each current partition, get a set of active scenarios.
        write_log(out, iter, "Computing active scenarios")
        points = active_glouton(solved, problem, partitions)

        # Solve deterministic model to obtain lower bound.
        write_log(out, iter, "Solving lower bound model")
        lwsolved = solve_lower_bound_model(problem, points[1])
        lowerB = getobjectivevalue(lwsolved)

        # Initialise list to store the new partitions.
        new_partitions = ConstraintSet[]

        # For each current partition and corresponding set of active scenarios...
        write_log(out, iter, "Building new partitions")
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
    end

    write_log(out, 0, "Procedure finished")
end

main()
