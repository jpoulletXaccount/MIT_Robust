include("src/read.jl")
include("src/constraints.jl")
include("src/graphModeling.jl")
include("src/active.jl")
include("src/partitions.jl")
include("src/problem.jl")
include("src/utils.jl")


function mainGraph()

    MAX_ITER = 5
    out = stdout

    # Load jobs and create problem data.
    jobs = build_jobs("/data/runway_1_shift_57.csv", [0, 200])
    problem = Problem(jobs, 10, 10, 1e5)

    # Box and norm constraints on u.
    base = ConstraintSet()
    add_constraint!(base, GeneralConstraint(u -> u - 50))
    add_constraint!(base, GeneralConstraint(u -> -u - 50))
    add_constraint!(base, NormConstraint(1, 100))

    # No other constraints in the single partition to start with.
    partitions = [base]

    # Initialise termination criterion variables.
    upperB = Inf
    lowerB = -Inf
    tol = 0.05
    finalListPoint = []

    # Start the partition and bound procedure.
    logUtils(out, 0, "Starting main loop")
    for iter in 1:MAX_ITER

        # Solve the robust model with partitioned uncertainty.
        logUtils(out, iter, "Solving partition model")
        solved = create_graph_model(problem, partitions)
        upperB = getobjectivevalue(solved)
        logUtils(out, iter, "upper bound obtained " * string(upperB))

        # Initialise list to store the new partitions.
        new_partitions = ConstraintSet[]

        # For each current partition and corresponding set of active scenarios...
        logUtils(out, iter, "Building new partitions")
        nump = 0
        for partition in partitions
            nump +=1
            # For each current partition, get a set of active scenarios.
            listConstraint,listPoint = active_naive_graph(solved, problem, partition,nump)

            # Add the new constraints to the old ones.
            for con in listConstraint.constraints
                new_partition = deepcopy(partition)
                add_constraint!(new_partition, con)
                push!(new_partitions, new_partition)
            end
            for point in listPoint
                push!(finalListPoint,point)
            end
        end

        # Update the partitions for next iteration.
        partitions = new_partitions
        logUtils(out,iter,"number of partitions built " * string(length(partitions)) * " and number of points to compute lower bound " * string(length(finalListPoint)))
        # Solve deterministic model to obtain lower bound.
        logUtils(out, iter, "Solving lower bound model")
        lwsolved = solve_lower_bound_model_graph(problem, finalListPoint)
        lowerB = getobjectivevalue(lwsolved)
        logUtils(out, iter, "lower bound obtained " * string(lowerB))

        # Check bound gap relative to tolerance, as the exit criterion.
        gap = (upperB - lowerB) / abs(lowerB)
        logUtils(out, iter, "Relative bound gap at iteration $iter is $gap")
        if gap < tol
            logUtils(out, iter, "Bound gap below tolerance at iteration $iter")
            break
        end

    end

    logUtils(out, 0, "Procedure finished")
end

mainGraph()
