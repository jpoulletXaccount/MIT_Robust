using CSV
using Dates

include("../../src/read.jl")
include("../../src/constraints.jl")
include("../../src/model.jl")
include("../../src/deterministic.jl")
include("../../src/active.jl")
include("../../src/partitions.jl")
include("../../src/problem.jl")
include("../../src/utils.jl")
include("../../src/uncertainty.jl")


"""
    Objective

On the smallest problem size, show how the number of partitions, variables and constraints
grow on each iteration when using the Voronoi splitting approach.
"""

function main()

    # Get CL arg.
    possible = [50, 100, 150, 200, 300, 400, 500]
    num_jobs = possible[parse(Int, ARGS[1])]

    # Logging.
    out = stdout

    # Define run parameters.
<<<<<<< HEAD
    #num_workers = 22
=======
>>>>>>> firstResultsJulie
    cost_backup = 1
    M = 1e5
    intervals = [0, 400]

    # Initialise termination criterion variables.
<<<<<<< HEAD
    max_iters = 1
=======
    max_iters = 100
>>>>>>> e7e5eb9683cecb558052fbb4181a18066ec06eb1
    upper = Inf
    lower = -Inf
    tol = 0.05

    # Load jobs and create problem data.
    write_log(out, 0, "Loading problem.")
<<<<<<< HEAD
<<<<<<< HEAD
    jobs = build_jobs(joinpath(@__DIR__, "../../data/test/test-2.csv"), intervals)
    problemInterm = Problem(jobs, -1, 10, 1e5)
    num_worker = solve_deterministic_model(problemInterm)
    write_log(out, 0, "Nb workers need to solve the deterministic problem " * string(num_worker))
    problem = Problem(jobs, num_worker, cost_backup, M)
=======
    jobs = build_jobs(joinpath(@__DIR__, "../../data/actual/spread_50.csv"), intervals)
=======
    jobs = build_jobs(joinpath(@__DIR__, "../../data/actual/spread_$num_jobs.csv"), intervals)

>>>>>>> e7e5eb9683cecb558052fbb4181a18066ec06eb1
    # Get deterministic number of workers.
    first_problem = Problem(jobs, 0, 0, 0)
    num_workers = solve_deterministic_model(first_problem)
    write_log(out, 0, "# workers needed to solve the deterministic problem is $num_workers.")

<<<<<<< HEAD
    problem = Problem(jobs, num_workers, cost_backup, M)
>>>>>>> firstResultsJulie
=======
    # Create actual problem.
    problem = Problem(jobs, num_workers, cost_backup, M)
>>>>>>> e7e5eb9683cecb558052fbb4181a18066ec06eb1

    # Standard constraints on u.
    write_log(out, 0, "Defining base uncertainty.")
    base = standard_uncertainty(problem)

    # No other constraints in the single partition to start with.
    partitions = [base]

    # Start the partition and bound procedure.
    write_log(out, 0, "Starting main loop.")
    for i in 1:max_iters

        # Get current time to record iteration duration.
        start = now()

        # Solve the robust model with partitioned uncertainty.
        write_log(out, i, "Solving partition model with $(length(partitions)) partition(s).")
        partition_model = build_partition_model(problem, partitions)
        partition_status = solve(partition_model)

        # Get upper bound from solution.
        upper = getobjectivevalue(partition_model)
        write_log(out, i, "Upper bound obtained: $upper")

        # For each current partition, get a set of active scenarios.
        write_log(out, i, "Computing active scenarios.")
        points = greedy(partition_model, problem, partitions, true)

        # Solve deterministic model to obtain lower bound.
<<<<<<< HEAD
        num_points = sum(length(points[p]) for p in 1:length(points))
        write_log(out, i, "Solving lower bound model with $(num_points) points.")
        points_model = build_lower_bound_model(problem, points)
=======
        write_log(out, i, "Solving lower bound model with $(size(points[1],1)) points.")
        points_model = build_lower_bound_model(problem, points[1])
>>>>>>> firstResultsJulie
        points_status = solve(points_model)

        # Get lower bound from solution.
        lower = getobjectivevalue(points_model)
        write_log(out, i, "Lower bound obtained: $lower")

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
        gap = (upper - lower) / abs(upper)
        write_log(out, i, "Relative bound gap is $gap.")
<<<<<<< HEAD
        # Collect results.
        push!(results, [i,
                        length(partitions),
                        partition_model.numCols,
                        length(partition_model.linconstr),
                        getobjectivevalue(partition_model),
                        partition_status,
                        num_points,
                        points_model.numCols,
                        length(points_model.linconstr),
                        getobjectivevalue(points_model),
                        points_status,
                        now() - start])
=======

        # Create DataFrame to store results.
        results = DataFrame(iteration = i,
                            num_partitions = length(partitions),
                            upper_variables = partition_model.numCols,
                            upper_constraints = length(partition_model.linconstr),
                            upper_bound = getobjectivevalue(partition_model),
                            upper_status = partition_status,
                            num_points = num_points,
                            lower_variables = points_model.numCols,
                            lower_constraints = length(points_model.linconstr),
                            lower_bound = getobjectivevalue(points_model),
                            lower_status = points_status,
                            duration = now() - start)

        # Record results in the script's directory.
        write_log(out, 0, "Writing results to file.")
        if i == 1
            CSV.write(joinpath(@__DIR__, "voronoi_$num_jobs.csv"), results)
        else
            CSV.write(joinpath(@__DIR__, "voronoi_$num_jobs.csv"), results, append=true)
        end
>>>>>>> e7e5eb9683cecb558052fbb4181a18066ec06eb1

        if gap < tol
            write_log(out, i, "Bound gap below tolerance. Exiting.")
            break
<<<<<<< HEAD
        elseif (partition_status == :TIME_LIMIT) | (points_status == :TIME_LIMIT)
            write_log(out, i, "Either problem exceeded time limit. Exiting.")
            break
        end

=======
        end
>>>>>>> e7e5eb9683cecb558052fbb4181a18066ec06eb1
    end

    write_log(out, 0, "Procedure finished.")
end

main()
