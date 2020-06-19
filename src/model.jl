using JuMP, JuMPeR, Gurobi
using LinearAlgebra

include("constraints.jl")
include("problem.jl")


"""
    build_partition_model(problem::Problem, partitions::Array{ConstraintSet})::Model

Builds a model from problem data, with partitions on uncertainty.

partitions is an array of ConstraintSet objects for the constraints to be applied to each
partition separately (these constraints should be mutually exclusive, so that each partition
is indeed a partition).
"""
function build_partition_model(problem::Problem, partitions::Array{ConstraintSet})::Model

    # Get problem data.
    num_partitions = length(partitions)
    jobs = problem.jobs
    num_jobs = problem.num_jobs
    num_workers = problem.num_workers
    num_intervals = problem.num_intervals
    num_flights = problem.num_flights
    cost_backup = problem.cost_backup
    M = problem.M

    # Create model.
    m = RobustModel(solver=GurobiSolver(problem.gurobi_env))

    ### Static

    # Epigraph variable.
    @variable(m, r)

    # Assignment variables for workers to jobs.
    @variable(m, x[1:num_workers, 1:num_jobs], Bin)

    # Backup agent variables.
    @variable(m, z[1:num_jobs], Bin)

    ## Adjustable

    # Workers to jobs.
    @variable(m, x_u[1:num_partitions, 1:num_workers, 2:num_intervals, 1:num_jobs], Bin)

    # Backup agents.
    @variable(m, z_u[1:num_partitions, 2:num_intervals, 1:num_jobs], Bin)

    # Variables to capture clashes.
    @variable(m, y_u[1:num_partitions, 1:num_jobs, 1:num_jobs], Bin)

    ## Uncertainty

    # One uncertain vector for each partition, and one uncertain parameter for each flight.
    @uncertain(m, u[1:num_partitions, 1:num_flights])
    for p in 1:num_partitions
        apply_constraints!(m, u[p, :], partitions[p])
    end

    # Adjustable jobs cannot be performed in the first stage.
    adjustable_jobs = jobs_not_in_interval(jobs, 1)
    @constraint(m, z[adjustable_jobs] .== 0)
    for i in 1:num_workers

        # Again, we can't perform adjustable jobs in the first stage.
        @constraint(m, x[i, adjustable_jobs] .== 0)

        # We can only perform a job in the interval it nominally falls into.
        for t in 2:num_intervals, p in 1:num_partitions
            non_interval_jobs = jobs_not_in_interval(jobs, t)
            @constraint(m, x_u[p, i, t, :][non_interval_jobs] .== 0)
            @constraint(m, z_u[p, t, :][non_interval_jobs] .== 0)
        end
    end

    # Objective function.
    @objective(m, Min, r + cost_backup * sum(z))
    for p in 1:num_partitions
        @constraint(m,  r >= cost_backup * sum(sum(z_u[p, t, :] for t in 2:num_intervals)))
    end

    # Coverage constraints for first stage.
    for j in jobs_in_interval(jobs, 1)
        @constraint(m, sum(x[i, j] for i in 1:num_workers) + z[j] == 1)
    end

    # Coverage constraints for adjustable stages.
    for t in 2:num_intervals, j in jobs_in_interval(jobs, t), p in 1:num_partitions
        @constraint(m, sum(x_u[p, i, t, j] for i in 1:num_workers) + z_u[p, t, j] == 1)
    end

    # Chaining constraints to capture when jobs clash.
    for j in 1:num_jobs, k in 1:num_jobs, p in 1:num_partitions

        # We aren't interested in constraining for the same job.
        if j != k

            # Get flight indices (for uncertainty).
            j_index = jobs[j, :FLIGHT_INDEX]
            k_index = jobs[k, :FLIGHT_INDEX]

            # Get start and end times.
            j_end = jobs[j, :END]
            k_start = jobs[k, :START]

            @constraint(m, - M * y_u[p, j, k] <= k_start + u[p, k_index] - (j_end + u[p, j_index]))
            # @constraint(m, k_start + u[k] - (j_end + u[j]) <= M * (1 - y_u[j, k]))
        end
    end

    # No clashes for initial jobs with any other jobs.
    for i in 1:num_workers
        for j in 1:(num_jobs - 1), k in (j + 1):num_jobs, p in 1:num_partitions

            # Clashes of initial jobs with other initial jobs.
            @constraint(m, x[i, j] + x[i, k] <= 3 - (y_u[p, j, k] + y_u[p, k, j]))

            # Clashes of initial jobs with adjustable jobs.
            for t in 2:num_intervals
                @constraint(m, x[i, j] + x_u[p, i, t, k] <= 3 - (y_u[p, j, k] + y_u[p, k, j]))
            end
        end
    end

    # No clashes for adjustable jobs.
    for i in 1:num_workers
        for t in 2:num_intervals, t_bar in t:num_intervals, p in 1:num_partitions
            for j in 1:(num_jobs - 1), k in (j + 1):num_jobs
                @constraint(m, x_u[p, i, t, j] + x_u[p, i, t_bar, k] <= 3 - (y_u[p, j, k] + y_u[p, k, j]))
            end
        end
    end

    return m
end


"""
    build_lower_bound_model(problem::Problem, points::Array{<:Number})::Model

Builds a model from problem data, with variables for each discrete scenario in a list.

points is a set of uncertain scenarios obtained through heuristics. This model is fully
deterministic and is used to compute a lower bound.
"""
function build_lower_bound_model(problem::Problem, points)::Model

    # Get problem data.
<<<<<<< HEAD
    num_points = num_points = sum(length(points[p]) for p in 1:length(points))
    final_points = []
    for p in 1:length(points)
        for n in 1:length(points[p])
            push!(final_points,points[p][n])
        end
    end
=======
    num_points = size(points,1)
    @show(num_points)
>>>>>>> firstResultsJulie
    jobs = problem.jobs
    num_jobs = problem.num_jobs
    num_workers = problem.num_workers
    num_intervals = problem.num_intervals
    num_flights = problem.num_flights
    cost_backup = problem.cost_backup
    M = problem.M

    # Create model.
    m = Model(solver=GurobiSolver(problem.gurobi_env))

    ### Static

    # Epigraph variable.
    @variable(m, r)

    # Assignment variables for workers to jobs.
    @variable(m, x[1:num_workers, 1:num_jobs], Bin)

    # Backup agent variables.
    @variable(m, z[1:num_jobs], Bin)

    ## Adjustable

    # Workers to jobs.
    @variable(m, x_u[1:num_points, 1:num_workers, 2:num_intervals, 1:num_jobs], Bin)

    # Backup agents.
    @variable(m, z_u[1:num_points, 2:num_intervals, 1:num_jobs], Bin)

    # Variables to capture clashes.
    @variable(m, y_u[1:num_points, 1:num_jobs, 1:num_jobs], Bin)

    # Constraint on the tasks not assigned during interval t
    adjustable_jobs = jobs_not_in_interval(jobs, 1)
    @constraint(m, z[adjustable_jobs] .== 0)
    for i in 1:num_workers

        # Jobs in interval 1 are static.
        @constraint(m, x[i, adjustable_jobs] .== 0)

        # Remaining jobs are adjustable.
        for t in 2:num_intervals, p in 1:num_points
            non_interval_jobs = jobs_not_in_interval(jobs, t)
            @constraint(m, x_u[p, i, t, :][non_interval_jobs] .== 0)
            @constraint(m, z_u[p, t, :][non_interval_jobs] .== 0)
        end
    end

    # Objective function.
    @objective(m, Min, r + cost_backup * sum(z))
    for p in 1:num_points
        @constraint(m,  r >= cost_backup * sum(sum(z_u[p, t, :]) for t in 2:num_intervals))
    end

    # Coverage constraints for first stage.
    for j in jobs_in_interval(jobs, 1)
        @constraint(m, sum(x[i, j] for i in 1:num_workers) + z[j] == 1)
    end

    # Coverage constraints for adjustable stages.
    for t in 2:num_intervals, j in jobs_in_interval(jobs, t), p in 1:num_points
        @constraint(m, sum(x_u[p, i, t, j] for i in 1:num_workers) + z_u[p, t, j] == 1)
    end

    println("chaining")
    # Chaining constraints for jobs.
    for j in 1:num_jobs, k in 1:num_jobs, p in 1:num_points

        # We aren't interested in a clash between a job and itself.
        if j != k

            # Get flight indices (for uncertainty).
            j_index = jobs[j, :FLIGHT_INDEX]
            k_index = jobs[k, :FLIGHT_INDEX]

            # Get start and end times.
            j_end = jobs[j, :END]
            k_start = jobs[k, :START]

            @constraint(m, - M * y_u[p, j, k] <= k_start + final_points[p][k_index] - (j_end + final_points[p][j_index]))
            # @constraint(m, k_start + u[k] - (j_end + u[j]) <= M * (1 - y_u[j, k]))
        end
    end

    println("chaining2")
    # No clash for initial jobs with any other jobs.
    for i in 1:num_workers
        for j in 1:(num_jobs - 1), k in (j + 1):num_jobs, p in 1:num_points

            # Clashes of initial jobs with other initial jobs.
            @constraint(m, x[i, j] + x[i, k] <= 3 - (y_u[p, j, k] + y_u[p, k, j]))

            # Clashes of initial jobs with adjustable jobs.
            for t in 2:num_intervals
                @constraint(m, x[i, j] + x_u[p, i, t, k] <= 3 - (y_u[p, j, k] + y_u[p, k, j]))
            end
        end
    end
    println("chaining3")
    # No clash for adjustable jobs.
    for i in 1:num_workers
        for t in 2:num_intervals, t_bar in t:num_intervals, p in 1:num_points
            for j in 1:(num_jobs - 1), k in (j + 1):num_jobs
                @constraint(m, x_u[p, i, t, j] + x_u[p, i, t_bar, k] <= 3 - (y_u[p, j, k] + y_u[p, k, j]))
            end
        end
    end
    println("returning")

    return m
end
