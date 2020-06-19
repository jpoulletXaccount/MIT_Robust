using JuMP
using Random

include("problem.jl")
include("constraints.jl")


"""
    Active Scenario Functions

All functions that retrieve active scenarios should have the following signature:

    active(solved::Model, problem::Problem, partitions::Array{ConstraintSet})

It should accept a solved model, problem parameters, and uncertainty partition constraints -- and
return an array of matrices which correspond to active scenarios in each partition.
"""


"""
    basic(solved::Model, problem::Problem, partitions::Array{ConstraintSet})

Use a silly basic strategy to find realisations of uncertainty where jobs overlap. Doesn't
take into account any information from the solved model.
"""
function basic(solved::Model, problem::Problem, partitions::Array{ConstraintSet})

    # Get problem data.
    jobs = problem.jobs
    num_jobs = problem.num_jobs
    num_workers = problem.num_workers
    num_flights = problem.num_flights

    # Create array to store points.
    points = []

    # Get a set of active points for each partition.
    for partition in partitions

        partition_points = []

        for j in 1:num_jobs, k in 1:num_jobs

            m = Model(solver=GurobiSolver(problem.gurobi_env, OutputFlag=0))

            # Uncertainty is now a variable, which we constrain to lie within the
            # relevant partition.
            @variable(m, u[1:num_flights])
            apply_constraints!(m, u, partition)

            # Get flight indicies (for using to reference uncertain vector).
            j_index = jobs[j, :FLIGHT_INDEX]
            k_index = jobs[k, :FLIGHT_INDEX]

            # Get start and end times for pair of jobs.
            j_end = jobs[j, :END]
            k_start = jobs[k, :START]

            # Look for scenarios within the current partition in which jobs overlap.
            @objective(m, Min, k_start + u[k_index] - (j_end + u[j_index]))

            # Solve the model, and save the scenario if the jobs overlap.
            solve(m)
            if getobjectivevalue(m) < 0
                push!(partition_points, getvalue(u))
            end
        end

        # Add a matrix of the points to the list, or an empty one if none were found.
        if length(partition_points) > 0
            partition_matrix = transpose(hcat(partition_points...))
            push!(points, unique(partition_matrix, dims=1))
        else
            push!(points, [])
        end
    end

    return points
end


"""
    follow_paper(solved::Model, problem::Problem, partitions::Array{ConstraintSet})

Follow exactly the reference paper of Bertsimas and Dunning, but does not make a lot of sense in our case since it
is equivalent to just divide any single option of delay.
"""
function follow_paper(solved::Model, problem::Problem, partitions::Array{ConstraintSet})

    # Get problem data.
    jobs = problem.jobs
    num_jobs = problem.num_jobs
    num_workers = problem.num_workers
    num_flights = problem.num_flights
    M = problem.M
    y_u = getindex(solved, :y_u)
    y_u = getvalue(y_u)

    # Create array to store points.
    points = []

    # Get a set of active points for each partition.
    comp = 0
    for partition in partitions
        comp +=1
        partition_points = []

        for j in 1:num_jobs, k in 1:num_jobs

            m = Model(solver=GurobiSolver(problem.gurobi_env,OutputFlag=0))

            # Uncertainty is now a variable, which we constrain to lie within the
            # relevant partition.
            @variable(m, u[1:num_flights])
            apply_constraints!(m, u, partition)

            # Look for scenarios within the current partition in which jobs overlap.
            j_index = jobs[j, :FLIGHT_INDEX]
            k_index = jobs[k, :FLIGHT_INDEX]
            j_end = jobs[j, :END]
            k_start = jobs[k, :START]
            @objective(m, Min, -M * y_u[comp, j, k] - (k_start + u[k_index]) + (j_end + u[j_index]))

            # Solve the model, and save the scenario if the jobs overlap.
            solve(m)
            push!(partition_points, getvalue(u))
        end

        # Add a matrix of the points to the list.
        if length(partition_points) > 0
            partition_matrix = transpose(hcat(partition_points...))
            push!(points, unique(partition_matrix, dims=1))
        end
    end

    return points

end


"""
    greedy(solved::Model, problem::Problem, partitions::Array{ConstraintSet})

Look at all the jobs currently performed by a backup agent, and compare these jobs to ones which
don't clash in the nominal case. These uncertain realisations may lie on the 'border' between the
optimal decision partitions.
"""
function greedy(solved::Model, problem::Problem, partitions::Array{ConstraintSet}, dichotomic::Bool)

    # Get problem data.
    jobs = problem.jobs
    num_jobs = problem.num_jobs
    num_workers = problem.num_workers
    num_intervals = problem.num_intervals
    num_flights = problem.num_flights
    M = problem.M

    # Create array to store points.
    points = []

    # Get results obtained with the MIP.
    y_u = getvalue(getindex(solved, :y_u))
    z = getvalue(getindex(solved, :z))
    z_u = getvalue(getindex(solved, :z_u))

    # Get a set of active points for each partition.
    for (p, partition) in enumerate(partitions)

        # Get all jobs performed by backup agents.
        jobs_performed_backup = findall(z .> 0.5)
        for t in 2:num_intervals
            jobs_performed_backup = vcat(jobs_performed_backup, findall(z_u[p, t, :] .> 0.5))
        end

        partition_points = []

        for j in jobs_performed_backup, k in 1:num_jobs

            # Get flight indicies (for indexing uncertainty).
            j_index = jobs[j, :FLIGHT_INDEX]
            k_index = jobs[k, :FLIGHT_INDEX]

            # Get start and end times of jobs.
            j_end = jobs[j, :END]
            k_start = jobs[k, :START]

            # Only consider the case where k starts after j finishes in the nominal case
            # (i.e. no nominal clash).
            if k_start >= j_end

                m = Model(solver=GurobiSolver(problem.gurobi_env, OutputFlag=0))

                # Uncertainty is now a variable, which we constrain to lie within the
                # relevant partition.
                @variable(m, u[1:num_flights])
                apply_constraints!(m, u, partition)

                # Look for scenarios within the current partition in which jobs overlap.
                @objective(m, Min, M * y_u[p, j, k] + (k_start + u[k_index]) - (j_end + u[j_index]))

                # Solve the model, and save the scenario.
                solve(m)
                push!(partition_points, getvalue(u))
            end
        end

        # Add the points to the list.
        if length(partition_points) > 0
            if dichotomic
                unique!(partition_points)
                shuffle!(partition_points)
                push!(points, [partition_points[1],partition_points[2]])
            else
                push!(points, unique(partition_points))
            end
        else
            push!(points, [])
        end
    end

    return points
end


"""
    active_naive_graph(solved::Model, problem::Problem, partition::ConstraintSet, partition_index::Int)

Try to take a look at every job currently done by agent, and compare the link to a task which doesn't
clash in the nominal case. Obtain two partitions corresponding to halfspaces.
"""
function active_naive_graph(solved::Model, problem::Problem, partition::ConstraintSet, partition_index::Int)

    # Get problem data.
    jobs = problem.jobs
    num_jobs = problem.num_jobs
    num_workers = problem.num_workers
    num_intervals = problem.num_intervals
    num_flights = problem.num_flights
    M = problem.M

    # Create list to hold points.
    points = []

    # Get variable values from solved MIP.
    x = getvalue(getindex(solved, :x))
    y = getvalue(getindex(solved, :y))

    # Get the set of jobs performed by backup.
    backup_jobs = findall(y[partition_index, :] .> 0.5)
    agent_jobs = setdiff(1:num_jobs, backup_jobs)

    for j in backup_jobs, k in agent_jobs

        j_index = jobs[j, :FLIGHT_INDEX]
        k_index = jobs[k, :FLIGHT_INDEX]
        j_end = jobs[j, :END]
        k_start = jobs[k, :START]

        # Only look at this pair of jobs if we have j -> k in the nominal case.
        if k_start >= j_end

            m = Model(solver=GurobiSolver(problem.gurobi_env, OutputFlag=0))

            # Uncertainty is now a variable, which we constrain to lie within the
            # relevant partition.
            @variable(m, u[1:num_flights])
            apply_constraints!(m, u, partition)

            # Look for scenarios within the current partition in which jobs overlap.
            @objective(m, Min, (k_start + u[k_index]) - (j_end + u[j_index]))

            # Solve the model, and save the scenario if the jobs overlap.
            solve(m)

            # If there is a realisation of uncertainty for which these two jobs overlap.
            if getobjectivevalue(m) < 0

                # Record the corresponding value of u, and create the two partitions.
                points = [getvalue(u)]
                left = ConstraintSet(GeneralConstraint(z -> k_start + z[k] - (j_end + z[j])))
                right = ConstraintSet(GeneralConstraint(z -> -(k_start + z[k]) + (j_end + z[j])))

                # Get another point for the lower bound model while we're at it.
                @objective(m, Max, (k_start + u[k]) - (j_end + u[j]))
                solve(m)
                push!(points, getvalue(u))

                return [left, right], points
            end
        end
    end

    # If we didn't find any new links to add, then return an empty ConstraintSet.
    return [ConstraintSet()], []
end
