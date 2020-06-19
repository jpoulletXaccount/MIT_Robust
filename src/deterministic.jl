using JuMP, Gurobi

include("problem.jl")

"""
function chaining_matrix(problem::Problem)::Array{Int}

Return a matrix of all chainging. If the value is zero then job k cannot
be done after job j. Otherwise value is 1

"""

function chaining_matrix(problem::Problem)::Array{Int}

    # Get relevant data.
    num_jobs = problem.num_jobs
    jobs = problem.jobs

    # Create model.
    m = Model(solver=GurobiSolver(problem.gurobi_env, OutputFlag=0))

    # Initial matrix is all zeros (no chaining).
    link_matrix = fill(0, (num_jobs, num_jobs))

    # Only look at links between two first stage jobs. A link_matrix between a first and a second
    # stage job can be a MAYBE, which updates depending on the partition.
    for j in 1:num_jobs, k in 1:num_jobs

        # Get nominal start and end times.
        j_end = jobs[j, :END]
        k_start = jobs[k, :START]

        # Check whether we *ever* have  k >= j.
        if k_start - j_end >= 0
            link_matrix[j, k] = 1
        end
    end
    return link_matrix
end

"""
solve_deterministic_model(problem::Problem)::Int

determine how many workers are needed to complete the tasks in the determinic cases
Return the number of workers
"""

function solve_deterministic_model(problem::Problem)::Int

    # Get problem data.
    cost_backup = problem.cost_backup
    num_jobs = problem.num_jobs
    jobs = problem.jobs
    link_matrix = chaining_matrix(problem)

    # Get indicies of dummy nodes.
    source = 0
    sink = num_jobs + 1

    # Create the model (corresponding to multiple graphs, one per partition).
    m = Model(solver=GurobiSolver(problem.gurobi_env, OutputFlag=0))

    # Worker flow variables (one for each arc).
    @variable(m, x[source:sink, source:sink], Bin)
    @variable(m,totalWorker)
    # Objective function. Minimize the number of workers
    @objective(m, Min, totalWorker)

    # No flow into source or out of sink.
    @constraint(m, [j = source:sink], x[j, source] == 0)
    @constraint(m, [j = source:sink], x[sink, j] == 0)

    # Set all arcs which can't be included in this partition to zero.
    discard = findall(map(x -> x in [0], link_matrix))
    for (j, k) in map(x -> x.I, discard)
        @constraint(m, x[j, k] == 0)
    end

    # Ensure all jobs are completed by either a worker or backup agent.
    for k in 1:num_jobs
        @constraint(m, sum(x[j, k] for j in source:sink) == 1)
    end

    # Flow constraints for workers at each node.
    for j in 1:num_jobs
        @constraint(m,
            sum(x[j, k] for k in source:sink) ==
            sum(x[i, j] for i in source:sink))
    end

    #Flow constraints at sink and source
    @constraint(m, sum(x[source, j] for j in 1:num_jobs) == totalWorker)
    @constraint(m, sum(x[j,sink] for j in 1:num_jobs) == totalWorker)

    solve(m)

    return getvalue(totalWorker)
end
