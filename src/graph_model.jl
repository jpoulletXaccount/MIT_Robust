using JuMP, Gurobi

include("problem.jl")

# Create constants for entries in the link matrix.
ALWAYS =  1
MAYBE  =  0
NEVER  = -1


"""
    Link Matrices

A link matrix is a (num_jobs, num_jobs) matrix where entries indicate the link status between
two jobs under the current partition. There are three statuses:

    ALWAYS: a worker can ALWAYS complete j -> k.
    MAYBE: depending on the realisation of uncertainty, a worker can complete j -> k.
    NEVER: for all uncertainty in the current partition, a worker cannot complete j -> k.

Statuses between two first stage jobs must be either ALWAYS or NEVER -- these are not allowed to
update based on a partition.
"""


"""
    initial_link_matrix(problem::Problem, base::ConstraintSet)::Array{Int}

Create an initial link matrix for the first stage jobs, which have to have their links
set to either 'ALWAYS' or 'NEVER' from the base uncertainty. This is slightly different
behaviour when compared with second stage jobs.

Then update the link matrix for second stage jobs.
"""
function initial_link_matrix(problem::Problem, base::ConstraintSet)::Array{Int}

    # Get relevant data.
    num_jobs = problem.num_jobs
    num_flights = problem.num_flights
    jobs = problem.jobs

    # Create model.
    m = Model(solver=GurobiSolver(problem.gurobi_env, OutputFlag=0))

    # Uncertainty is constrained to lie within the base partition.
    @variable(m, u[1:num_flights])
    apply_constraints!(m, u, base)

    # Initial matrix is all MAYBEs.
    link_matrix = fill(MAYBE, (num_jobs, num_jobs))

    # Only look at links between two first stage jobs. A link_matrix between a first and a second
    # stage job can be a MAYBE, which updates depending on the partition.
    println("Nb jobs in first interval " * string(length(jobs_in_interval(jobs,1))))
    for j in jobs_in_interval(jobs, 1), k in jobs_in_interval(jobs, 1)

        # Get flight indices (for uncertainty).
        j_index = jobs[j, :FLIGHT_INDEX]
        k_index = jobs[k, :FLIGHT_INDEX]

        # Get start and end times.
        j_end = jobs[j, :END]
        k_start = jobs[k, :START]

        # Check whether we *ever* have j -\> k. We need to classify each first stage link_matrix as
        # either NEVER or ALWAYS.
        @objective(m, Min, (k_start + u[k_index]) - (j_end + u[j_index]))
        solve(m)
        if getobjectivevalue(m) < 0
            link_matrix[j, k] = NEVER
        else
            link_matrix[j, k] = ALWAYS
        end
    end

    # Update links for second stage.
    update_link_matrix_set!(problem, base, link_matrix)

    return link_matrix
end


"""
    update_link_matrix_set!(problem::Problem, partition::ConstraintSet, link_matrix::Array{Int})

Update a link matrix based on a particular partition of uncertainty.

Only the 'MAYBE' elements in the link matrix are updated (for speed), so we need to be careful that
link matrices are only updated based on successive partitions for which the second is a subset of the
first.
"""
function update_link_matrix_set!(problem::Problem, partition::ConstraintSet, link_matrix::Array{Int})

    # Get relevant data.
    num_jobs = problem.num_jobs
    num_flights = problem.num_flights
    jobs = problem.jobs

    # Create model.
    m = Model(solver=GurobiSolver(problem.gurobi_env, OutputFlag=0))

    # Uncertainty is constrain to lie within the partition.
    @variable(m, u[1:num_flights])
    apply_constraints!(m, u, partition)

    # Only update pairs of jobs where we are unsure of a clash.
    check = findall(link_matrix .== MAYBE)
    nbLinkInitial = length(check)
    # Loop over by retrieving indices from the CartesianIndex objects returned previously.
    for (j, k) in map(x -> x.I, check)

        # Get flight indices (for uncertainty).
        j_index = jobs[j, :FLIGHT_INDEX]
        k_index = jobs[k, :FLIGHT_INDEX]

        # Get start and end times.
        j_end = jobs[j, :END]
        k_start = jobs[k, :START]

        # Check whether the jobs can ALWAYS be chained.
        @objective(m, Min, (k_start + u[k_index]) - (j_end + u[j_index]))
        solve(m)
        if getobjectivevalue(m) >= 0
            link_matrix[j, k] = ALWAYS
        end

        # Check whether the jobs can NEVER be chained.
        @objective(m, Max, (k_start + u[k_index]) - (j_end + u[j_index]))
        solve(m)
        if getobjectivevalue(m) < 0
            link_matrix[j, k] = NEVER
        end
    end
    check = findall(link_matrix .== MAYBE)
    nbLinkFinal = length(check)
    if nbLinkFinal >= nbLinkInitial
        println("Warning we didn't decrease the number of maybe links")
    end
end


"""
    update_link_matrix_point!(problem::Problem, u::Array{Float64}, link_matrix::Array{Int})

Update a link matrix for a particular realisation of uncertainty.

This is used to compute the link matrices for realisations of uncertainty used in obtaining lower
bounds.
"""
function update_link_matrix_point!(problem::Problem, u::Array{Float64}, link_matrix::Array{Int})

    # Get relevant data.
    num_jobs = problem.num_jobs
    num_flights = problem.num_flights
    jobs = problem.jobs

    # Only update pairs of jobs where we are unsure of a clash.
    clashes = findall(link_matrix .== MAYBE)

    # Loop over by retrieving indices from the CartesianIndex objects returned previously.
    for (j, k) in map(x -> x.I, clashes)

        # Get flight indices (for uncertainty).
        j_index = jobs[j, :FLIGHT_INDEX]
        k_index = jobs[k, :FLIGHT_INDEX]

        # Get start and end times.
        j_end = jobs[j, :END]
        k_start = jobs[k, :START]

        # Check whether the jobs can ALWAYS be chained.
        if (k_start + u[k_index]) - (j_end + u[j_index]) >= 0
            link_matrix[j, k] = ALWAYS
        else
            # jobs can never be chained
            link_matrix[j, k] = NEVER
        end
    end

end


"""
    get_outgoing(link_matrix, j)

Return a list of job IDs, [k], where the type of arc from j -> k is ALWAYS.
"""
function get_outgoing(link_matrix, j, sink)
    return vcat(findall(link_matrix[j, :] .== ALWAYS), [sink])
end


"""
    get_incoming(link_matrix, k)

Return a list of job IDs, [j], where the type of arc from j -> k is ALWAYS.
"""
function get_incoming(link_matrix, k, source)
    return vcat([source], findall(link_matrix[:, k] .== ALWAYS))
end


"""
    solve_graph_model(problem::Problem, link_matrices)::Model

Builds and solves a graph model from problem data and link matrices (one for each partition or
discrete uncertain point).
"""
function solve_graph_model(problem::Problem, link_matrices, zero_statuses, relaxed)::Model

    # Get problem data.
    cost_backup = problem.cost_backup
    num_workers = problem.num_workers
    num_jobs = problem.num_jobs
    jobs = problem.jobs
    num_partitions = length(link_matrices)

    # Get indicies of dummy nodes.
    source = 0
    sink = num_jobs + 1

    # Create the model (corresponding to multiple graphs, one per partition).
    m = Model(solver=GurobiSolver(problem.gurobi_env))

    # Worker flow variables (one for each arc).
    if relaxed
        @variable(m, x[1:num_partitions, source:sink, source:sink] >= 0)
    else
        @variable(m, x[1:num_partitions, source:sink, source:sink], Bin)
    end

    # Backup agent variables (one for each node, excluding source and sink).
    if relaxed
        @variable(m, y[1:num_partitions, 1:num_jobs] >= 0)
    else
        @variable(m, y[1:num_partitions, 1:num_jobs], Bin)
    end

    # Epigraph variable for objective.
    @variable(m, z)

    # Objective function.
    @objective(m, Min, cost_backup * z)

    # No flow into source or out of sink.
    @constraint(m, [p = 1:num_partitions, j = source:sink], x[p, j, source] == 0)
    @constraint(m, [p = 1:num_partitions, j = source:sink], x[p, sink, j] == 0)

    # Ensure that all first stage decisions are consistent across partitions.
    for p in 2:num_partitions, j in source:sink, k in jobs_in_interval(jobs, 1)
            @constraint(m, x[p, j, k] == x[1, j, k])
    end

    # Loop through partitions and links to create second stage constraints.
    for (p, link_matrix) in enumerate(link_matrices)

        # Set all arcs which can't be included in this partition to zero.
        discard = findall(map(x -> x in zero_statuses, link_matrix))
        for (j, k) in map(x -> x.I, discard)
            @constraint(m, x[p, j, k] == 0)
        end

        # Add epigraph constraint to objective.
        @constraint(m, z >= sum(y[p, j] for j in 1:num_jobs))

        # Add the worker constraints on source and sink dummy nodes.
        @constraint(m, sum(x[p, source, k] for k in 1:num_jobs) <= num_workers)
        @constraint(m, sum(x[p, j, sink] for j in 1:num_jobs) <= num_workers)  # Not necessary, but interpretable.

        # Ensure all jobs are completed by either a worker or backup agent.
        for k in 1:num_jobs
            @constraint(m, sum(x[p, j, k] for j in source:sink) + y[p, k] == 1)
        end

        # Flow constraints for workers at each node.
        for j in 1:num_jobs
            @constraint(m, sum(x[p, j, k] for k in source:sink) == sum(x[p, i, j] for i in source:sink))
        end

    end

    solve(m)

    return m
end


function build_new_partitions_link_matrices_naive(problem, solved, partitions, link_matrices,lb_points,lb_link_matrices)

    # Initialise lists to store all the new things.
    new_partitions = ConstraintSet[]
    new_link_matrices = []


    i = 1
    for (partition, link_matrix) in zip(partitions, link_matrices)

        # Create subpartitions and obtain points for lower bounding. Need to update this to be
        # more sensible.
        subpartitions, point_list = active_naive_graph(solved, problem, partition, i)

        # For each of the subpartitions...
        for subpartition in subpartitions

            # HACK: to deal with the case where the subpartitions contain no constraints. Need
            # to fix this.
            new_partition = deepcopy(partition)

            # Add the new constraints to the ones from the parent partition.
            for constraint in subpartition.constraints
                new_partition = deepcopy(partition)
                add_constraint!(new_partition, constraint)
            end
            push!(new_partitions, new_partition)

            # Now that we have a new partition, use this to update the corresponding link matrix.
            new_link_matrix = deepcopy(link_matrix)
            update_link_matrix_set!(problem, new_partition, new_link_matrix)
            push!(new_link_matrices, new_link_matrix)
        end

        # Now go through the points returned, and create link matrices for each of them.
        for point in point_list
            push!(lb_points, point)
            new_lb_link_matrix = deepcopy(link_matrix)
            update_link_matrix_point!(problem, point, new_lb_link_matrix)
            push!(lb_link_matrices, new_lb_link_matrix)
        end

    end
    if length(new_partitions) != length(new_link_matrices)
        println("Attention error on the new partitions or new_link_matrix")
        @show(length(new_partitions),length(new_link_matrices))
        return
    end
    return new_partitions, new_link_matrices, lb_points, lb_link_matrices
end


"""
function build_new_partitions_link_matrices_relax_MAYBE(problem, solved, partitions, link_matrices,lb_points,lb_link_matrices)

solve the graph with all links maybe and create subpartitions and points
"""

function build_new_partitions_link_matrices_relax_MAYBE(problem, solvedFewLink, partitions, link_matrices,lb_points,lb_link_matrices)

    # Initialise lists to store all the new things.
    new_partitions = ConstraintSet[]
    new_link_matrices = []
    cost_backup = problem.cost_backup
    num_jobs = problem.num_jobs

    # Relax the requirement that MAYBEs cannot be inclded as links.
    solved = solve_graph_model(problem, link_matrices, (NEVER),true)
    lb = getobjectivevalue(solved)

    # Examine the worst case partition (so that we don't partition too many times at
    # each iteration).

    for p in 1:length(partitions)
        partition = partitions[p]
        link_matrix = link_matrices[p]
        # Create subpartitions and obtain points for lower bounding. Need to update this to be
        # more sensible.
        subpartitions, criticalPoint = relax_maybes_Dual(solved, problem, link_matrix, p)

        # For each of the subpartitions...
        for subpartition in subpartitions

            # HACK: to deal with the case where the subpartitions contain no constraints. Need
            # to fix this.
            new_partition = deepcopy(partition)

            # Add the new constraints to the ones from the parent partition.
            for constraint in subpartition.constraints
                new_partition = deepcopy(partition)
                add_constraint!(new_partition, constraint)
            end
            if check_feasible_partition(new_partition,problem)
                push!(new_partitions, new_partition)

                # Now that we have a new partition, use this to update the corresponding link matrix.
                new_link_matrix = deepcopy(link_matrix)
                update_link_matrix_set!(problem, new_partition, new_link_matrix)
                push!(new_link_matrices, new_link_matrix)

                # Draw clever points from this partition
                listPoint = draw_points(new_partition,problem,criticalPoint[1],criticalPoint[2])

                # Now go through the points returned, and create link matrices for each of them.
                for point in listPoint
                    push!(lb_points, point)
                    new_lb_link_matrix = deepcopy(new_link_matrix)
                    update_link_matrix_point!(problem, point, new_lb_link_matrix)
                    push!(lb_link_matrices, new_lb_link_matrix)
                end
            end
        end

    end
    if length(new_partitions) != length(new_link_matrices)
        println("Attention error on the new partitions or new_link_matrix")
        @show(length(new_partitions),length(new_link_matrices))
        return
    end

    return new_partitions, new_link_matrices, lb_points, lb_link_matrices, lb
end


function build_new_partitions_link_matrices_relax_MAYBE_worstPartition(problem, solvedFewLink, partitions, link_matrices,lb_points,lb_link_matrices)

    # Initialise lists to store all the new things.
    new_partitions = ConstraintSet[]
    new_link_matrices = []
    partitions_to_be_added = ConstraintSet[]
    link_matrices_to_be_added = []
    cost_backup = problem.cost_backup
    num_jobs = problem.num_jobs

    # Relax the requirement that MAYBEs cannot be inclded as links.
    solved = solve_graph_model(problem, link_matrices, (NEVER),true)
    lb = getobjectivevalue(solved)

    # Examine the worst case partition (so that we don't partition too many times at
    # each iteration).
    worst_partition = []
    maxPartition = -1
    maxTest = -1
    for p in 1:length(partitions)
        solvedTest = solve_graph_model(problem, [link_matrices[p]], (NEVER,MAYBE),true)
        y = getindex(solvedTest, :y)

        test = sum(getvalue(y[1, j]) for j in 1:num_jobs)
        if cost_backup * test == getobjectivevalue(solvedFewLink)
            push!(worst_partition, p)
            partition = partitions[p]
            link_matrix = link_matrices[p]

            # Create subpartitions and obtain points for lower bounding. Need to update this to be
            # more sensible.
            find_link_divide(problem, partition, link_matrix,new_partitions,new_link_matrices,lb_points,lb_link_matrices,0)

        else
            push!(new_partitions,partitions[p])
            push!(new_link_matrices, link_matrices[p])
            if test > maxTest
                maxTest = test
                maxPartition = length(new_partitions)
            end
        end
    end

    println("Number of worst partition " * string(length(worst_partition)))
    if length(worst_partition) == 0
        """
        partition = new_partitions[maxPartition]
        link_matrix = new_link_matrices[maxPartition]
        deleteat!(new_partitions,maxPartition)
        deleteat!(new_link_matrices,maxPartition)
        # Create subpartitions and obtain points for lower bounding. Need to update this to be
        # more sensible.
        find_link_divide(problem, partition, link_matrix,new_partitions,new_link_matrices,lb_points,lb_link_matrices,0)
        """
        return build_new_partitions_link_matrices_relax_MAYBE(problem, solved, new_partitions, new_link_matrices,lb_points,lb_link_matrices)
    end

    if length(new_partitions) != length(new_link_matrices)
        println("Attention error on the new partitions or new_link_matrix")
        @show(length(new_partitions),length(new_link_matrices))
        return
    end

    return new_partitions, new_link_matrices, lb_points, lb_link_matrices, lb
end

"""
function find_link_divide(problem, partition, link_matrix,new_partitions,new_link_matrices,lb_points,lb_link_matrices)

    find the link to divide and update everything
"""

function find_link_divide(problem, partition, link_matrix,new_partitions,new_link_matrices,lb_points,lb_link_matrices,iter)
    # Create subpartitions and obtain points for lower bounding. Need to update this to be
    # more sensible.

    solvedTestMaybe = solve_graph_model(problem, [link_matrix], (NEVER),true)
    solvedTest = solve_graph_model(problem, [link_matrix], (NEVER,MAYBE),true)
    initialValue = getobjectivevalue(solvedTest)
    @show(initialValue)

    iter += 1

    pos_partition = -1

    subpartitions, criticalPoint = relax_maybes_Dual(solvedTestMaybe, problem, link_matrix, 1)

    # For each of the subpartitions...
    for subpartition in subpartitions

        # HACK: to deal with the case where the subpartitions contain no constraints. Need
        # to fix this.
        new_partition = deepcopy(partition)

        # Add the new constraints to the ones from the parent partition.
        for constraint in subpartition.constraints
            new_partition = deepcopy(partition)
            add_constraint!(new_partition, constraint)
        end
        if check_feasible_partition(new_partition,problem)
            push!(new_partitions, new_partition)

            # Now that we have a new partition, use this to update the corresponding link matrix.
            new_link_matrix = deepcopy(link_matrix)
            update_link_matrix_set!(problem, new_partition, new_link_matrix)
            push!(new_link_matrices, new_link_matrix)

            # Test if we have improve or not
            solvedTest = solve_graph_model(problem, [new_link_matrix], (NEVER,MAYBE),true)
            finalValue = getobjectivevalue(solvedTest)
            @show(finalValue)
            if finalValue == initialValue
                pos_partition = length(new_partitions)
            end

            # Draw clever points from this partition
            listPoint = draw_points(new_partition,problem,new_link_matrix)

            # Now go through the points returned, and create link matrices for each of them.
            for point in listPoint
                push!(lb_points, point)
                new_lb_link_matrix = deepcopy(new_link_matrix)
                update_link_matrix_point!(problem, point, new_lb_link_matrix)
                push!(lb_link_matrices, new_lb_link_matrix)
            end
        end
    end

    if iter <=5 && pos_partition != -1
        partition = new_partitions[pos_partition]
        link_matrix = new_link_matrices[pos_partition]
        deleteat!(new_partitions,pos_partition)
        deleteat!(new_link_matrices,pos_partition)
        find_link_divide(problem, partition, link_matrix,new_partitions,new_link_matrices,lb_points,lb_link_matrices,iter)
    end

end
