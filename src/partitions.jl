using LinearAlgebra
using JuMP

include("constraints.jl")

EPS = 0.1
"""
    Partition Functions

By design, a valid partition function should have the following signature:

    partition(points::Array{<:Number})::Array{ConstraintSet}

It should accept an array of points, and return an array of ConstraintSet objects that define
partitions of the space of uncertainties.

The input array should have size (m x n), where m is the number of points and n is the dimension
of the space.
"""


"""
    random(points::Array{<:Number})::Array{ConstraintSet}

Selects two random points and returns partitions defined by the maximum margin separating
hyperplane.
"""
function random(points::Array{<:Number})::Array{ConstraintSet}

    # Get number of points and dimension of space.
    m, n = size(points)

    # Get random indices.
    i, j = rand(1:m, 2)

    # Extract actual points.
    x = points[i, :]
    y = points[j, :]

    partitions = [
        ConstraintSet([GeneralConstraint(max_margin(x, y))]),
        ConstraintSet([GeneralConstraint(max_margin(y, x))])
    ]

    return partitions
end


"""
    voronoi(points::Array{<:Number})::Array{ConstraintSet}

Constructs Voronoi regions around individual points, returning the same number of partitions
as input points. Not efficient, and unlikely to be practical with many points.
"""
function voronoi(points)::Array{ConstraintSet}

    # Get number of points and dimension of space.
    m = size(points,1)
    # One partition for each point.
    partitions = [ConstraintSet() for i in 1:m]

    for i in 1:(m - 1), j in (i + 1):m

        # Get current pair of points.
        x = points[i]
        y = points[j]

        # Add opposing constraints for the relevant partitions.
        add_constraint!(partitions[i], GeneralConstraint(max_margin(x, y)))
        add_constraint!(partitions[j], GeneralConstraint(max_margin(y, x)))
    end

    return partitions
end


"""
    max_margin(x, y)

Returns a function f, such that f(z) <= 0 iff z is closer to x than y according to the
Euclidean norm.
"""
function max_margin(x, y)

    # Compute value of affine function at midpoint between x and y.
    alpha = dot(y - x, 0.5 * (x + y))

    # Return function in the required form.
    return z -> dot(y - x, z) - alpha
end


"""
function relax_maybes(solved, problem, link_matrix, p)

Relax every MAYBE link, and then pick some into the solution
Currently random, but we may use dual to choose which one
"""

function relax_maybes(solved, problem, link_matrix, p)
    num_jobs = problem.num_jobs
    jobs = problem.jobs
    x = getindex(solved, :x)

    potentialLink = findall(link_matrix[:, :] .== MAYBE)

    # Select arc (j, k) to remove. Maybe look at duals.
    # Todo currently a bit random
    maxLink = -1
    for (j, k) in map(x -> x.I, potentialLink)
        if getvalue(x[p,j,k])> 0.5
            maxLink = (j,k)
            break
        end
    end

    # check if a MAYBE link has been used
    if maxLink == -1
        j = rand(1:num_jobs)
        k = rand(1:num_jobs)
        return [ConstraintSet()],(j,k)
    else
        (j,k) = maxLink

        # Get start and end times for these two jobs.
        j_end = jobs[j, :END]
        k_start = jobs[k, :START]

        # Create left and right halfspaces.
        left_set = ConstraintSet(GeneralConstraint(z -> k_start + z[k] - (j_end + z[j])+ EPS))
        right_set = ConstraintSet(GeneralConstraint(z -> -(k_start + z[k]) + (j_end + z[j]) ))

        return [left_set, right_set],(j,k)
    end
end

"""
function relax_maybes_Dual(solved, problem, link_matrix, p)

For every partition p, is going to force to be zero the links maybe currentluy
used and then the dual cost is used to select the most useful
"""

function relax_maybes_Dual(solved, problem, link_matrix, p)
    num_jobs = problem.num_jobs
    jobs = problem.jobs
    x = getindex(solved, :x)

    potentialLink = findall(link_matrix[:, :] .== MAYBE)
    # Select arc (j, k) to remove to then look at dual
    linkMaybe = []
    dualConstraint = []
    for (j, k) in map(x -> x.I, potentialLink)
        if getvalue(x[p,j,k])> 0.5
            push!(linkMaybe,(j,k))
            push!(dualConstraint, @constraint(solved, x[p,j,k] <=0.0))
            break
        end
    end

    solve(solved)

    # check which one has the best dual cost
    maxLink = -1
    minShadowPrice = 1
    for (i,(j,k)) in enumerate(linkMaybe)
        if getdual(dualConstraint[i]) !=0
            println("a dual not at zero value " * string(getdual(dualConstraint[i])))
        end
        if getdual(dualConstraint[i]) < minShadowPrice
            minShadowPrice = getdual(dualConstraint[i])
            maxLink = (j,k)
        end
        # we remove the constraint of solved
        JuMP.setRHS(dualConstraint[i],1)
    end

    if maxLink == -1
        j = rand(1:num_jobs)
        k = rand(1:num_jobs)
        return [ConstraintSet()],(j,k)
    else
        (j,k) = maxLink

        # Get start and end times for these two jobs.
        j_index = jobs[j, :FLIGHT_INDEX]
        k_index = jobs[k, :FLIGHT_INDEX]
        j_end = jobs[j, :END]
        k_start = jobs[k, :START]

        # Create left and right halfspaces.
        left_set = ConstraintSet(GeneralConstraint(z -> k_start + z[k_index] - (j_end + z[j_index])  + EPS))
        right_set = ConstraintSet(GeneralConstraint(z -> -(k_start + z[k_index]) + (j_end + z[j_index])))

        return [left_set, right_set] ,(j,k)
    end
end

"""
function draw_points(partition)::Array{<:Number}

take a partition and draw some point from it
"""
function draw_points(partition,problem::Problem,link_matrix)
    num_jobs = problem.num_jobs
    num_flights = problem.num_flights
    jobs = problem.jobs

    potentialLink = findall(link_matrix[:, :] .== MAYBE)
    linkMaybe = [j for (j, k) in map(x -> x.I, potentialLink)]
    shuffle!(linkMaybe)
    unique!(linkMaybe)
    num_points = []

    # Create model.
    m = Model(solver=GurobiSolver(problem.gurobi_env, OutputFlag=0))

    # Uncertainty is constrain to lie within the partition.
    @variable(m, u[1:num_flights])
    apply_constraints!(m, u, partition)

    @objective(m, Min, 0)
    status = solve(m)
    point = getvalue(u)

    i = 1
    while (status == :Optimal) && (i <= length(linkMaybe))
        comp = 1
        nbJ = linkMaybe[i]
        j_end = jobs[nbJ, :END]
        j_index = jobs[nbJ, :FLIGHT_INDEX]
        while (status == :Optimal) && (comp <= num_jobs)
            point = getvalue(u)
            if link_matrix[nbJ,comp] == MAYBE
                job_start = jobs[comp, :START]
                job_index = jobs[comp, :FLIGHT_INDEX]
                @constraint(m,job_start + u[job_index] - (j_end + u[j_index]) <= -0.5)
            end
            status = solve(m)
            comp +=1
        end
        i +=1
    end
    push!(num_points, point)
    #@objective(m, Min, sum(-u[i] for i in 1:num_flights))
    #status = solve(m)
    #if(status == :Optimal)
        #push!(num_points, getvalue(u))
    #end

    return num_points
end

"""
function check_feasible_partition(partition,problem)
 return true is the parition is non empty
 """

function check_feasible_partition(partition,problem)
    m = Model(solver=GurobiSolver(problem.gurobi_env, OutputFlag=0))
    num_flights = problem.num_flights
    # Uncertainty is constrain to lie within the partition.
    @variable(m, u[1:num_flights])
    apply_constraints!(m, u, partition)
    @objective(m, Max,0)

    status = solve(m)

    return status == :Optimal
end
