using JuMP, JuMPeR, Gurobi

include("constraints.jl")


"""
    build_model(base::ConstraintSet, partitions::Array{ConstraintSet})

Builds a toy example to test how a model can be built with partitioned uncertainty.
"""
function build_model(base::ConstraintSet, partitions::Array{ConstraintSet})

    # Create base model.
    model = RobustModel(solver=GurobiSolver(OutputFlag=0))

    ## Static

    @variable(model, x >= 0)
    @variable(model, s)

    @objective(model, Min, s)

    ## Adjustable

    num_partitions = length(partitions)

    # Define copies of uncertainty and adjustable variables, indexed by partition.
    @uncertain(model, z[1:num_partitions, 1:2])
    @variable(model, y[1:num_partitions], Bin)

    for (i, partition) in enumerate(partitions)

        @constraint(model, s >= x + y[i])

        # Apply both base and partition cosntraints to get the final uncertainty set.
        apply_constraints(model, z[i, :], base)
        apply_constraints(model, z[i, :], partition)

        @constraint(model, x + y[i] >= z[i, 1] - z[i, 2])
        @constraint(model, x - y[i] <= z[i, 1] + z[i, 2])
    end

    return model
end


##### SCRIPT #####

# Box constraints on z.
base = ConstraintSet()
add_constraint!(base, GeneralConstraint(z -> -z))
add_constraint!(base, GeneralConstraint(z -> z - 1))

# Split the box in two.
low = ConstraintSet(GeneralConstraint(z -> z[1] + z[2] - 1))
high = ConstraintSet(GeneralConstraint(z -> 1 - z[1] - z[2]))

partitions = [low, high]
model = build_model(base, partitions)
