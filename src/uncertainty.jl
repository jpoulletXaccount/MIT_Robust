include("problem.jl")
include("constraints.jl")


function standard_uncertainty(problem::Problem)

	# Create constraints that define a standard uncertainty set.
	base = ConstraintSet()

	# u <= 60 minutes.
    add_constraint!(base, GeneralConstraint(u -> u - 100))

    # u >= -30 minutes.
    add_constraint!(base, GeneralConstraint(u -> -u - 50))

    # norm(u, 1) <= 10 minutes * num_flights.
    add_constraint!(base, NormConstraint(1, (30 / 6) * problem.num_flights))
	#add_constraint!(base, NormConstraint(1, 150))

    return base
end
