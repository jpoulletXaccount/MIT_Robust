using JuMP
using JuMPeR
using LinearAlgebra


##### TYPES #####

"""
    Constraint

Supertype for the different constraints that can be applied to either uncertain parameters
or variables through a ConstraintSet.
"""
abstract type Constraint end


"""
    ConstraintSet

Holds constraints on uncertain parameters and variables.

A ConstraintSet should be defined for a particular uncertain parameter or variable (or array of
these), and all constraints should be designed to act on just this object.
"""
struct ConstraintSet

    # Constraints are stored in an array.
    constraints::Array{Constraint}

    # Basic constructors for convenience.
    ConstraintSet() = new([])
    ConstraintSet(constraint::C) where {C<:Constraint} = new([constraint])
    ConstraintSet(constraints::Array{<:Constraint}) = new(constraints)
end


"""
    GeneralConstraint

Defines a general constraint with a function of the form f(x) <= 0.
"""
struct GeneralConstraint <: Constraint
    f::Function
end


"""
    NormConstraint

Defines a norm constraint. We need this as a wrapper since JuMP doesn't support norm constraints
in the form @constraint(m, norm(u, p) <= Ω) like JuMPeR does.
"""
struct NormConstraint <: Constraint
    norm::Union{Int,Float64}
    bound::Float64
end


##### METHODS ####

"""
    add_constraint(cs::ConstraintSet, c::C) where {C<:Constraint}

Add a constraint to a ConstraintSet.
"""
function add_constraint!(cs::ConstraintSet, c::C) where {C<:Constraint}
    push!(cs.constraints, c)
end


"""
    apply_constraints!(m::Model,
                       u::Union{AS, Array{AS}},
                       cs::ConstraintSet) where {AS<:JuMP.AbstractJuMPScalar}

Apply the constraints contained in a ConstraintSet to a particular uncertain parameter or
variable (or array of these).
"""
function apply_constraints!(m::Model,
                            u::Union{AS, Array{AS}},
                            cs::ConstraintSet) where {AS<:JuMP.AbstractJuMPScalar}
    for constraint in cs.constraints
        apply_constraint!(m, u, constraint)
    end
end


"""
    apply_constraint!(m::Model,
                      u::Union{AS, Array{AS}},
                      g::GeneralConstraint) where {AS<:JuMP.AbstractJuMPScalar}

Apply a GeneralConstraint.
"""
function apply_constraint!(m::Model,
                           u::Union{AS, Array{AS}},
                           g::GeneralConstraint) where {AS<:JuMP.AbstractJuMPScalar}
    apply_general_constraint!(m, g.f(u))
end


"""
    apply_constraint!(m::Model, u::Array{JuMP.Variable}, n::NormConstraint)

Apply a NormConstraint to an array of variables (norms only make sense for arrays).
"""
function apply_constraint!(m::Model, u::Array{JuMP.Variable}, n::NormConstraint)

    # Define auxiliary variable and relevant constraints.
    v = @variable(m, [1:length(u)])
    @constraint(m, v .>= u)
    @constraint(m, v .>= -u)
    @constraint(m, sum(v) <= n.bound ^ n.norm)
end


"""
    apply_constraint!(m::Model, u::Array{JuMPeR.Uncertain}, n::NormConstraint)

Add a NormConstraint for an array of uncertain parameters (norms only make sense for arrays).
"""
function apply_constraint!(m::Model, u::Array{JuMPeR.Uncertain}, n::NormConstraint)

    # Directly use JuMPeR's functionality.
    @constraint(m, norm(u, n.norm) <= n.bound)
end


"""
    apply_general_constraint!(m::Model, lhs::AS) where {AS<:JuMP.AbstractJuMPScalar}

Add a GeneralConstraint where the function evaluates to a scalar expression.
"""
function apply_general_constraint!(m::Model, lhs::AS) where {AS<:JuMP.AbstractJuMPScalar}
    @constraint(m, lhs <= 0)
end


"""
    apply_general_constraint!(m::Model, lhs::Array{<:JuMP.AbstractJuMPScalar})

Add a GeneralConstraint where the function evaluates to an array expression.
"""
function apply_general_constraint!(m::Model, lhs::Array{<:JuMP.AbstractJuMPScalar})
    @constraint(m, lhs .<= 0)
end
