using DataFrames
using Gurobi


"""
	Problem

Defines a data instance for a problem, which helps to group everything together and
avoid passing lists of parameters around functions.
"""
struct Problem

	# Supplied data and parameters.
	jobs::DataFrame
	num_workers::Int
	cost_backup::Float64
	M::Float64
	gurobi_env::Gurobi.Env

	# Need to be computed from jobs DataFrame.
	num_jobs::Int
	num_intervals::Int
	num_flights::Int

	function Problem(jobs, num_workers, cost_backup, M)

		env = Gurobi.Env()
<<<<<<< HEAD
		setparam!(env, "TimeLimit", 15 * 60)
		#setparam!(env, "OutputFlag", 0)
=======
		# setparam!(env, "TimeLimit", 15 * 60)
		setparam!(env, "OutputFlag", 0)
>>>>>>> e7e5eb9683cecb558052fbb4181a18066ec06eb1

		return new(jobs, num_workers, cost_backup, M, env, nrow(jobs),
			maximum(jobs[:INTERVAL]), length(unique(jobs[:FLIGHT_NUM])))
	end
end


"""
	jobs_in_interval(jobs::DataFrame, interval::Int)::Array{Int}

Return the list of job IDs which fall into a particular interval.
"""
function jobs_in_interval(jobs::DataFrame, interval::Int)::Array{Int}
    return filter(row -> row[:INTERVAL] == interval, jobs)[:ID]
end


"""
	jobs_not_in_interval(jobs::DataFrame, interval::Int)::Array{Int}

Return the list of job IDs which do not fall into a particular interval.
"""
function jobs_not_in_interval(jobs::DataFrame, interval::Int)::Array{Int}
    return filter(row -> row[:INTERVAL] != interval, jobs)[:ID]
end
