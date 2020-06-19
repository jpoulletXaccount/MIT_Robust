include("../../src/read.jl")
include("../../src/problem.jl")
include("../../src/utils.jl")
include("../../src/deterministic.jl")


"""
    Objective

On any particular problem size, determine the number of workers required to eliminate all
backup agents in the nominal case.
"""

function main()

    # Logging.
    out = stdout

    # Load jobs and create problem data.
    write_log(out, 0, "Loading problem.")
    jobs = build_jobs("/../../data/actual/spread_60.csv", [0])
    problem = Problem(jobs, 0, 0, 0)

    write_log(out, 0, "Solving deterministic model.")
    num_workers = solve_deterministic_model(problem)

    write_log(out, 0, "Number of workers required is $num_workers.")    

    write_log(out, 0, "Procedure finished.")
end

main()
