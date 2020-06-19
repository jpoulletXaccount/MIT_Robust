using CSV
using DataFrames
using LinearAlgebra


"""
    build_jobs(path::String, intervals::Array{<:Number})::DataFrame

Load a DataFrame with jobs from file, and transforms it in the following ways:
    - Adds ID column to uniquely identify jobs.
    - Changes START and END columns into clicks for easier computation.
    - Computes an interval into which each jobs falls.
"""
function build_jobs(path::String, intervals::Array{<:Number})::DataFrame

    # Read in raw data.
    jobs = CSV.read(path)

    # Expand each job according to the number of workers it requires.
    result = DataFrame()
    for job in eachrow(jobs)
        if nrow(result) == 0
            result = repeat(DataFrame(job), job[:NUMBER])
        else
            result = vcat(result, repeat(DataFrame(job), job[:NUMBER]))
        end
    end

    # Convert to clicks.
    result[:START] = date_to_clicks.(result[:START])
    result[:END] = date_to_clicks.(result[:END])

    # Job IDs and interval labelling.
    result[:ID] = 1:nrow(result)
    result[:INTERVAL] = get_interval.(Ref(intervals), result[:START])

    # Update flight number indexing.
    flight_num_list = Int.(unique(result[:FLIGHT_NUM]))
    flight_nums = Int.(result[:FLIGHT_NUM])
    result[:FLIGHT_INDEX] = get_flight_index.(flight_nums, Ref(flight_num_list))

    return result
end


"""
    get_interval(intervals::Array{<:Number}, start::Float64)::Int

Compute an interval index based on the divisions supplied.
"""
function get_interval(intervals::Array{<:Number}, start::Float64)::Int

    # Ensure intervals are sorted.
    sorted = sort(intervals)

    # Find first index of interval start time.
    return findlast(sorted .<= start)
end


"""
    date_to_clicks(date::String, sep::String=":")::Float64

Convert a date in the format HH:MM:SS into clicks, for example:
    - 8:00:00 -> 800
    - 8:30:00 -> 850
    - 12:45:30 -> 1275.5
"""
function date_to_clicks(date::String, sep::String=":")::Float64

    # Split and convert to numeric type.
    parts = split(date, sep)
    numbers = map(s -> parse(Float64, s), parts)

    # Convert minutes and seconds to clicks.
    numbers[2:3] *= (10 / 6)

    return dot([100, 1, 0.01], numbers)
end


"""
   get_flight_index(flight_num, flight_nums)

Get the index into the list of unique flight numbers of a particular flight (makes it easier)
to index uncertainty vector.
"""
function get_flight_index(flight_num, flight_nums)
    return findfirst(isequal(flight_num), flight_nums)
end
