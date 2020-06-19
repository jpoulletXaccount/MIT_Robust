using Dates


function write_log(stream, iteration::Int, message)
    if iteration < 1
        it_string = ""
    else
        it_string = "(Iteration $iteration) "
    end
    write(stream, "$(Dates.format(now(), "HH:MM:SS")) | $it_string$message\n")
    return nothing
end
