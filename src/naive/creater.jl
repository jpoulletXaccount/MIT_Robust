include("object.jl")

function createJobs(nameFile)
    dictJob = Dict()
    open(nameFile) do f
        for line in enumerate(eachline(f))
            if line[1] > 1
                token = split(line[2],",")

                # We need a few modifications before being able to compute the clicks
                newtoken = split(token[6], ":")
                beginT = 100 * parse(Int,newtoken[1]) + floor((10/6) * parse(Int,newtoken[2]))
                newtoken = split(token[7], ":")
                endT = 100 * parse(Int,newtoken[1]) + floor((10/6) * parse(Int,newtoken[2]))

                nbAgents = parse(Int,token[5])
                for i in 1:nbAgents
                    stringId = string(line[1]) * "_" * string(i)
                    jobCreated = JobAF(stringId, beginT, endT ,parse(Int,token[2]))
                    dictJob[jobCreated.id] = jobCreated
                end

            end
        end
    end
    println(" We have finished to create all jobs: " * string(length(dictJob)))
    return(dictJob)
end


function createDictIntervalJob(listInterval,dictJob)
    beginInter = 0
    dictInterJob = Dict()
    for inter in listInterval
        listJob = []
        for (key,value) in dictJob
            considerBegin = value.beginTime
            if beginInter <= considerBegin < inter
                push!(listJob,key)
            end
        end
        dictInterJob[inter] = listJob
        println(" For interval: [" * string(beginInter) * "," * string(inter) * " ], we have " * string(length(listJob)) * " jobs ")
        beginInter = inter
    end

    return dictInterJob
end
