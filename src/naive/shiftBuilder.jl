EPS = 0.001
# File use to build the shift and print then
function recreateAndPrint(x,x_u,dictRelationshipJob,dictRelationshipInter,dictJob)
    dictWorkerJob = builtShift(x,x_u,dictRelationshipJob,dictRelationshipInter,dictJob)
    #printShift(dictWorkerJob,dictJob)
    covertnessCheck(z,z_u,x,x_u)
end

function builtShift(x,x_u,dictRelationshipJob,dictRelationshipInter,dictJob)
    nbWorker = length(x[:,1])
    nbInter = length(dictRelationshipInter)
    nbJobs = length(dictJob)
    dictWorkerJob = Dict()

    nbTasksDoneAgent = 0

    for i in 1:nbWorker
        listJob = []
        for j in 1:nbJobs
            if x[i,j] >1 - EPS
                nbTasksDoneAgent +=1
                jobId = findJobId(j,dictRelationshipJob)
                push!(listJob,jobId)
            end
            for t in 2:nbInter
                if x_u[i,t,j] > 1 - EPS
                    nbTasksDoneAgent +=1
                    jobId = findJobId(j,dictRelationshipJob)
                    push!(listJob,jobId)
                end
            end
        end
        dictWorkerJob[i] = listJob
    end
    println("Nb taks done by agents  " * string(nbTasksDoneAgent))
    return dictWorkerJob
end


function findJobId(j,dictRelationshipJob)
    for (key,value) in dictRelationshipJob
        if value == j
            return key
        end
    end
    println("Seems to be a big mistake")
    return ""
end

function printShift(dictWorkerJob,dictJob)
    for i in keys(dictWorkerJob)
        dictJobBeginTime = Dict()
        listBeginTime = []
        for jobId in dictWorkerJob[i]
            dictJobBeginTime[dictJob[jobId].beginTime] = jobId
            push!(listBeginTime,dictJob[jobId].beginTime)
        end
        sort!(listBeginTime)
        txt = ""
        for beginT in listBeginTime
            txt *= dictJobBeginTime[beginT] * " -> "
        end
        println(string(i) * " does job " * txt)
    end
end


function analyzeBackup(z,z_u)
    nbJobs = length(z)
    nbInter = length(z_u[:,1])
    nbBackup = 0
    for j in 1:nbJobs
        if z[j] > 1 - EPS
            nbBackup +=1
        end
    end
    for t in 2:nbInter
        for j in 1:nbJobs
            if z_u[t,j] > 1 - EPS
                nbBackup +=1
            end
        end
    end
    println("Nb tasks performed by backup " * string(nbBackup))
end

function covertnessCheck(z,z_u,x,x_u)
    nbJobs = length(z)
    nbWorker = length(x[:,1])
    nbInter = length(z_u[:,1])
    for j in 1:nbJobs
        hasBeenDone = false
        for i in 1:nbWorker
            if x[i,j] > 1 - EPS
                hasBeenDone= true
            end

            for t in 2:nbInter
                if x_u[i,t,j] > 1 - EPS
                    hasBeenDone = true
                end
            end
        end
        if z[j] > 1 - EPS
            hasBeenDone =true
        end
        for t in 2:nbInter
            if z_u[t,j] > 1 - EPS
                hasBeenDone = true
            end
        end

        if hasBeenDone == false
            println("Error one task is not covered ... " * string(j))
        end
    end
end
