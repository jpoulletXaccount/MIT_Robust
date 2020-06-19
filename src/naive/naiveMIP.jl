using JuMP, Gurobi, JuMPeR, LinearAlgebra

# solve MIP and return solution

function solveNaiveMip(dictJob, dictInterJob, nbWorker, costBackUp)
    dictRelationshipInter = findRelationInter(dictInterJob)     # dict[comp] = inter
    dictRelationshipJob = findRelationJob(dictJob)              # dict[jobId] = comp

    #modelNaive = RobustModel(solver=GurobiSolver(OutputFlag=0))
    modelNaive = RobustModel(solver=GurobiSolver())

    #variables
    nbInter = length(dictInterJob)
    nbJobs = length(dictJob)
    @variable(modelNaive,r)         # represents the max in the two echelons
    @variable(modelNaive,x[1:nbWorker, 1:nbJobs],Bin)
    @variable(modelNaive,z[1:nbJobs],Bin)
    @variable(modelNaive,x_u[1:nbWorker,2:nbInter,1:nbJobs],Bin)
    @variable(modelNaive,z_u[2:nbInter,1:nbJobs],Bin)
    @variable(modelNaive,y_u[1:nbJobs,1:nbJobs],Bin)

    @uncertain(modelNaive,u[1:nbJobs])
    #@variable(modelNaive,u[1:nbJobs])

    # uncertainty set
    @constraint(modelNaive, norm(u,1) <= 100)
    @constraint(modelNaive, norm(u,Inf) <= 10)
    for i in 1:nbJobs
        @constraint(modelNaive,u[i]>=0)
    end
    println("All variables created")

    inter1 = dictRelationshipInter[1]
    # Constraint of the tasks not assigned during interval t
    for i in 1:nbWorker
        for t in 2:nbInter
            inter = dictRelationshipInter[t]
            listJobDone = [dictRelationshipJob[jobId] for jobId in dictInterJob[inter]]
            for j in 1:nbJobs
                if j in listJobDone

                else
                    @constraint(modelNaive,x_u[i,t,j] ==0)
                    @constraint(modelNaive,z_u[t,j] ==0)
                end
            end
        end

        listJobDone = [dictRelationshipJob[jobId] for jobId in dictInterJob[inter1]]
        for j in 1:nbJobs
            if j in listJobDone
            else
                @constraint(modelNaive,x[i,j] ==0)
                @constraint(modelNaive,z[j] ==0)
            end
        end
    end

    # objective function
    @objective(modelNaive,Min, r + costBackUp * sum(z[dictRelationshipJob[jobId]] for jobId in dictInterJob[inter1] ) )
    @constraint(modelNaive, r >= costBackUp * sum(z_u[t,dictRelationshipJob[jobId]] for t in 2:nbInter for jobId in dictInterJob[dictRelationshipInter[t]] ))
    println("Objective function done")

    # constraints covertness t =1
    for jobId in dictInterJob[inter1]
        j = dictRelationshipJob[jobId]
        @constraint(modelNaive, sum(x[i,j] for i in 1:nbWorker) + z[j] == 1)
    end
    println("Constraints covertness half done")

    # constraint covertness t >1
    for t in 2:nbInter
        inter = dictRelationshipInter[t]
        for jobId in dictInterJob[inter]
            j = dictRelationshipJob[jobId]
            @constraint(modelNaive, sum(x_u[i,t,j] for i in 1:nbWorker) + z_u[t,j] ==1)
        end
    end

    println("Constraints covertness done")

    # constraint chaining
    M = 100000
    for jId in keys(dictJob)
        jPos = dictRelationshipJob[jId]
        jEndTime = dictJob[jId].endTime
        for kId in keys(dictJob)
            kPos = dictRelationshipJob[kId]
            kBeginTime = dictJob[kId].beginTime
            @constraint(modelNaive, - M * y_u[jPos,kPos] <= kBeginTime + u[kPos] - (jEndTime + u[jPos]))
            #@constraint(modelNaive, kBeginTime + u[kPos] - (jEndTime + u[jPos])<= M*(1 - y_u[jPos,kPos]))
        end
    end

    println("Constraints chaining done")

    # no clash init
    for i in 1:nbWorker
        for jId in dictInterJob[inter1]
            jPos = dictRelationshipJob[jId]
            for kId in dictInterJob[inter1]
                kPos = dictRelationshipJob[kId]
                if kPos != jPos
                    @constraint(modelNaive, x[i,jPos] + x[i,kPos] <= 3 - (y_u[jPos,kPos] + y_u[kPos,jPos]))
                end
            end

            for t in 2:nbInter
                inter = dictRelationshipInter[t]
                for kId in dictInterJob[inter]
                    kPos = dictRelationshipJob[kId]
                    if kPos != jPos
                        @constraint(modelNaive,x[i,jPos] + x_u[i,t,kPos] <= 3 - (y_u[jPos,kPos] + y_u[kPos,jPos]))
                    end
                end
            end
        end
    end

    println("Constraints clash init done")

    # no clash inter
    for i in 1:nbWorker
        nbKeys = length(dictRelationshipInter)
        for comp in 2:nbKeys
            inter = dictRelationshipInter[comp]
            for compTilde in comp:nbKeys
                interTilde = dictRelationshipInter[compTilde]

                for j in dictInterJob[inter]
                    jPos = dictRelationshipJob[j]
                    for k in dictInterJob[interTilde]
                        if k!= j
                            kPos = dictRelationshipJob[k]
                            @constraint(modelNaive, x_u[i,comp,jPos] + x_u[i,compTilde,kPos] <= 3 - (y_u[jPos,kPos] + y_u[kPos,jPos]))
                        end
                    end
                end
            end
        end
    end

    println("Constraints clash inter done")

    solve(modelNaive)
    println("Model solved successfully")
    return(getvalue(z),getvalue(z_u),getvalue(x),getvalue(x_u),dictRelationshipJob,dictRelationshipInter)
end

function findRelationInter(dictInterJob)
    listInter = []
    for inter in keys(dictInterJob)
        push!(listInter,inter)
    end
    sort!(listInter)

    dictRelationship = Dict()
    comp = 1
    for inter in listInter
        dictRelationship[comp] = inter
        comp += 1
    end

    return dictRelationship
end

function findRelationJob(dictJob)
    dictRelationship = Dict()
    comp = 1
    for jobId in keys(dictJob)
        dictRelationship[jobId] = comp
        comp += 1
    end

    return dictRelationship
end
