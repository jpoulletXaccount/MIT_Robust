include("creater.jl")
include("naiveMIP.jl")
include("shiftBuilder.jl")


dictJob = createJobs("../../data/runway_1_shift_57.csv")
listInterval = [800,1200,1800,2400]   # correspond to the interval that we are going to consider
dictInterJob = createDictIntervalJob(listInterval , dictJob)

z,z_u,x,x_u,dictRelationshipJob,dictRelationshipInter = solveNaiveMip(dictJob, dictInterJob, 10, 1)
analyzeBackup(z,z_u)
recreateAndPrint(x,x_u,dictRelationshipJob,dictRelationshipInter,dictJob)

println("End of run")
