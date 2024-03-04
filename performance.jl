#
# Spectral Projected Gradient Method - SPG Method
#
using CUTEst, NLPModels, LinearAlgebra, DataFrames, Random, Printf, Plots, BenchmarkProfiles, JLD2

include("spgperformance.jl")

sproblem= "problem"

V = Float64[]
T = Float64[]
S = Float64[]

# Defina uma matriz de nomes de problemas CUTEst
#problems = ["BDEXP", "EXPLIN", "EXPLIN2"]#, "EXPQUAD", "PROBPENL", "S368",
#"HADAMALS", "LINVERSE", "NONSCOMP",
#"TORSION1","TORSION2","TORSION3","TORSION4","TORSION5","TORSION6",
#"TORSIONA","TORSIONB","TORSIONC","TORSIOND","TORSIONE","TORSIONF"]

#dimension = ["5000", "120", "120"]#, "120", "500", "100",
#"1024", "1999", "10000",
#"14884", "14884", "14884", "14884", "14884", "14884", 
#"14884", "14884", "14884", "14884", "14884", "14884"]

# Defina uma matriz de nomes de problemas CUTEst
#problems = ["TORSION1","TORSION2","TORSION3","TORSION4","TORSION5","TORSION6",
#"TORSIONA","TORSIONB","TORSIONC","TORSIOND","TORSIONE","TORSIONF"]

#dimension = ["14884", "14884", "14884", "14884", "14884", "14884",
#"14884", "14884", "14884", "14884", "14884", "14884"]

#custom_filter = x -> x["origin"] == "real"
#custom_filter = x -> x["origin"] == "academic"
#custom_filter = x -> x["origin"] == "modelling"
#problems = CUTEst.select(objtype="quadratic", contype="bounds", custom_filter=custom_filter)
#problems = CUTEst.select(objtype = "sum_of_squares", contype = "bounds", custom_filter=custom_filter)

problems1 = CUTEst.select(objtype="quadratic", contype="bounds")#, custom_filter=custom_filter)
problems2 = CUTEst.select(objtype = "sum_of_squares", contype = "bounds")#, custom_filter=custom_filter)

problems = vcat(problems1, problems2)

for B in 1:2 
    if B == 1 
        Ls = "SPG1"
        println("SPG1")
        linesearch = backtracking1
    else 
        Ls = "SPG2"
        println("SPG2")
        linesearch = backtracking2
    end

    for ip in 1:length(problems)

        # Inicialize um vetor de modelos CUTEst
        nlp = CUTEstModel(problems[ip])#, "-param", "N="*dimension[i])
    
        # Initial guess from CUTEst
        x0 = nlp.meta.x0
    
        # Objective functions
    
        global function f(x)
            return obj(nlp,x)  # Objective function from CUTEst
        end
    
        # Gradient of Objective function
    
        global function gradf(x)
            return grad(nlp,x) # Gradient of objective function (CUTEst)
        end
    
        #
        #  Upper and lower bounds seting
        #
        l = Array{Float64}(undef,size(x0))
        u = Array{Float64}(undef,size(x0))
        for i in 1 : size(x0,1)
            global l[i] = -100.0
            global u[i] = 50.0
        end
        #
        # Orthogonal projection
        #
        global function proj(x)
            n = size(x,1)
            z = Array{Float64}(undef,size(x0))
    
            for i in 1:n
                z[i] = max(l[i],min(x[i],u[i]))
            end
    
            return z
        end
    
        #
        #  Solver parameter seting and calling
        #
        tol = 1.e-5
        maxiter = 10000
        alpha_min = 1.e-30
        alpha_max = 1.e+30
        M = 10
        sigma1 = 0.1
        sigma2 = 0.9
        gamma = 1.e-4

        println(problems[ip])
        println(length(x0))
        
        (x,error,info,seqx,etime,evalf) = spgperformance(x0, f, gradf, proj, tol, maxiter, alpha_min, alpha_max, M, sigma1, sigma2, gamma, linesearch);
        
        filename = "echo/" * sproblem * string(ip) * Ls * ".jld2"
        @save filename info 

        if error > 0
            push!(V, Inf)
            push!(T, Inf)
            push!(S, Inf)
        else
            iters = size(seqx, 2)
            push!(V, iters)
            push!(T, etime)
            push!(S, evalf)
        end   

        finalize(nlp)
    end
end

h=length(problems);
W=[V[1:h] V[h+1:2h]]; #Matriz com iteradas
Z=[T[1:h] T[h+1:2h]]; #Matriz com os tempos
R=[S[1:h] S[h+1:2h]]; #Matriz com as avaliações de função

colors=[:royalblue1, :green2]

X = performance_profile(PlotsBackend(), W, ["SPG1", "SPG2"], xlabel="Performance ratio: # iterations", ylabel="Solved problems [%]", legend=:bottomright, palette=colors, lw=2, dpi=1000)
Y = performance_profile(PlotsBackend(), Z, ["SPG1", "SPG2"], xlabel="CPU time ratio", ylabel="Solved problems [%]", legend=:bottomright, palette=colors, lw=2, dpi=1000)
Q = performance_profile(PlotsBackend(), R, ["SPG1", "SPG2"], xlabel="Function evaluation", ylabel="Solved problems [%]", legend=:bottomright, palette=colors, lw=2, dpi=1000)

plot(X)

savefig("performanceprofileiters")

plot(Y)

savefig("performanceprofiletime")

plot(Q)

savefig("performanceprofileevalf")
