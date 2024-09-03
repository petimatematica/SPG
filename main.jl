#
# Perfomance profiles of Spectral Projected Gradient Method (SPG)
#

using CUTEst, NLPModels, LinearAlgebra, DataFrames, Random, Printf, Plots, BenchmarkProfiles, JLD2

include("spg.jl")

V = Float64[]
T = Float64[]
S = Float64[]
G = Float64[]

# Parameters
gamma = 1.e-4 
ε = 1.e-5 
min_step = 1.e-5
max_iter = 15000
lambda_min = 1.e-30
lambda_max = 1.e+30
M = 10
sigma1 = 0.1
sigma2 = 0.9

problems = ["BDEXP", "EXPLIN", "EXPLIN2", "EXPQUAD", "PROBPENL", "S368",
"HADAMALS", "NONSCOMP", "DECONVB", "BIGGSB1", "BQPGABIM", "BQPGASIM", 
"BQPGAUSS", "JNLBRNG1", "JNLBRNG2", "JNLBRNGA", "NCVXBQP1", "NCVXBQP2",
"NCVXBQP3", "OBSTCLAL", "OBSTCLBL", "OBSTCLBM", "OBSTCLBU", "PENTDI",
"LINVERSE", "NOBNDTOR", "TORSION1", "TORSION2", "TORSION3", "TORSION4", 
"TORSION5", "TORSION6", "TORSIONA", "TORSIONB", "TORSIONC", "TORSIOND", 
"TORSIONE", "TORSIONF"]

dimension = ["5000", "120", "120", "120", "500", "100", 
"1024", "10000", "63", "1000", "50", "50", 
"2003", "10000", "10000", "10000","10000", "10000", 
"10000", "10000", "10000", "10000", "10000", "1000",
"1000", "61", "61", "61", "61", "61", 
"61", "61", "61", "61", "61", "61", 
"61", "61"] 


for B in 1:2 
    if B == 1 
        Ls = "SPG1"
        println("SPG1")
        linesearch = spg1
    else 
        Ls = "SPG2"
        println("SPG2")
        linesearch = spg2
    end

    for ip in 1:length(problems)

             if ip in 1:24
             nlp = CUTEstModel(problems[ip], "-param", "N="*dimension[ip])
             elseif ip in 25:length(problems)
                 nlp = CUTEstModel(problems[ip], "-param", "Q="*dimension[ip])
             end

            println(problems[ip])
            println(dimension[ip])
            
            # Initial guess from CUTEst
            x0 = nlp.meta.x0
            global x0
        
            # Objective functions from CUTEst
            global function f(x)
                return obj(nlp, x) 
            end
            
            # Gradient of Objective function from CUTEst
            global function ∇f(x)
                return grad(nlp, x)
            end

            # Upper and lower bounds setting
            l = Array{Float64}(undef,size(x0))
            u = Array{Float64}(undef,size(x0))
            for i in 1 : size(x0,1)
                global l[i] = -100.0
                global u[i] = 50.0
            end

            # Orthogonal projection
            global function proj(x)
                n = size(x,1)
                z = Array{Float64}(undef,size(x0))
            
                for i in 1:n
                    z[i] = max(l[i],min(x[i],u[i]))
                end
                return z
            end

            (x,error,info,seqx,et,evalf,evalsproj) = spg(x0, f, ∇f, proj, ε, max_iter, lambda_min, lambda_max, M, sigma1, sigma2, gamma, linesearch)

            filename = "echo/" * problems[ip] * Ls * ".jld2"
            @save filename info 

            if error > 0
                push!(V, Inf)
                push!(T, Inf)
                push!(S, Inf)
                push!(G, Inf)
            else
                iters = size(seqx, 2)
                push!(V, iters)
                push!(T, et)
                push!(S, evalf)
                push!(G, evalsproj)
            end   

            finalize(nlp)
        end
end

ENV["GKSwstype"] = "100"

h = length(problems)
W=[V[1:h] V[h+1:2h]]; #Matrix which stores iterations
Z=[T[1:h] T[h+1:2h]]; #Matrix which stores CPU time
R=[S[1:h] S[h+1:2h]]; #Matrix which stores function evaluation
E=[G[1:h] G[h+1:2h]]; #Matrix which stores projection evaluation

colors=[:blue, :green2]

X = performance_profile(PlotsBackend(), W, ["SPG1", "SPG2"], xlabel="Performance ratio: # iterations", ylabel="Solved problems [%]", legend=:bottomright, palette=colors, lw=1.5, dpi=1000)
Y = performance_profile(PlotsBackend(), Z, ["SPG1", "SPG2"], xlabel="Performance ratio: CPU time", ylabel="Solved problems [%]", legend=:bottomright, palette=colors, lw=1.5, dpi=1000)
Q = performance_profile(PlotsBackend(), R, ["SPG1", "SPG2"], xlabel="Performance ratio: # function evaluations", ylabel="Solved problems [%]", legend=:bottomright, palette=colors, lw=1.5, dpi=1000)
N = performance_profile(PlotsBackend(), E, ["SPG1", "SPG2"], xlabel="Performance ratio: # projection evaluations", ylabel="Solved problems [%]", legend=:bottomright, palette=colors, lw=1.5, dpi=1000)

p = plot(X)
savefig(p, "performanceprofileiters.png") 

q = plot(Y)
savefig(q,"performanceprofiletime.png")

r = plot(Q)
savefig(r,"performanceprofileevalf.png")

n = plot(N)
savefig(n,"performanceprofileevalproj.png")
