#
# Perfomance profiles of Spectral Projected Gradient Method (SPG)
#

using CUTEst, NLPModels, LinearAlgebra, DataFrames, Random, Printf, Plots, BenchmarkProfiles, JLD2

include("spg.jl")

sproblem= "problem"

V = Float64[]
T = Float64[]
S = Float64[]

problems = ["BDEXP", "EXPLIN", "EXPLIN2", "EXPQUAD", "PROBPENL", "S368", 
"HADAMALS", "CHEBYQAD", "HS110", "LINVERSE", "NONSCOMP", "QR3DLS", 
"DECONVB", "BIGGSB1", "BQPGABIM", "BQPGASIM", "JNLBRNG1", "JNLBRNGA",  
"NCVXBQP1", "NOBNDTOR", "PENTDI", "TORSION1", "TORSION2", "TORSION3", 
"TORSION4", "TORSION5", "TORSION6", "TORSIONA", "TORSIONB", "TORSIONC", 
"TORSIOND", "TORSIONE", "TORSIONF"]


dimension = ["5000", "120", "120", "120", "500", "100",
"1024", "50", "50", "1999", "10000", "610", 
"61", "1000", "50", "50", "15625", "15625",
"10000", "14884","1000", "14884", "14884", "14884",
"14884", "14884", "14884", "14884", "14884", "14884", 
"14884","14884", "14884"]

#length(problems)

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

        nlp = CUTEstModel(problems[ip], "-param", "N="*dimension[ip])
    
        # Initial guess from CUTEst
        x0 = nlp.meta.x0
    
        # Objective functions from CUTEst
        global function f(x)
            return obj(nlp,x) 
        end
        
        # Gradient of Objective function from CUTEst
        global function gradf(x)
            return grad(nlp,x)
        end
        

        #  Upper and lower bounds seting
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
        
        #  Solver parameter seting and calling
        tol = 1.e-5
        maxiter = 10000
        lambda_min = 1.e-30
        lambda_max = 1.e+30
        M = 10
        sigma1 = 0.1
        sigma2 = 0.9
        gamma = 1.e-4

        println(problems[ip])
        println(length(x0))

        (x,error,info,seqx,etime,evalf) = spg(x0, f, gradf, proj, tol, maxiter, lambda_min, lambda_max, M, sigma1, sigma2, gamma, linesearch);
            
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
W=[V[1:h] V[h+1:2h]]; #Matrix which stores iterations
Z=[T[1:h] T[h+1:2h]]; #Matrix which stores CPU time
R=[S[1:h] S[h+1:2h]]; #Matrix which stores function evaluation

colors=[:royalblue1, :green2]

X = performance_profile(PlotsBackend(), W, ["SPG1", "SPG2"], xlabel="Number of iterations", ylabel="Solved problems [%]", legend=:bottomright, palette=colors, lw=2, dpi=1000)
Y = performance_profile(PlotsBackend(), Z, ["SPG1", "SPG2"], xlabel="CPU time ratio", ylabel="Solved problems [%]", legend=:bottomright, palette=colors, lw=2, dpi=1000)
Q = performance_profile(PlotsBackend(), R, ["SPG1", "SPG2"], xlabel="Function evaluation", ylabel="Solved problems [%]", legend=:bottomright, palette=colors, lw=2, dpi=1000)

plot(X)
savefig("performanceprofileiters")

plot(Y)
savefig("performanceprofiletime")

plot(Q)
savefig("performanceprofileevalf")