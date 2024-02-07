#
# Test with Price-Dixon function
#

include("spectralprojectedgradient.jl");
include("grippolampariellolucidi.jl");
include("spgarm.jl");
include("armijo.jl");

using LinearAlgebra, DataFrames, Random, Printf, Plots, CUTEst, NLPModels

#nlp=CUTEstModel("PRICE3")
#x0 = [1.0; 2.0]

#nlp=CUTEstModel("ROSENBR")
#x0 = nlp.meta.x0

nlp=CUTEstModel("FLETCHCR")
x0 = nlp.meta.x0

function f(x)
    return obj(nlp,x)
end

function gradf(x)
    return grad(nlp,x)
end

function ortogonal(n)
    a = Float64[]
    push!(a, -1/sqrt(2))
    push!(a, 1)
    for i in 3:n
        push!(a, 0)
    end 
    return a
end

function projf(x)
    n=length(x)
    a=ortogonal(n)
    b=0
    px=x+(b-dot(a,x))*a/norm(a)^2
    return px
end

(x,ierror,info1,etime,seqx) = spg(x0, f, gradf, projf, 1.e-30, 1.e+30, 1.e-6, 1.e-4, 10, 3000, gll)

(x,ierror,info2,etime,seqx) = spgarm(x0, f, gradf, projf, 1.e-30, 1.e+30, 1.e-8, 1.e-4, 3000, armijo, 1.e-6)

#println("$info")

finalize(nlp)
