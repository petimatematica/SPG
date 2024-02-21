#
# Spectral Projected Gradient Method - SPG Method
#
using CUTEst, NLPModels, LinearAlgebra
include("spg.jl")

# Loading objective function from CUTEst
nlp = CUTEstModel("CVXBQP1", "-param", "N=10000")

# Initial guess from CUTEst
x0 = nlp.meta.x0

# Objective function
#
function f(x)
    return obj(nlp,x)  # Objective function from CUTEst
end
#
# Gradient of Objective function
#
function gradf(x)
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
function proj(x)
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
linesearch = backtracking1

sol,error = spg()

finalize(nlp)