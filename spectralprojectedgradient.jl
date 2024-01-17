#
# Spectral Projected Gradient Method with Linesearch
#

using LinearAlgebra, DataFrames

function spg(x, f, gradf, projf, alpha_min, alpha_max, epsilon, eta, M, maxiter, linesearch)
    
    Gnorm = Float64[];
    fvals = Float64[];
    stplen = Float64[];

    seqx=x 

    info = DataFrame()

    ierror = 0 # Number of errors
    iter = 0 # Number of iterations

    push!(stplen, NaN)

    t0 = time()

    sigma1 = 0.1
    sigma2 = 0.9

    
    while true
        
        fx = f(x)
        push!(fvals, fx)

        gradfx = gradf(x)
        normgradfx = norm(gradfx, 2)
        push!(Gnorm, normgradfx)

        proj_grad_fk = projf(x - gradfx)
        
        # Compute spectral steplength
        #if iter==0
        #alpha = min(alpha_max, max(alpha_min, 1.0 / norm(proj_grad_fk - x)))
        #end

        distance = norm(proj_grad_fk - x)
        # Verifies if convergence was achieved
        if distance ≤ epsilon
            #println("Solution has found!")

            info.fvals = fvals
            info.gradnorm = Gnorm
            info.stplen = stplen 
            et = time() - t0
            return(x,ierror,info,et,seqx)
        end


        #dk = projf(x - alpha * gradfx) - x 

        # Update iter count
        iter += 1

        #Verifies if maxiter was achieved
        if iter > maxiter
            ierror = 1
            #println("Maximum of iterations was achieved! Stopping...")

            info.fvals = fvals
            info.gradnorm = Gnorm
            info.stplen = stplen
            et = time() - t0
            return(x,ierror,info,et,seqx)
        end
        
         # Update sequence
        if iter == 1
            (alpha,iet) = linesearch(x, f, gradfx, projf, eta, fx, sigma1, sigma2, alpha_min, alpha_max, fvals, M, iter, seqx[end])
         push!(stplen,alpha)
        else
            (alpha,iet) = linesearch(x, f, gradfx, projf, eta, fx, sigma1, sigma2, alpha_min, alpha_max, fvals, M, iter, seqx[:, end-1])
            push!(stplen,alpha)
        end
         

        xnproj = x - alpha * gradfx
        newx = projf(x - alpha * gradfx)
        proj_grad_fk = projf(x - gradfx)

        x=newx

        seqx=[seqx x];
    end

end
