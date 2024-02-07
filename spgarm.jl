#
# Spectral Projected Gradient Method with Linesearch GLL
#

using LinearAlgebra, DataFrames

function spgarm(x, f, gradf, projf, alpha_min, alpha_max, epsilon, eta, maxiter, linesearch, minstep)
    
    Gnorm = Float64[];
    fvals = Float64[];
    stplen = Float64[];

    seqx=x 

    info2 = DataFrame()

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

            info2.fvals = fvals
            info2.gradnorm = Gnorm
            info2.stplen = stplen 
            et = time() - t0
            return(x,ierror,info2,et,seqx)
        end


        #dk = projf(x - alpha * gradfx) - x 

        # Update iter count
        iter += 1

        #Verifies if maxiter was achieved
        if iter > maxiter
            ierror = 1
            #println("Maximum of iterations was achieved! Stopping...")

            info2.fvals = fvals
            info2.gradnorm = Gnorm
            info2.stplen = stplen
            et = time() - t0
            return(x,ierror,info2,et,seqx)
        end
        
         # Update sequence
        if iter == 1
            (alpha,ilserror,iet) = linesearch(x, f, gradfx, projf, eta, fx, alpha_min, alpha_max, iter, seqx[end], minstep)
         push!(stplen,alpha)
        else
            (alpha,ilserror,iet) = linesearch(x, f, gradfx, projf, eta, fx, alpha_min, alpha_max, iter, seqx[:, end-1], minstep)
            push!(stplen,alpha)
        end
         
         if ilserror > 0
            ierror = 2
            #println("Step lenght toot small!")

            push!(fvals,NaN)
            info2.fvals = fvals
            push!(Gnorm,NaN)
            info2.gradnorm = Gnorm
            info2.stplen = stplen
            #info.inttime = inttime
            et = time() - t0
            return(x,ierror,info2,et,seqx)
        end

        xnproj = x - alpha * gradfx
        newx = projf(x - alpha * gradfx)
        proj_grad_fk = projf(x - gradfx)

        x=newx

        seqx=[seqx x];
    end

end
