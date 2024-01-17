#
#   Linesearch nonmonotone Grippo-Lampariello-Lucidi with spectral steplength
#

function gll(x, f, gradfx, projf, eta, fx, sigma1, sigma2, alpha_min, alpha_max, fvals, M, iter, prev_x)
    t0 = time()

    proj_grad_fk = projf(x - gradfx)

    # Compute spectral steplength para o SPG2
    #if iter==1
    #    alpha = min(alpha_max, max(alpha_min, 1.0 / norm(proj_grad_fk - x)))
    #else 
    #    alpha = 1.0
    #end

    if iter==1
        alpha = min(alpha_max, max(alpha_min, 1.0 / norm(proj_grad_fk - x)))
    else 
        # Calculate new direction
        sk = x - prev_x
        yk = gradf(x) - gradf(prev_x)

        # Update alpha for the next iteration
        if dot(sk, yk) ≤ 0
           alpha = alpha_max
        else
            alpha = min(alpha_max, max(alpha_min, dot(sk, sk) / dot(sk, yk)))
        end
    end
    
    while true
        dk = projf(x - alpha * gradfx) - x # ATENÇÃO
        GD = eta * dot(dk, gradfx)
        q = projf(x - alpha * gradfx)
        fq = f(q)

        m = min(iter-1,M-1)
        #m = min(iter-1,M)
        f_max = maximum(fvals[end - m : end])
        
        if fq ≤ f_max + GD
        #if fq ≤ f_max + alpha * GD
            et = time() - t0
            return (alpha, et)
        else 
            # Quadratic interpolation
            if alpha <= 0.1
                alpha = alpha / 2.0
                et = time() - t0
                return (alpha, et)
            else
                at = -alpha^2 * dot(gradfx, dk) / (2 * (fq - fx - alpha * dot(gradfx, dk)))
                at < sigma1 || at > sigma2 * alpha && (at = at / 2.0)
                alpha = at
                et = time() - t0
                return (alpha, et)
            end
        end
    end
end
