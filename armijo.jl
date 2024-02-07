#
#   Armijo's Linesearch (for problem constrained)
#
function armijo(x, f, gradfx, projf, eta, fx, alpha_min, alpha_max, iter, prev_x, minstep)
    ierror = 0
    t0 = time()

    proj_grad_fk = projf(x - gradfx)

    #α = 1.0  # Inicialize α aqui

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
        GD = eta * dot(gradfx, projf(x - alpha * gradfx) - x)

        q = projf(x - alpha * gradfx)
        fq = f(q)
        stptest = fq - fx + GD

        if stptest > 0.0
            alpha = 0.5 * alpha
            if alpha < minstep
                ierror = 1
                et = time() - t0
                return (alpha, ierror, et)
            end
        else
            et = time() - t0
            return (alpha, ierror, et)
        end
    end
end