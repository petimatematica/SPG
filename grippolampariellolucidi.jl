#
#   Linesearch nonmonotone Grippo-Lampariello-Lucidi
#

function gll(x, f, gradfx, projf, eta, fx, sigma1, sigma2,fvals,M,iter)
    t0 = time()

    alpha = 1.0  # Inicialize α aqui

    while true
        dk = projf(x - alpha * gradfx) - x # ATENÇÃO
        GD = eta * dot(dk, gradfx)
        q = projf(x - alpha * gradfx)
        fq = f(q)
      # f_max = f(x) # Não consegui uma forma de calcular a f_max como está no algoritmo
        m = min(iter-1,M-1)
        f_max = maximum(fvals[end - m : end])
        
        if fq ≤ f_max + GD
            et = time() - t0
            return (alpha, et)
            break
        end
        
        # Quadratic interpolation
        if alpha ≤ 0.1
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
