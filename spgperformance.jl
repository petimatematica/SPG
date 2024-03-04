#
# Spectral Projected Gradient Method
#
function spgperformance(x0, f, gradf, proj, tol, maxiter, alpha_min, alpha_max, M, sigma1, sigma2, gamma, linesearch)

    stplen = Float64[];
    fvals = Float64[];
    Gnorm = Float64[];

    info = DataFrame()

    f_hist = Float64[]
    feval = Float64[]

    x = copy(x0)
    seqx=x
    
    push!(f_hist,f(x))

    # Set number of iterations
    iter = 0

    push!(stplen,NaN)

    gradf_x = gradf(x)
    
    alpha = 1.0 / norm(gradf_x,2)

    t0 = time()

    while true
        P_x_m_grad = proj(x-gradf_x)
        norm_P = norm(P_x_m_grad - x,2)
        norm_G1_inf =  norm(P_x_m_grad - x,Inf)
        
        fx = f(x)
        push!(fvals,fx)

        normgradfx = norm(gradf_x,2)
        push!(Gnorm,normgradfx)

        # print info code 
        #println("iter = $iter  norm_G1_inf = $norm_G1_inf  f(x) = $(f(x))")

        # Detect whether the current point is stationary
        if norm_G1_inf < tol
            # print info code 
            #println("iter = $iter  norm_G1_inf = $norm_G1_inf  f(x) = $(f(x))")
            #println("Solutions has found!")

            info.stplen = stplen 
            info.fvals = fvals
            info.gradnorm = Gnorm

            error=0

            et = time() - t0
            
            #evalf_k = feval[end:end]
            #evalf_k = float(feval[end])
            evalf_k = sum(feval)

            println("iter = $iter  norm_G1_inf = $norm_G1_inf  f(x) = $(f(x)) tempo = $et evalf = $evalf_k")
            println("Solutions has found!")

            #return x,0
            return(x,error,info,seqx,et,evalf_k)
        end

        # Update number of iterations
        iter += 1

        #Detect whether the maximum of iterations was achieved
        if iter > maxiter
            println("Maximum of iterations was achieved! Stoping...")

            info.stplen = stplen 
            info.fvals = fvals
            info.gradnorm = Gnorm

            error=1 

            et = time() - t0
            
            #evalf_k = feval[end:end]
            #evalf_k = float(feval[end])
            evalf_k = sum(feval)

            #return x,1
            return(x,error,info,seqx,et,evalf_k)
        end


        # Backtrackin routine
        (x,gradf_x,s,y,f_hist,lambda,et,evalf) = linesearch(iter,alpha,x,gradf_x,f_hist, M, sigma1, sigma2, gamma)
        push!(stplen,lambda)
        push!(feval,evalf)

        # Step 3
        b = dot(s,y)
        if b > 0.0
            a = dot(s,s)
            alpha = min(alpha_max,max(alpha_min,a/b))
        else
            alpha = alpha_max               
        end

        seqx=[seqx x];
    end
end

#
# Backtrackin routine
#
function backtracking1(k,alpha_k,x_k,gradf_x_k,f_hist, M, sigma1, sigma2, gamma)
    lambda = copy(alpha_k)
    t0 = time()
    evalf = 0  # Initialize evalf counter

    while true
        x_plus = proj(x_k - lambda * gradf_x_k)
        #println("$x_plus")
        m_k = min(k-1,M-1)
        f_max = maximum(f_hist[end-m_k:end])
        f_x_plus = f(x_plus)
        test = f_x_plus > f_max + gamma * dot(x_plus - x_k,gradf_x_k)
        if ~test
            s_k = x_plus - x_k
            gradf_x_kp1 = gradf(x_plus)
            y_k = gradf_x_kp1 - gradf_x_k
            push!(f_hist,f_x_plus)
            et = time() - t0
            evalf += 1  # Increment evalf counter
            return (x_plus,gradf_x_kp1,s_k,y_k,f_hist,lambda,et,evalf)
        else
            fxx = f(x_k + lambda * (x_plus - x_k))
            #lambda = lambda / 2.0 
            if lambda <= 0.1
                lambda = lambda / 2.0
            else
                atemp = (- dot(x_plus - x_k,gradf_x_k) * lambda^2) / (2.0 * (fxx - f(x_k) - lambda * dot(x_plus - x_k,gradf_x_k)))
                #if atemp < sigma1 *lambda || atemp > sigma2 * lambda
                if atemp < sigma1 || atemp > sigma2 * lambda
                    atemp = lambda / 2.0
                end
                lambda = atemp
            end
        end 
    end
end

#
# Backtrackin routine SPG2
#
function backtracking2(k,alpha_k,x_k,gradf_x_k,f_hist, M, sigma1, sigma2, gamma)
    alpha = copy(alpha_k)
    lambda = 1
    t0 = time()
    evalf = 0  # Initialize evalf counter

    while true
        d_k = proj(x_k - alpha * gradf_x_k) - x_k
        x_plus = x_k + lambda * d_k
        #println("$x_plus")
        m_k = min(k-1,M-1)
        f_max = maximum(f_hist[end-m_k:end])
        f_x_plus = f(x_plus)
        test = f_x_plus > f_max + gamma * lambda * dot(d_k,gradf_x_k)
        if ~test
            s_k = x_plus - x_k
            gradf_x_kp1 = gradf(x_plus)
            y_k = gradf_x_kp1 - gradf_x_k
            push!(f_hist,f_x_plus)
            et = time() - t0
            evalf += 1  # Increment evalf counter
            return (x_plus,gradf_x_kp1,s_k,y_k,f_hist,lambda,et,evalf)
        else
            fxx = f(x_k + lambda * (x_plus - x_k))
            #lambda = lambda / 2.0 
            if lambda <= 0.1
                lambda = lambda / 2.0
            else
                atemp = (- dot(x_plus - x_k,gradf_x_k) * lambda^2) / (2.0 * (fxx - f(x_k) - lambda * dot(x_plus - x_k,gradf_x_k)))
                #if atemp < sigma1 *lambda || atemp > sigma2 * lambda
                if atemp < sigma1 || atemp > sigma2 * lambda
                    atemp = lambda / 2.0
                end
                lambda = atemp
            end
        end 
    end
end