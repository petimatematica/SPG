#
# Spectral Projected Gradient Method
#
function spg()

    f_hist = Float64[]
    x = copy(x0)
    
    push!(f_hist,f(x))

    # Set number of iterations
    iter = 0
    gradf_x = gradf(x)
    alpha = 1.0 / norm(gradf_x,2)
    while true
        P_x_m_grad = proj(x-gradf_x)
        norm_P = norm(P_x_m_grad - x,2)
        norm_G1_inf =  norm(P_x_m_grad - x,Inf)

        # print info code 
        println("iter = $iter  norm_G1_inf = $norm_G1_inf  f(x) = $(f(x))")

        # Detect whether the current point is stationary
        if norm_G1_inf < tol
            println("Solutions has found!")
            return x,0
        end

        # Update number of iterations
        iter += 1

        #Detect whether the maximum of iterations was achieved
        if iter > maxiter
            println("Maximum of iterations was achieved! Stoping...")
            return x,1
        end

        # Backtrackin routine
        x,gradf_x,s,y,f_hist = linesearch(iter,alpha,x,gradf_x,f_hist)

        # Step 3
        b = dot(s,y)
        if b > 0.0
            a = dot(s,s)
            alpha = min(alpha_max,max(alpha_min,a/b))
        else
            alpha = alpha_max               
        end
    end
end

#
# Backtrackin routine
#
function backtracking1(k,alpha_k,x_k,gradf_x_k,f_hist)
    lambda = copy(alpha_k)

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
            return x_plus,gradf_x_kp1,s_k,y_k,f_hist
        else
            #lambda = lambda / 2.0 
            if lambda <= 0.1
                lambda = lambda / 2.0
            else
                atemp = (- dot(x_plus - x_k,gradf_x_k) * lambda^2) / (2.0 * (f_x_plus - f(x_k) - lambda * dot(x_plus - x_k,gradf_x_k)))
                if atemp < sigma1 || atemp > sigma2 * lambda
                    atemp = lambda / 2.0
                end
                lambda = atemp
            end
        end 
    end
end


#
# Backtrackin routine avarege
#
function backtracking2(k,alpha_k,x_k,gradf_x_k,f_hist)
    lambda = copy(alpha_k)

    while true
        x_plus = proj(x_k - lambda * gradf_x_k)
        #println("$x_plus")
        m_k = min(k-1,M-1)
        f_med = sum(f_hist[end-m_k:end])/m_k
        f_x_plus = f(x_plus)
        test = f_x_plus > f_med + gamma * dot(x_plus - x_k,gradf_x_k)
        if ~test
            s_k = x_plus - x_k
            gradf_x_kp1 = gradf(x_plus)
            y_k = gradf_x_kp1 - gradf_x_k
            push!(f_hist,f_x_plus)
            return x_plus,gradf_x_kp1,s_k,y_k,f_hist
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
# Backtrackin routine Armijo
#
function armijo(k, alpha_k, x_k, gradf_x_k, f_hist)
    #lambda = 1.0 
    lambda = copy(alpha_k)

    minstep = 1.e-6
    
    while true
        GD = gamma * dot(gradf_x_k, proj(x_k - lambda * gradf_x_k) - x_k)

        x_plus = proj(x_k - lambda * gradf_x_k)
        f_x_plus = f(x_plus)
        gradf_x_kp1 = gradf(x_plus)

        stptest = f_x_plus - f(x_k) + GD
        
        if stptest > 0.0
            #lambda = 0.5 * lambda

            atemp = (- dot(x_plus - x_k,gradf_x_k) * lambda^2) / (2.0 * (f_x_plus - f(x_k) - lambda * dot(x_plus - x_k,gradf_x_k)))
            if atemp < sigma1 || atemp > sigma2 * lambda
                atemp = lambda / 2.0
            end
            lambda = atemp

            if lambda < minstep
                println("minstep")
                return x_plus, gradf_x_kp1
                break
            end
        else
            s_k = x_plus - x_k
            gradf_x_kp1 = gradf(x_plus)
            y_k = gradf_x_kp1 - gradf_x_k
            push!(f_hist,f_x_plus)
            return x_plus,gradf_x_kp1,s_k,y_k,f_hist
        end        
    end
end
