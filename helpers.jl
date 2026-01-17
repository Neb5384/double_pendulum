#long functions I want to invoke in main here to make main less clustered

#comparison gif------
function create_comparison_gif(positions1_sim, positions2_sim, positions1_ana, positions2_ana, fps=50, filename="pendulum_comparison.gif")
    println("Creating comparison gif...")
    
    steps = min(length(positions1_sim), length(positions1_ana))
    
    # Calculate max length
    anchor = position(0, 0)
    L1_sim = distance(anchor, positions1_sim[1])
    L2_sim = distance(positions1_sim[1], positions2_sim[1])
    max_length = L1_sim + L2_sim
    
    fig = Figure(size=(1200, 600)) 
    
    # Left panel: Simulation
    ax1 = GLMakie.Axis(fig[1, 1], title="Simulation", 
                       limits=((-max_length, max_length), (-max_length, max_length)))
    
    # Right panel: Analysis (from video)
    ax2 = GLMakie.Axis(fig[1, 2], title="Video Analysis", 
                       limits=((-max_length, max_length), (-max_length, max_length)))
    
    # Simulation observables
    line1_sim_obs = Observable([Point2f(0.0, 0.0), Point2f(positions1_sim[1].x, positions1_sim[1].y)])
    line2_sim_obs = Observable([Point2f(positions1_sim[1].x, positions1_sim[1].y), Point2f(positions2_sim[1].x, positions2_sim[1].y)])
    dot1_sim_obs = Observable([Point2f(positions1_sim[1].x, positions1_sim[1].y)])
    dot2_sim_obs = Observable([Point2f(positions2_sim[1].x, positions2_sim[1].y)])
    
    # Analysis observables
    line1_ana_obs = Observable([Point2f(0.0, 0.0), Point2f(positions1_ana[1].x, positions1_ana[1].y)])
    line2_ana_obs = Observable([Point2f(positions1_ana[1].x, positions1_ana[1].y), Point2f(positions2_ana[1].x, positions2_ana[1].y)])
    dot1_ana_obs = Observable([Point2f(positions1_ana[1].x, positions1_ana[1].y)])
    dot2_ana_obs = Observable([Point2f(positions2_ana[1].x, positions2_ana[1].y)])
    
    # Plot simulation
    lines!(ax1, line1_sim_obs, color=:blue, linewidth=2)
    lines!(ax1, line2_sim_obs, color=:red, linewidth=2)
    scatter!(ax1, dot1_sim_obs, color=:blue, markersize=15)
    scatter!(ax1, dot2_sim_obs, color=:red, markersize=15)
    
    # Plot analysis
    lines!(ax2, line1_ana_obs, color=:blue, linewidth=2)
    lines!(ax2, line2_ana_obs, color=:red, linewidth=2)
    scatter!(ax2, dot1_ana_obs, color=:blue, markersize=15)
    scatter!(ax2, dot2_ana_obs, color=:red, markersize=15)
    
    @time record(fig, filename, 1:steps; framerate=fps) do i
        # Update simulation
        p1_sim = positions1_sim[i]
        p2_sim = positions2_sim[i]
        line1_sim_obs[] = [Point2f(0.0, 0.0), Point2f(p1_sim.x, p1_sim.y)]
        line2_sim_obs[] = [Point2f(p1_sim.x, p1_sim.y), Point2f(p2_sim.x, p2_sim.y)]
        dot1_sim_obs[] = [Point2f(p1_sim.x, p1_sim.y)]
        dot2_sim_obs[] = [Point2f(p2_sim.x, p2_sim.y)]
        
        # Update analysis
        p1_ana = positions1_ana[i]
        p2_ana = positions2_ana[i]
        line1_ana_obs[] = [Point2f(0.0, 0.0), Point2f(p1_ana.x, p1_ana.y)]
        line2_ana_obs[] = [Point2f(p1_ana.x, p1_ana.y), Point2f(p2_ana.x, p2_ana.y)]
        dot1_ana_obs[] = [Point2f(p1_ana.x, p1_ana.y)]
        dot2_ana_obs[] = [Point2f(p2_ana.x, p2_ana.y)]
    end
    
    println("Comparison gif created: $filename")
end

#calculate-time-accuracy----
function time_accuracy(positions_sim, positions_ana, max_tolerance)
    for i in eachindex(positions_sim)
        diff = distance(positions_sim[i], positions_ana[i])
        if diff > max_tolerance
            return i-1
        end
    end
    return i
end

#calculate-mse-----
function compute_mse(positions_sim, positions_ana)
    total_error = 0.0
    
    for i in eachindex(positions_sim)
        total_error += squared_distance(positions_sim[i], positions_ana[i])
    end
    
    return total_error / length(positions_sim)
end
#create a hybrid scoring system--
function hybrid_score(positions_sim, positions_ana, max_tolerance)
    time_score = time_accuracy(positions_sim, positions_ana, max_tolerance)
    mse_score = compute_mse(positions_sim, positions_ana)
    final_score = time_score - mse_score
    return final_score
end

#optimise-weights-andinitial-speed------------
function run_simulation(p)
    m1_test = p[1]
    w1_test = p[2]
    w2_test = p[3]
    
    positions1_sim, positions2_sim = simulate(
        anchor, point1, point2, 
        w1_test, w2_test, m1_test, m2, 
        dt, time, method
    )
    # Downsample to match analysis frame rate
    ratio = Int(floor(length(positions1_sim) / length(positions1_ana)))
    positions1_sim_framed = positions1_sim[1:ratio:end]
    positions2_sim_framed = positions2_sim[1:ratio:end]
    
    # Convert to position type for comparison
    positions1_sim_f = [position(p.x, p.y) for p in positions1_sim_framed]
    positions2_sim_f = [position(p.x, p.y) for p in positions2_sim_framed]
    
    return positions2_sim_f
end


function compute_numerical_gradient(p)
    grad = zeros(length(p))
    base_score = hybrid_score(run_simulation(p), positions2_ana_f,max_tolerance)

    epsilons = [0.5, 0.05, 0.05]# I need mass to change more to feel its effects so it gets its own epsilon
    for i in 1:length(p)
        p_plus = copy(p)
        p_plus[i] += epsilons[i]
        
        score_plus = hybrid_score(run_simulation(p_plus), positions2_ana_f,max_tolerance)
        grad[i] = (score_plus - base_score) / epsilons[i]
    end

    return grad
end

function adam_optimize(point1, point2, w1_init, w2_init, m1, m2, dt, time, method,
        max_iterations = 50,
        learning_rate = 0.3,
        beta1 = 0.9,
        beta2 = 0.999,
        epsilon_adam = 1e-8)

    params = [m1, w1_init, w2_init]


    m = zeros(length(params))  # First moment 
    v = zeros(length(params))  # Second moment 
    best_params = copy(params)
    best_score = -Inf

    println("Starting Adam Optimization")

    for iter in 1:max_iterations

        gradient = compute_numerical_gradient(params)
        
        # Update biased first moment 
        m = beta1 * m + (1 - beta1) * gradient
        # Update biased second moment 
        v = beta2 * v + (1 - beta2) * (gradient .^ 2)
        
        # Compute bias-corrected moment 
        m_hat = m / (1 - beta1^iter)
        v_hat = v / (1 - beta2^iter)
        
        # Update prams
        params = params + learning_rate * m_hat ./ (sqrt.(v_hat) .+ epsilon_adam)
        

        # Evaluate current params
        current_score = hybrid_score(run_simulation(params), positions2_ana_f,max_tolerance)
        
        # Track best params
        if current_score > best_score
            best_score = current_score
            best_params = copy(params)
            println("[V]Iter $iter: score = $(round(current_score, digits=4)), params = $(round.(params, digits=3))")
        else
            println("   Iter $iter: score = $(round(current_score, digits=4)), params = $(round.(params, digits=3))")
        end
    end
    return best_params
end