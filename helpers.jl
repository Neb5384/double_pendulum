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


#functions needed to calculate speed at end of video
function calculate_angle(position1, position2)
    x = position2.x - position1.x
    y = - position2.y + position1.y

    return atan(x, y)
end

function calculate_w(pivot_pos1, pivot_pos2, pos1, pos2, fps, step)
    angle_1 = calculate_angle(pivot_pos1, pos1)
    angle_2 = calculate_angle(pivot_pos2, pos2)
    angle_diff = angle_2 - angle_1
    angle_diff = mod(angle_diff + pi, 2pi) - pi
    w = angle_diff * fps / step

    return w
end


#function to show pendulum video extension
function create_extended_overlay_gif(positions1_extended, positions2_extended, video_end_frame, fps=50, filename="extended_overlay.gif")
    println("Creating extended overlay gif...")
    
    steps = length(positions1_extended)
    
    # Calculate max length
    anchor = position(0, 0)
    L1 = distance(anchor, positions1_extended[1])
    L2 = distance(positions1_extended[1], positions2_extended[1])
    max_length = L1 + L2
    
    # Create title observable FIRST
    title_obs = Observable("Tracking Video (1/$(video_end_frame))")
    
    fig = Figure(size=(600, 600)) 
    # Pass the observable when creating the axis
    ax = GLMakie.Axis(fig[1, 1], title=title_obs, 
                      limits=((-max_length, max_length), (-max_length, max_length)))
    
    # Observables for positions
    line1_obs = Observable([Point2f(0.0, 0.0), Point2f(positions1_extended[1].x, positions1_extended[1].y)])
    line2_obs = Observable([Point2f(positions1_extended[1].x, positions1_extended[1].y), Point2f(positions2_extended[1].x, positions2_extended[1].y)])
    dot1_obs = Observable([Point2f(positions1_extended[1].x, positions1_extended[1].y)])
    dot2_obs = Observable([Point2f(positions2_extended[1].x, positions2_extended[1].y)])
    
    # Create plots
    lines!(ax, line1_obs, color=:blue, linewidth=2)
    lines!(ax, line2_obs, color=:red, linewidth=2)
    scatter!(ax, dot1_obs, color=:blue, markersize=15)
    scatter!(ax, dot2_obs, color=:red, markersize=15)
    
    @time record(fig, filename, 1:steps; framerate=fps) do i
        p1 = positions1_extended[i]
        p2 = positions2_extended[i]
        
        # Update positions
        line1_obs[] = [Point2f(0.0, 0.0), Point2f(p1.x, p1.y)]
        line2_obs[] = [Point2f(p1.x, p1.y), Point2f(p2.x, p2.y)]
        dot1_obs[] = [Point2f(p1.x, p1.y)]
        dot2_obs[] = [Point2f(p2.x, p2.y)]
        
        # Update title when transitioning from video to prediction
        if i <= video_end_frame
            title_obs[] = "Tracking Video ($(i)/$(video_end_frame))"
        else
            title_obs[] = "Prediction ($(i - video_end_frame)/$(steps - video_end_frame))"
        end
    end
    
    println("Extended overlay gif created: $filename")
end


#calculate energies
function kinetic_energy(m1, m2, w1, w2, pos1, pos2)

    # velocity of mass 1
    v1x = -w1 * pos1.y
    v1y =  w1 * pos1.x
    v1_squared = v1x^2 + v1y^2

    # relative position of mass 2 to mass 1
    rx = pos2.x - pos1.x
    ry = pos2.y - pos1.y

    # velocity of mass 2
    v2x = v1x - w2 * ry
    v2y = v1y + w2 * rx
    v2_squared = v2x^2 + v2y^2

    return 0.5*m1*v1_squared + 0.5*m2*v2_squared
end

function potential_energy(m1,m2,pos1,pos2,g) # mgh
    potential_energy = m1 * g * pos1.y + m2 * g * pos2.y
    return potential_energy
end

function total_energy(m1,m2,w1,w2,L1,L2,pos1,pos2,g)
    total_energy = kinetic_energy(m1,m2,w1,w2, pos1, pos2) + potential_energy(m1,m2,pos1,pos2,g)
    return total_energy
end

#apply energy calculations on list
function compute_energy_trajectory(positions1, positions2, w1s, w2s, L1, L2, m1, m2, g)

    kinetic_energies = []
    potential_energies = []
    total_energies = []
    
    for i in eachindex(positions1)
        ke = kinetic_energy(m1,m2,w1s[i],w2s[i],positions1[i],positions2[i])
        pe = potential_energy(m1,m2,positions1[i],positions2[i],g)
        te = ke + pe
        push!(kinetic_energies, ke)
        push!(potential_energies, pe)
        push!(total_energies, te)
    end
    
    return kinetic_energies, potential_energies, total_energies
end


#plot energies
function plot_energy(time,dt, kinetic_energies, potential_energies, total_energies, 
                     filename="energy_plot.png")

    println("Creating energy plot...")
    times = collect(0:dt:(time-dt))

    fig = Figure(size=(1000, 600))
    ax = GLMakie.Axis(fig[1, 1], 
              xlabel="Time (s)", 
              ylabel="Energy (J)",
              title="Double Pendulum Energy")
    
    lines!(ax, times, kinetic_energies, label="Kinetic Energy", color=:blue, linewidth=2)
    lines!(ax, times, potential_energies, label="Potential Energy", color=:red, linewidth=2)
    lines!(ax, times, total_energies, label="Total Energy", color=:green, linewidth=2)
    
    axislegend(ax, position=:rt)
    
    save(filename, fig)
    println("Energy plot saved: $filename")
    
    return fig
end