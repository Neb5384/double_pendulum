#comparison gif------
function create_comparison_gif(positions1_sim, positions2_sim, positions1_ana, positions2_ana, fps=50, filename="pendulum_comparison.gif")
    println("Creating comparison gif...")
    
    steps = min(length(positions1_sim), length(positions1_ana))
    
    # Calculate max length
    anchor = position(0, 0)
    L1_sim = distance(anchor, positions1_sim[1])
    L2_sim = distance(positions1_sim[1], positions2_sim[1])
    max_length = L1_sim + L2_sim
    println(max_length,L1_sim,L2_sim)
    
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