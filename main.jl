include("double_pendulum.jl")
include("analysis.jl")
include("helpers.jl")

#analyse video
path = "First_video_2s.mp4"
positions1_ana, positions2_ana = analyse_video(path)

#adapt analysis output to be compatible 
first_pos1_ana = positions1_ana[1]
first_pos2_ana = positions2_ana[2]


L1_pix = sqrt( (first_pos1_ana[1])^2 + (first_pos2_ana[2])^2)
L1_m = 0.09174
ppm = L1_pix/L1_m

positions1_ana_resized = [(x/ppm, y/ppm) for (x, y) in positions1_ana]
positions2_ana_resized = [(x/ppm, y/ppm) for (x, y) in positions2_ana]

positions1_ana_f = [position(p[1], p[2]) for p in positions1_ana_resized]
positions2_ana_f = [position(p[1], p[2]) for p in positions2_ana_resized]

#define start variables-----
dt = 0.0005          # it is strongly recommended that 1/fps of reference video /dt be an Int to be able to compare on perfect time 
time = 2    #s

anchor = position(0,0)
point1 = position(first_pos1_ana[1]/ppm,first_pos1_ana[2]/ppm)
point2 = position(first_pos2_ana[1]/ppm,first_pos2_ana[2]/ppm)

m1 = 5.0
m2 = 1.0

w1_init = 0
w2_init = 0

method = "rk4"



#compute best parameters---------
max_tolerance = 0.03
best_params = adam_optimize(point1,point2,w1_init,w2_init,m1, m2, dt, time, method)

m1 = best_params[1]
w1_init = best_params[2]
w2_init = best_params[3]


#launch simulation-----
positions1_sim, positions2_sim, w1s, w2s = simulate(anchor, point1, point2, w1_init, w2_init, m1, m2, dt, time, method)

ratio = Int(floor(length(positions1_sim)/length(positions1_ana)))


#adapt simulated positions to be compatible with timeframe
positions1_sim_f = positions1_sim[1:ratio:end]
positions2_sim_f = positions2_sim[1:ratio:end]


println("MSE: ", compute_mse(positions2_sim_f,positions2_ana_f))
println("time-accuracy: ",time_accuracy(positions2_sim_f,positions2_ana_f,max_tolerance)," out of ",length(positions1_sim_f))

#calculate energies
kinetic_energies, potential_energies, total_energies = compute_energy_trajectory(positions1_sim, positions2_sim, w1s, w2s, m1, m2, g)
energy_drift = total_energies[end]-total_energies[1]
println("energy drift: ", energy_drift)
max_energy_diff = maximum(total_energies)-minimum(total_energies)
energy_accuaracy = energy_diff/total_energies[1]*100
println("worst energy accuracy percentage : ", energy_accuaracy)
energy_standard_dev = standard_deviation(total_energies)
println("energy standard deviation : ", energy_standard_dev)

plot_energy(time,dt, kinetic_energies, potential_energies, total_energies)

#compute start parameters from video end----
fps = 100
halfstep = 8
stepsize = 2*halfstep
last_pos1_ana = positions1_ana_f[end]
last_pos2_ana = positions2_ana_f[end]
tolast_pos1_ana = positions1_ana_f[end-stepsize]
tolast_pos2_ana = positions2_ana_f[end-stepsize]
midlast_pos1_ana = positions1_ana_f[end-halfstep]
midlast_pos2_ana = positions2_ana_f[end-halfstep]

w1_init = calculate_w(anchor,anchor,tolast_pos1_ana,last_pos1_ana,fps,stepsize)
w2_init = calculate_w(tolast_pos1_ana,last_pos1_ana,tolast_pos2_ana,last_pos2_ana,fps,stepsize)

positions1_extend, positions2_extend = simulate(anchor, midlast_pos1_ana, midlast_pos2_ana, w1_init, w2_init, m1, m2, dt, 1+halfstep/fps, method)


#adapt simulated positions to be compatible with timeframe
positions1_extend_f = positions1_extend[ratio*halfstep:ratio:end]
positions2_extend_f = positions2_extend[ratio*halfstep:ratio:end]

#extend existing tracking
positions1_extended_f = vcat(positions1_ana_f, positions1_extend_f)
positions2_extended_f = vcat(positions2_ana_f, positions2_extend_f)

#create gifs---------

#comparison
create_comparison_gif(positions1_sim_f,positions2_sim_f,positions1_ana_f,positions2_ana_f)

#video_extension
video_end_frame = length(positions1_ana_f)
create_extended_overlay_gif(positions1_extended_f, positions2_extended_f, video_end_frame)

println("done")

