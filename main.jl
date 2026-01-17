include("double_pendulum.jl")
include("analysis.jl")
include("helpers.jl")

#analyse video
path = "First_video_2s.mp4"
positions1_ana, positions2_ana = analyse_video(path)
#adapt analysis output to be compatible 
positions1_ana_resized = [(x/ppm, y/ppm) for (x, y) in positions1_ana]
positions2_ana_resized = [(x/ppm, y/ppm) for (x, y) in positions2_ana]

positions1_ana_f = [position(p[1], p[2]) for p in positions1_ana_resized]
positions2_ana_f = [position(p[1], p[2]) for p in positions2_ana_resized]

#define start variables-----
dt = 0.0005          # it is strongly recommended that 1/fps of reference video /dt be an Int to be able to compare on perfect time 
time = 2    #s

first_pos1_ana = positions1_ana[1]
first_pos2_ana = positions2_ana[2]

L1_pix = sqrt( (first_pos1_ana[1]-first_pos2_ana[1])^2 + (first_pos1_ana[2] - first_pos2_ana[2])^2)
L1_m = 0.09174
ppm = L1_pix/L1_m

anchor = position(0,0)
point1 = position(first_pos1_ana[1]/ppm,first_pos1_ana[2]/ppm)
point2 = position(first_pos2_ana[1]/ppm,first_pos2_ana[2]/ppm)

m1 = 10.0
m2 = 1.0

w1_init = 0
w2_init = 0

method = "euler"



#compute best parameters---------
max_tolerance = 0.03
best_params = adam_optimize(point1,point2,w1_init,w2_init,m1, m2, dt, time, method)

m1 = best_params[1]
w1_init = best_params[2]
w2_init = best_params[3]


#launch simulation-----
positions1_sim, positions2_sim = simulate(anchor, point1, point2, w1_init, w2_init, m1, m2, dt, time, method)

ratio = Int(floor(length(positions1_sim)/length(positions1_ana)))


#adapt simulated positions to be compatible with timeframe
positions1_sim_f = positions1_sim[1:ratio:end]
positions2_sim_f = positions2_sim[1:ratio:end]


println("MSE: ", compute_mse(positions2_sim_f,positions2_ana_f))
println("time-accuracy: ",time_accuracy(positions2_sim_f,positions2_ana_f,max_tolerance))


#create gifs


create_comparison_gif(positions1sim_framed,positions2sim_framed,positions1_ana_final,positions2_ana_final)

println("done")

