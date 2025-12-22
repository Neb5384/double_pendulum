include("double_pendulum.jl")
include("analysis.jl")


#analyse video
path = "First_video_2s.mp4"
positions1_ana, positions2_ana = analyse_video(path)



#define start variables-----
dt = 0.0000005
time = 2


first_pos1_ana = positions1_ana[1]
first_pos2_ana = positions2_ana[2]

L1_pix = sqrt( (first_pos1_ana[1]-first_pos2_ana[1])^2 + (first_pos1_ana[2] - first_pos2_ana[2])^2)
L1_m = 0.09174
pixelsperm = L1_pix/L1_m

println(pixelsperm)

anchor = position(0,0)
point1 = position(first_pos1_ana[1]/pixelsperm,first_pos1_ana[2]/pixelsperm)
point2 = position(first_pos2_ana[1]/pixelsperm,first_pos2_ana[2]/pixelsperm)

m1 = 30
m2 = 2

w1_init = 0
w2_init = 0

method = "euler"

#launch simulation-----
positions1_sim, positions2_sim = simulate(anchor, point1, point2, w1_init, w2_init, m1, m2, dt, time,method)


#create gifs

create_gif_sim(positions1_sim,positions2_sim)
create_gif_ana(positions1_ana, positions2_ana)

println("done")

