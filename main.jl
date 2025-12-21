include("double_pendulum.jl")
#include("analysis.jl")


#define start variables-----
dt = 0.0000005
time = 5

anchor = position(0,0)
point1 = position(0.5,0.5)
point2 = position(01.2,1.2)

m1 = 2
m2 = 3

w1_init = 0
w2_init = 0

method = "euler"

#launch simulation-----
positions1, positions2 = simulate(anchor, point1, point2, w1_init, w2_init, m1, m2, dt, time,method)




create_gif(positions1,positions2)

