#imports-------------
# using Pkg
# Pkg.add("GLMakie")
using GLMakie

#constants-----------------------
g = 9.80665

#structs-------------------
mutable struct position
	x::Float64
	y::Float64
end

function distance(pos1::position,pos2::position)
    return sqrt((pos1.x - pos2.x)^2 + (pos1.y -pos2.y)^2)
end

function squared_distance(pos1::position,pos2::position)
    return (pos1.x - pos2.x)^2 + (pos1.y -pos2.y)^2
end

#helper funcs----------------
function get_alpha1(g,L1,L2,m1,m2,w1,w2,theta1,theta2)
    delta = theta2 - theta1
    numerator = m2*L1*w1^2*sin(delta)*cos(delta) + m2*g*sin(theta2)*cos(delta) + m2*L2*w2^2*sin(delta) - (m1+m2)*g*sin(theta1)
    denominator = (m1+m2)*L1 -m2*L1*cos(delta)^2
    return alpha1 = numerator / denominator
end

function get_alpha2(g,L1,L2,m1,m2,w1,w2,theta1,theta2)
    delta = theta2 - theta1
    numerator = -1*m2*L2*w2^2*sin(delta)*cos(delta) + (m1+m2)*( g*sin(theta1)*cos(delta) - L1*w1^2*sin(delta) - g * sin(theta2))
    denominator = (m1+m2)*L2 -m2*L2*cos(delta)^2
    return alpha2 = numerator / denominator
end

function cartesian_pos(theta1, theta2, L1, L2)
    pos1_x = L1 * sin(theta1)
    pos1_y = - L1 * cos(theta1)
    pos2_x = pos1_x + L2 * sin(theta2)
    pos2_y = pos1_y - L2 * cos(theta2)

    return position(pos1_x,pos1_y), position(pos2_x, pos2_y)
end

function derivatives(g,L1,L2,m1,m2,theta1,theta2,w1,w2)
    alpha1 = get_alpha1(g,L1,L2,m1,m2,w1,w2,theta1,theta2)
    alpha2 = get_alpha2(g,L1,L2,m1,m2,w1,w2,theta1,theta2)
    return w1, w2, alpha1, alpha2
end

#main funcs--------
function euler_step(g,L1,L2,m1,m2,w1,w2,theta1,theta2,dt)

    alpha1 = get_alpha1(g,L1,L2,m1,m2,w1,w2,theta1,theta2)
    alpha2 = get_alpha2(g,L1,L2,m1,m2,w1,w2,theta1,theta2)

    w1_next = w1 + alpha1 * dt
    w2_next = w2 + alpha2 * dt

    theta1_next = theta1 + w1_next * dt
    theta2_next = theta2 + w2_next * dt

    return w1_next, w2_next, theta1_next, theta2_next

end

function rk4_step(g,L1,L2,m1,m2,w1,w2,theta1,theta2,dt)
    # k1
    dtheta1_1, dtheta2_1, dw1_1, dw2_1 = derivatives(g,L1,L2,m1,m2,theta1,theta2,w1,w2)
    
    # k2
    dtheta1_2, dtheta2_2, dw1_2, dw2_2 = derivatives(
        g,L1,L2,m1,m2,
        theta1 + 0.5*dt*dtheta1_1,
        theta2 + 0.5*dt*dtheta2_1,
        w1 + 0.5*dt*dw1_1,
        w2 + 0.5*dt*dw2_1
    )
    
    # k3
    dtheta1_3, dtheta2_3, dw1_3, dw2_3 = derivatives(
        g,L1,L2,m1,m2,
        theta1 + 0.5*dt*dtheta1_2,
        theta2 + 0.5*dt*dtheta2_2,
        w1 + 0.5*dt*dw1_2,
        w2 + 0.5*dt*dw2_2
    )
    
    # k4
    dtheta1_4, dtheta2_4, dw1_4, dw2_4 = derivatives(
        g,L1,L2,m1,m2,
        theta1 + dt*dtheta1_3,
        theta2 + dt*dtheta2_3,
        w1 + dt*dw1_3,
        w2 + dt*dw2_3
    )
    
    # combine increments
    theta1_next = theta1 + dt*(dtheta1_1 + 2*dtheta1_2 + 2*dtheta1_3 + dtheta1_4)/6
    theta2_next = theta2 + dt*(dtheta2_1 + 2*dtheta2_2 + 2*dtheta2_3 + dtheta2_4)/6
    w1_next = w1 + dt*(dw1_1 + 2*dw1_2 + 2*dw1_3 + dw1_4)/6
    w2_next = w2 + dt*(dw2_1 + 2*dw2_2 + 2*dw2_3 + dw2_4)/6
    
    return w1_next, w2_next, theta1_next, theta2_next
end

#simulation-------------------


function simulate(anchor, point1, point2, w1, w2, m1, m2, dt, time,method, verbose = false)
    
    steps = Int(time / dt)
    verbose && println(steps)

    L1 = distance(anchor,point1)
    L2 = distance(point1,point2)

    theta1 = atan(point1.x - anchor.x, - point1.y + anchor.y)
    theta2 = atan(point2.x - point1.x, - point2.y + point1.y)

    p1check, p2check = cartesian_pos(theta1, theta2, L1, L2)

    theta1s = [theta1]
    theta2s = [theta2]

    verbose && println("simulating physics....")

    if method == "rk4"
        for i in 2:steps  #Julia counting starts at 1. and we already have qst step, hence 2
            w1, w2, theta1, theta2 = rk4_step(g,L1,L2,m1,m2,w1,w2,theta1,theta2,dt)
            append!(theta1s,theta1)
            append!(theta2s,theta2)
        end
    else
        for i in 2:steps
            w1, w2, theta1, theta2 = euler_step(g,L1,L2,m1,m2,w1,w2,theta1,theta2,dt)
            append!(theta1s,theta1)
            append!(theta2s,theta2)
        end
    end


    positions1 = [] 
    positions2 = []  

    verbose && println("changing to cartesian coordinates...")

    for (th1, th2) in zip(theta1s, theta2s)
        p1, p2 = cartesian_pos(th1, th2, L1, L2)
        push!(positions1, p1)
        push!(positions2, p2)
    end

    verbose && println("simulation done ")

    return positions1, positions2

end

function create_gif_sim(positions1,positions2)
    println("creating gif...")

    steps = length(positions1)

    anchor = position(0,0)
    L1 = distance(anchor,positions1[1])
    L2 = distance(point1,positions2[2])
    max_length = L1 + L2

    fig = Figure(size=(600,600)) 
    ax = GLMakie.Axis(fig[1,1], title="Pendulum Simulation", limits = ((-max_length, max_length), (-max_length, max_length)))

    # Create observables for the plot data
    line1_obs = Observable([Point2f(0.0, 0.0), Point2f(positions1[1].x, positions1[1].y)])
    line2_obs = Observable([Point2f(positions1[1].x, positions1[1].y), Point2f(positions2[1].x, positions2[1].y)])
    dot1_obs = Observable([Point2f(positions1[1].x, positions1[1].y)])
    dot2_obs = Observable([Point2f(positions2[1].x, positions2[1].y)])

    # Create plots using the observables
    lines!(ax, line1_obs, color=:blue)
    lines!(ax, line2_obs, color=:red)
    scatter!(ax, dot1_obs, color=:blue, markersize=10)
    scatter!(ax, dot2_obs, color=:red, markersize=10)

    # Create a GIF
    framerate = 50
    interval = max(floor(Int, 1/dt/framerate),1)


    @time record(fig, "double_pendulum.gif", 1:interval:steps; framerate) do i
        p1 = positions1[i]
        p2 = positions2[i]
        
        # Update observables
        line1_obs[] = [Point2f(0.0, 0.0), Point2f(p1.x, p1.y)]
        line2_obs[] = [Point2f(p1.x, p1.y), Point2f(p2.x, p2.y)]
        dot1_obs[] = [Point2f(p1.x, p1.y)]
        dot2_obs[] = [Point2f(p2.x, p2.y)]
    end
end


