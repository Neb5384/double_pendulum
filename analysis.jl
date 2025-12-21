#using Pkg
# Pkg.add("ImageShow")
# Pkg.add("VideoIO")
# Pkg.add("Images")
#Pkg.add("JLD2")
#Pkg.add("FileIO")


using VideoIO, Images,ImageShow, ImageMorphology, Statistics, GLMakie, JLD2, FileIO

function get_frames(path_name)
    vid = VideoIO.openvideo(path_name)

    fps = VideoIO.framerate(vid)

    frames = []

    while !eof(vid)
        frame = VideoIO.read(vid)
        push!(frames, frame)
    end

    VideoIO.close(vid)

    return fps, frames
end

function orange_mask(img)
    img_hsv = HSV.(img)

    mask = map(c ->
        (10 <= c.h <= 35) &&   # orange hue range
        (c.s >= 0.5) &&      # saturation
        (0.4 <= c.v <= 0.7)     # brightness
    , img_hsv)


    mask_int = Int.(mask)
    mask_clean = area_opening(mask_int, min_area = 200)

    return mask_clean
end

function centroids_from_labels(labels)
    centroids = []

    for label in 1:maximum(labels)
        inds = findall(labels .== label)
        if !isempty(inds)
            rows = [Tuple(I)[1] for I in inds]
            cols = [Tuple(I)[2] for I in inds]
            push!(centroids, (mean(cols), mean(rows)))
        end
    end

    return centroids
end

function get_positions(frames)
    positions = []

    for frame in frames
        mask = orange_mask(frame)
        labels = label_components(mask)
        centroids = centroids_from_labels(labels)

        push!(positions, centroids)
    end

    return positions
end

function dist_2(pos1,pos2)
    return (pos1[1]-pos2[1])^2 + (pos1[2]-pos2[2])^2
end

function reorder_positions(positions)
    reordered = [positions[1]]  # first frame as reference

    for i in 2:length(positions)
        prev = reordered[end]
        cur = positions[i]

        # Compute distances and match by nearest neighbor
        if length(cur) == length(prev) == 2
            if dist_2(cur[1], prev[1]) + dist_2(cur[2], prev[2]) >
               dist_2(cur[2], prev[1]) + dist_2(cur[1], prev[2])
                cur = [cur[2], cur[1]]  # swap if this minimizes total distance
            end
        end
        push!(reordered, cur)
    end

    return reordered
end

function circle_center(p1, p2, p3)
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3

    denom = 2 * (x1*(y2 - y3) - x2*(y1 - y3) + x3*(y1 - y2))
    if denom == 0
        error("Points are collinear; circle center is undefined")
    end

    Cx = ((x1^2 + y1^2)*(y2 - y3) + (x2^2 + y2^2)*(y3 - y1) + (x3^2 + y3^2)*(y1 - y2)) / denom
    Cy = ((x1^2 + y1^2)*(x3 - x2) + (x2^2 + y2^2)*(x1 - x3) + (x3^2 + y3^2)*(x2 - x1)) / denom

    return (Cx, Cy)
end

function create_gif_ana(positions1,positions2)
    println("creating gif...")

    steps = length(positions1)

    max_length = 1000

    fig = Figure(size=(600,600)) 
    ax = GLMakie.Axis(fig[1,1], title="Pendulum tracker", limits = ((-max_length, max_length), (-max_length, max_length)))

    # Create observables for the plot data
    line1_obs = Observable([Point2f(0.0, 0.0), Point2f(positions1[1][1], positions1[1][2])])
    line2_obs = Observable([Point2f(positions1[1][1], positions1[1][2]), Point2f(positions2[1][1], positions2[1][2])])
    dot1_obs = Observable([Point2f(positions1[1][1], positions1[1][2])])
    dot2_obs = Observable([Point2f(positions2[1][1], positions2[1][2])])

    # Create plots using the observables
    lines!(ax, line1_obs, color=:blue)
    lines!(ax, line2_obs, color=:red)
    scatter!(ax, dot1_obs, color=:blue, markersize=10)
    scatter!(ax, dot2_obs, color=:red, markersize=10)

    # Create a GIF
    framerate = 50


    @time record(fig, "double_pendulum_fron_video.gif", 1:2:steps; framerate) do i
        p1 = positions1[i]
        p2 = positions2[i]
        
        # Update observables
        line1_obs[] = [Point2f(0.0, 0.0), Point2f(p1[1], p1[2])]
        line2_obs[] = [Point2f(p1[1], p1[2]), Point2f(p2[1], p2[2])]
        dot1_obs[] = [Point2f(p1[1], p1[2])]
        dot2_obs[] = [Point2f(p2[1], p2[2])]
    end
end


function analyse_video(path)
    println("analysing video ...")

    @time if isfile("positions.jld2")
        @load "positions.jld2" positions fps
    else
        fps, frames = get_frames(path)
        positions = get_positions(frames)
        positions = reorder_positions(positions)
        @save "positions.jld2" positions fps
    end

    center = circle_center(positions[1][2],positions[40][2],positions[60][2])

    positions_centered_flipped = map(frame -> map(pos -> (pos[1] - center[1], -(pos[2] - center[2])), frame), positions)

    positions1 = [frame[2] for frame in positions_centered_flipped]
    positions2 = [frame[1] for frame in positions_centered_flipped]

    println("analysis done")
    return positions1, positions2
end

