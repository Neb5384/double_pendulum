# using Pkg

# Pkg.add("ImageShow")
# Pkg.add("VideoIO")
# Pkg.add("Images")


using VideoIO, Images,ImageShow, ImageMorphology, Statistics

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


fps, frames = get_frames("First_video_2s.mp4")
positions = get_positions(frames)

println(fps)
println(positions)
