
# using Pkg

# Pkg.add("ImageShow")
# Pkg.add("VideoIO")
# Pkg.add("Images")
# Pkg.add("Colors")
# Pkg.add("ImageMorphology")
# Pkg.add("ImageSegmentation") 
# Pkg.add("Statistics") 


using VideoIO, Images, Colors,ImageMorphology, ImageSegmentation,ImageShow, Statistics

# Open video
vid = VideoIO.openvideo("First_video_2s.mp4")

frames = []

while !eof(vid)
    frame = VideoIO.read(vid)
    push!(frames, frame)
end


VideoIO.close(vid)

function orange_mask(img)
    img_hsv = HSV.(img)

    mask = map(c ->
        (10 ≤ c.h ≤ 50) &&   # orange hue range
        (c.s ≥ 0.2) &&      # saturation
        (c.v ≥ 0.2)         # brightness
    , img_hsv)

    #mask = BitMatrix(mask)
    return mask
end

mask = orange_mask(frames[1])
overlay = colorview(RGB,
    channelview(RGB.(frames[1])) .* reshape(mask, 1, size(mask)...)
)


function find_orange_centroids(mask)
    labels = label_components(mask)
    num_labels = maximum(labels)

    centroids = []

    for label in 1:num_labels
        coords = findall(labels .== label)
        isempty(coords) && continue

        ys = getindex.(coords, 1)
        xs = getindex.(coords, 2)

        push!(centroids, (mean(xs), mean(ys)))
    end

    return centroids
end

all_positions = []

for frame in frames
    mask = orange_mask(frame)
    centroids = find_orange_centroids(mask)
    push!(all_positions, centroids)
end

print(all_positions)


display(overlay)