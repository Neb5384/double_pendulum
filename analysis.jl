# using Pkg

# Pkg.add("ImageShow")
# Pkg.add("VideoIO")
# Pkg.add("Images")


using VideoIO, Images,ImageShow, ImageMorphology

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
        (10 <= c.h <= 35) &&   # orange hue range
        (c.s >= 0.5) &&      # saturation
        (0.4 <= c.v <= 0.7)         # brightness
    , img_hsv)

    return mask
end



test_frame = frames[60]

mask = orange_mask(test_frame)
mask_int = Int.(mask)
mask_clean = area_opening(mask_int, min_area = 200)

overlay = colorview(RGB,
    channelview(RGB.(test_frame)) .* reshape(mask_clean, 1, size(mask)...)
)
display(overlay)

labels = label_components(mask_clean)
print(maximum(labels)) 

