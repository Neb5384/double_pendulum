#long functions I want to invoke in main here to make main less clustered

#functions needed to calculate speed at end of video
function calculate_angle(position1, position2)
    x = position2.x - position1.x
    y = - position2.y + position1.y

    return atan(x, y)
end

function calculate_w(pivot_pos1, pivot_pos2, pos1, pos2, fps, step)
    angle_1 = calculate_angle(pivot_pos1, pos1)
    angle_2 = calculate_angle(pivot_pos2, pos2)
    angle_diff = angle_2 - angle_1
    angle_diff = mod(angle_diff + pi, 2pi) - pi
    w = angle_diff * fps / step

    return w
end
