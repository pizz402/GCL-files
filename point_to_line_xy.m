function d = point_to_line_xy(pt)

    pt2 = [sum(pt')/2;sum(pt')/2]';
    d = vecnorm(pt2' - pt');
    
end