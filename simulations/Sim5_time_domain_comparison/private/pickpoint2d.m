function point = pickpoint2d(vx,vy,vertices,metric_param)
pointvx = find(abs(vertices(:,1) - vx) < metric_param);
pointvy = find(abs(vertices(:,2) - vy) < metric_param);
point    = intersect(pointvx,pointvy); 
if isempty(point)
    point    = intersect(pointvx,pointvy);
end
end