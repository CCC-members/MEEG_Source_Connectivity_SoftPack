function point = pickpoint(vx,vy,vz,vertices,metric_param)
pointvx = find(abs(vertices(:,1) - vx) < metric_param);
pointvy = find(abs(vertices(:,2) - vy) < metric_param);
pointvz = find(abs(vertices(:,3) - vz) < metric_param);
point    = intersect(intersect(pointvx,pointvy),pointvz); 
if isempty(point)
    point    = intersect(pointvx,pointvy);
    if isempty(point)
        point    = intersect(pointvx,pointvz);
        if isempty(point)
            point    = intersect(pointvy,pointvz);
        end
    end
end
end