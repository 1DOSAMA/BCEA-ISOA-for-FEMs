function p_edge = findby2(p,edge,step,e)
p_edge_m=max(max(p));
p_edge_0=min(min(p));
r = 0;
while r < 500
    p_edge=(p_edge_m+p_edge_0)/2;
    temp_index = find(p>p_edge);
    p_sum = sum(p(temp_index)*step^2);
    if abs(p_sum-edge)<e
        break
    elseif p_sum<edge
        p_edge_m = p_edge;
    elseif p_sum>edge
        p_edge_0 = p_edge;
    end
    r = r+1;
    %disp(r);
end
if r == 500
    disp('二分法到达最大值，当前误差');
    disp((p_sum-edge)/e);
end
end