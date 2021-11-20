function p_edge = findby2(p,edge,mesh_area,e)
% find the critical value by binary search
% 2021.11.20 by CH Du.

p_edge_m = max(max(p));
p_edge_0 = min(min(p));
r = 0;
r_max = 200;
while r < r_max
    p_edge = (p_edge_m+p_edge_0)/2;
    p_sum = sum(p(p>p_edge)*mesh_area);
    %disp(p_sum-edge);
    if abs(p_sum-edge)<e
        break
    elseif p_sum<edge
        p_edge_m = p_edge;
    elseif p_sum>edge
        p_edge_0 = p_edge;
    end
    r = r+1;
end
if r == r_max
    disp('binary search error rate =');
    disp((p_sum-edge)/e);
end
end