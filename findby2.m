function p_edge = findby2(p,edge,mesh_area,e)
% findby2 find the probability density boundary by binary search
%
% INPUTS:    p - n by n matrix containing the density values over the n by n grid
%            edge - percent of the encompassed highest density points,
%            usually at 0.682 or 0.95
%            mesh_area - area of every grid(probability density)
%            e - allowable error to find the probability density boundary in
%               binary search
% OUTPUT:   p_edge - probability value at the boundary(edge) of the density
%
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