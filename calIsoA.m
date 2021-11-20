function [index_IsoA,density] = calIsoA(xy_deg,n,e,isPlot)
% calIsoA  Compute IsolineArea of fixational eye movements 
%
% INPUTS:    xy_deg - an 2 by N array with continuous data
%            n - size of the n by n grid over which the density is computed
%                n has to be a power of 2, otherwise n=2^ceil(log2(n))
%            e - allowable error to find the probability density boundary in
%               binary search
% OUTPUT:   index_IsoA - the Isoline Area
%           density - an n by n matrix containing the density values over the n by n grid
%
% 2021.11.20 by CH Du.

%The dafault limits of the bounding box over which the density is computed are computed as: 
%(or you can define it by yourself)
MAX=max(xy_deg,[],2); MIN=min(xy_deg,[],2); Range=MAX-MIN;
MAX_XY=MAX+Range/2; MIN_XY=MIN-Range/2;
x_deg = xy_deg(1,:);
y_deg = xy_deg(2,:);
X1 = linspace(MIN_XY(1),MAX_XY(1),n);
Y1 = linspace(MIN_XY(2),MAX_XY(2),n);
[x1,x2] = meshgrid(X1, Y1);
xi = [x1(:) x2(:)];
[f,xr,bw] = ksdensity(xy_deg',xi);
p_ksd = reshape(f,[n,n]);
density = p_ksd;
% Isoline Area
mesh_area = prod((MAX_XY-MIN_XY)/n); 
p_edge_IsA1 = findby2(p_ksd,0.682,mesh_area,e);
p_edge_IsA2 = findby2(p_ksd,0.95,mesh_area,e);
index_IsoA = numel(find(p_ksd>p_edge_IsA1))*mesh_area;
%% plot
if isPlot==1
    pcolor(X1,Y1,p_ksd);
    shading interp; 
    colorbar; colormap(flipud(hot));alpha(0.7);hold on
    xlabel('Position X(degree)');ylabel('Position Y(degree)');title('ISOA')
    scatter(x_deg',y_deg',10,'k','MarkerEdgeAlpha',0.3);hold on
    a=contour(X1,Y1,p_ksd, [p_edge_IsA1,p_edge_IsA2]);
    plot(a(1,2:a(2,1)+1),a(2,2:a(2,1)+1),'b','Linewidth',2);hold on;
    plot(a(1,a(2,1)+3:end),a(2,a(2,1)+3:end),'r','Linewidth',2);
end

end