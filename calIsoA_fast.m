function [index_IsoA] = calIsoA_fast(xy_deg,n,e,isPlot)
% INPUTS:    xy_deg - an 2 by N array with continuous data
%            n - size of the n by n grid over which the density is computed
%                n has to be a power of 2, otherwise n=2^ceil(log2(n));
%                the default value is 2^8;
% OUTPUT:   index_IsoA - the Isoline Area
%
% Zdravko Botev (2021). kernel density estimation
% (https://www.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation),
% MATLAB Central File Exchange. Retrieved November 20, 2021.
% 2021.11.20 by CH Du.

[~,density,X,Y]=kde2d(xy_deg',n);
p_ksd = density;
mesh_area = (max(max(X))-min(min(X)))*(max(max(Y))-min(min(Y)))/n/n;
p_edge_IsA1 = findby2(p_ksd,0.682,mesh_area,e);
p_edge_IsA2 = findby2(p_ksd,0.95,mesh_area,e);
index_IsoA = numel(find(p_ksd>p_edge_IsA1))*mesh_area;
X1 = X(1,:);
Y1 = Y(:,1);
%% plot
if isPlot==1
    pcolor(X1,Y1,p_ksd);
    shading interp; 
    colorbar; colormap(flipud(hot));alpha(0.7);hold on
    xlabel('Position X(degree)');ylabel('Position Y(degree)');title('ISOA Fast')
    scatter(xy_deg(1,:)',xy_deg(2,:)',10,'k','MarkerEdgeAlpha',0.3);hold on
    a=contour(X1,Y1,p_ksd, [p_edge_IsA1,p_edge_IsA2]);
    plot(a(1,2:a(2,1)+1),a(2,2:a(2,1)+1),'b','Linewidth',2);hold on;
    plot(a(1,a(2,1)+3:end),a(2,a(2,1)+3:end),'r','Linewidth',2);
end
end