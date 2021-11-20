function [index_BCEA_f,index_BCEA,density] = calBCEA(xy_deg,n,e,isPlot)
% calIsoA_fast  Compute bivariate contour ellipse Area of fixational eye movements 
%
% INPUTS:   xy_deg - an 2 by N array with continuous data
%           n - size of the n by n grid over which the density is computed
%               grid by n*n, but infact you can also compute density by n*m, make no difference
%           e - allowable error to find the probability density boundary in
%               binary search
% OUTPUT:   index_BCEA_f - calculate BCEA by formula
%           index_BCEA - calculate BCEA by probability density 
%           density - an n by n matrix containing the density values over the n by n grid
% Reference：
% Zdravko Botev (2021). kernel density estimation
% (https://www.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation),
% MATLAB Central File Exchange. Retrieved November 20, 2021.

% 2021.11.20 by CH Du.

%The dafault limits of the bounding box over which the density is computed are computed as: 
%(or you can define it by yourself)
MAX=max(xy_deg,[],2); MIN=min(xy_deg,[],2); Range=MAX-MIN;
MAX_XY=MAX+Range/2; MIN_XY=MIN-Range/2;
x_deg = xy_deg(1,:);
y_deg = xy_deg(2,:);

% calculate BCEA by formula
index.rH = std(x_deg);
index.rV = std(y_deg);
pHV =corrcoef(x_deg,y_deg);
pHV=pHV(1,2);
index_BCEA_f = chi2inv(0.682,2)*pi*index.rH*index.rV*(1-pHV^2)^0.5;
% calculate BCEA by probability density (same with ISOA)
X1 = linspace(MIN_XY(1),MAX_XY(1),n);
Y1 = linspace(MIN_XY(2),MAX_XY(2),n);
[x1,y1] = meshgrid(X1, Y1);
[p_Gaus,mu_x,mu_y,sigma_x,sigma_y,rho] = PvalueXY(x_deg,y_deg,x1,y1);
density = p_Gaus;
mesh_area = prod((MAX_XY-MIN_XY)/n); 
p_edge_BCEA1 = findby2(p_Gaus,0.682,mesh_area,e);
p_edge_BCEA2 = findby2(p_Gaus,0.95,mesh_area,e);
index_BCEA = numel(find(p_Gaus>p_edge_BCEA1))*mesh_area;
if isPlot==1
    pcolor(X1,Y1,p_Gaus);
    shading interp; 
    colorbar; colormap(flipud(hot));alpha(0.7);hold on
    xlabel('Position X(degree)');ylabel('Position Y(degree)');title('BCEA Method')
    scatter(x_deg',y_deg',10,'k','MarkerEdgeAlpha',0.3);hold on
    a=contour(X1,Y1,p_Gaus, [p_edge_BCEA1,p_edge_BCEA2]);
    plot(a(1,2:a(2,1)+1),a(2,2:a(2,1)+1),'b','Linewidth',2);hold on;
    plot(a(1,a(2,1)+3:end),a(2,a(2,1)+3:end),'r','Linewidth',2);
end
end