function [index_BCEA_f,index_BCEA] = calBCEA(xy_deg,n,e,isPlot)

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
% calculate BCEA
mesh_area = prod((MAX_XY-MIN_XY)/n); 
p_edge_BCEA1 = findby2(p_Gaus,0.682,mesh_area,e);
p_edge_BCEA2 = findby2(p_Gaus,0.95,mesh_area,e);
index_BCEA = numel(find(p_Gaus>p_edge_BCEA1))*mesh_area;
if isPlot==1
    pcolor(X1,Y1,p_Gaus);
    shading interp; 
    colorbar; colormap(flipud(hot));alpha(0.7);hold on
    xlabel('X');ylabel('Y');title('BCEA Method')
    scatter(x_deg',y_deg',10,'k','MarkerEdgeAlpha',0.3);hold on
    a=contour(X1,Y1,p_Gaus, [p_edge_BCEA1,p_edge_BCEA2]);
    plot(a(1,2:a(2,1)+1),a(2,2:a(2,1)+1),'b','Linewidth',2);hold on;
    plot(a(1,a(2,1)+3:end),a(2,a(2,1)+3:end),'r','Linewidth',2);
end
end