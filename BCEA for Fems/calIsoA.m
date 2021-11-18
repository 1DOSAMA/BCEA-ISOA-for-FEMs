function [index_IsoA] = calIsoA(xy_deg,mesh_step,e,isPlot)
% bcea_deg=[x0_deg;y0_deg];
margin = 0.5;
x0_deg = xy_deg(1,:);
y0_deg = xy_deg(2,:);
X1=[min(x0_deg)-margin:mesh_step:max(x0_deg)+margin];
Y1=[min(y0_deg)-margin:mesh_step:max(y0_deg)+margin];
[x1,x2] = meshgrid(X1, Y1);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
% ksdensity([x_bcea_deg;y_bcea_deg]',xi);
[f,xr,bw] = ksdensity([x0_deg;y0_deg]',xi);
p_ksd = reshape(f,[numel(Y1),numel(X1)]);
% 计算Isoline区域
p_edge_IsA1 = findby2(p_ksd,0.682,mesh_step,e);
p_edge_IsA2 = findby2(p_ksd,0.95,mesh_step,e);
index_IsoA = numel(find(p_ksd>p_edge_IsA1))*mesh_step*mesh_step;
%% 画图
if isPlot==1
    pcolor(X1,Y1,p_ksd);
    shading interp; 
    colorbar; colormap(flipud(hot));alpha(0.7);hold on
    xlabel('X');ylabel('Y');title('IsolineArea 方法')
    scatter(x0_deg',y0_deg',10,'k','MarkerEdgeAlpha',0.3);hold on
    plot_IsA1 = find(abs(p_ksd-p_edge_IsA1)<50*e);
    raw=mod(plot_IsA1,size(p_ksd,1)); col=ceil(plot_IsA1/size(p_ksd,1));
    scatter(X1(col),Y1(raw),'r.');hold on
    plot_IsA2 = find(abs(p_ksd-p_edge_IsA2)<30*e);
    raw=mod(plot_IsA2,size(p_ksd,1)); col=ceil(plot_IsA2/size(p_ksd,1));
    scatter(X1(col),Y1(raw),'b.');
end

end