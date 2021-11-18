function [index_BCEA1,index_BCEA2] = calBCEA(xy_deg,mesh_step,e,isPlot)
margin = 0.5;
% bcea_deg=[x0_deg;y0_deg];
x0_deg = xy_deg(1,:);
y0_deg = xy_deg(2,:);
index.rH = std(x0_deg);
index.rV = std(y0_deg);
pHV =corrcoef(x0_deg,y0_deg);
pHV=pHV(1,2);
index_BCEA1 = 2.291*pi*index.rH*index.rV*(1-pHV^2)^0.5;
% 求概率密度函数
%mesh_step=0.0005;
%e=0.0001;
X1=[min(x0_deg)-margin:mesh_step:max(x0_deg)+margin];Y1=[min(y0_deg)-margin:mesh_step:max(y0_deg)+margin];
[xL,yL]=meshgrid(X1,Y1);
[p_Gaus,mu_x,mu_y,sigma_x,sigma_y,rho] = PvalueXY(x0_deg,y0_deg,xL,yL);
% 计算BCEA区域椭圆位置
p_edge_BCEA1 = findby2(p_Gaus,0.682,mesh_step,e);
p_edge_BCEA2 = findby2(p_Gaus,0.95,mesh_step,e);
index_BCEA2 = numel(find(p_Gaus>p_edge_BCEA1))*mesh_step*mesh_step;
if isPlot==1
    pcolor(X1,Y1,p_Gaus);
    shading interp; 
    colorbar; colormap(flipud(hot));alpha(0.7);hold on
    xlabel('X');ylabel('Y');title('BCEA 方法')
    scatter(x0_deg',y0_deg',10,'k','MarkerEdgeAlpha',0.3);hold on
    plot_BCEA1 = find(abs(p_Gaus-p_edge_BCEA1)<50*e);
    raw=mod(plot_BCEA1,size(p_Gaus,1)); col=ceil(plot_BCEA1/size(p_Gaus,1));
    scatter(X1(col),Y1(raw),'r.');hold on
    plot_BCEA2 = find(abs(p_Gaus-p_edge_BCEA2)<30*e);
    raw=mod(plot_BCEA2,size(p_Gaus,1)); col=ceil(plot_BCEA2/size(p_Gaus,1));
    scatter(X1(col),Y1(raw),'b.');
end
end