%% 2021.11.8 对FEMs进行数据处理
% 此程序对 BCEA 与 ISA 方法进行了复现和分析
% 仅对单次数据分析
% 重要笔记：
% 1、 BCEA和ISA 本质上都是概率密度上68%的面积,一个是二元正态分布,一个是基于核函数的概率密度估计
% 2、 BCEA可以直接用公式计算,公式的推导在【Retinal Location of the Preferred Retinal Locus...】中有详细的推导过程
%     用matlab的meshgrid划分网格,二分法确定临界概率密度,进行数值积分得到的结果是一样的【本程序中的  index.bcea2 与 index.BCEA_Cal】
% 3、 主要超参是meshgrid划分精细度和允许误差【mesh_step，e】,数值计算对精细度要求不高,在0.01/0.0001上有千分之一的误差。
%     误差主要用于二分法的结束判断，同时在画图时，也利用68%附近±50(该值影响画图效果)*e的红点构成边界的线
%     画图需要较高精度，但计算时间长。其中二分法程序占用时间不长，主要是meshgrid后的核函数概率密度估计占用很长时间
%     推荐参数：BCEA：0.0005-0.0001
%程序对应数据路径：
%load('F:\科研\基于SSVEP的稳定性分析\实验范式\2021 FEMs\DataSave\DCH_20211108\DCH_1108155215stimevent.mat')
%load('F:\科研\基于SSVEP的稳定性分析\实验范式\2021 FEMs\DataSave\DCH_20211108\DCH_1108155215GazeData.mat')
%%
clc;clear;
EYE_EEG_loadData;
% 将眼动信号分5个通道整理到GazeDataTrail中
GazeData_All=collected_gaze_data;
GazeDataLength = size(GazeData_All,1);
%TimeStampStar = GazeData_All(1, 1).SystemTimeStamp;
GazeDataTrail = zeros(5,GazeDataLength);%左眼 右眼 时间戳 
for i=1:GazeDataLength
GazeDataTrail(1:2,i) = GazeData_All(i,1).LeftEye.GazePoint.OnDisplayArea';
GazeDataTrail(3:4,i) = GazeData_All(i,1).RightEye.GazePoint.OnDisplayArea';
GazeDataTrail(5,i) = GazeData_All(i, 1).SystemTimeStamp;%-TimeStampStar;
end
EEG_stamp = squeeze(stimevent.stamp); % EEG信号自带时标

for trail = 1:9
    % 根据时标  眼动数据分段                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
    [~,trail_start] = min(abs(GazeDataTrail(5,:)*10^-6-EEG_stamp(trail,1)));
    [~,trail_end] = min(abs(GazeDataTrail(5,:)*10^-6-EEG_stamp(trail,2)));
    GazeData{trail} = GazeDataTrail(:,trail_start:trail_start+719);
    disp(trail_end-trail_start);
    EEGData{trail} = DATA_online{1, trail}{1, 1}(1:7200,:);
end
%% 数据处理
% 整体
trail =4;
x1 = GazeData{trail}(1,:)*3840;y1 = GazeData{trail}(2,:)*2160;
x2 = GazeData{trail}(3,:)*3840;y2 = GazeData{trail}(4,:)*2160;% 转换成像素值
x = mean([x1;x2]);y = mean([y1;y2]);
x0 = x-stimevent.StimLocations(1,trail);y0 = y-stimevent.StimLocations(2,trail);%转换成相对像素值
x0_deg = atan(x0*0.16/2/600)/pi*180*2;
y0_deg = atan(y0*0.16/2/600)/pi*180*2; %转换成相对角度
distance = (x0.^2+y0.^2).^0.5;  %转换成中心点距离值
figure;
plot([x;y]');hold on;
plot(repmat(stimevent.StimLocations(:,trail),1,720)','LineWidth',1);

% 方法1 P1 P2
pixel_1 = round(tan(1/180*pi)*600/0.16);
pixel_2 = round(tan(2/180*pi)*600/0.16);%°转换成像素值
index.p1 = numel(find(distance<pixel_1))/720
index.p2 = numel(find(distance<pixel_2))/720
figure;
plot(distance);hold on;
plot(repmat([pixel_1;pixel_2],1,720)','LineWidth',1);

%% 方法2 BCEA 假设是高斯分布的
% 卡方分布 卡方值的求解 
a=0.5:0.001:0.9;
X2=chi2inv(a,2);%卡方分布的卡方值与p值的关系
%figure;plot(X2);
%X2_bcea=chi2inv(0.682,2); %可知取68%的点的时候  X2 = 2.2914

% 这种求法是错误的，并不是取68%的中心数据点，而是取正态分布上的68%的面积
% [~,i]=sort(distance);
% i=i(1:round(720*0.68));                                      
% x_bcea = x0(i);y_bcea = y0(i);
% x_bcea_de = atan(x_bcea*0.16/2/600)/pi*180*2;
% y_bcea_de = atan(y_bcea*0.16/2/600)/pi*180*2;
% index.rH = std(x_bcea_de);
% index.rV = std(y_bcea_de);
% pHV =corrcoef(x_bcea_de,y_bcea_de);
% pHV=pHV(1,2);
% index.bcea1 = 2.28*pi*index.rH*index.rV*(1-pHV^2)^0.5;
%% BCEA更正解法
bcea_deg=[x0_deg;y0_deg];
index.rH = std(x0_deg);
index.rV = std(y0_deg);
pHV =corrcoef(x0_deg,y0_deg);
pHV=pHV(1,2);
index.bcea2 = 2.291*pi*index.rH*index.rV*(1-pHV^2)^0.5;
% 求概率密度函数
mesh_step=0.0005;
e=0.0001;
X1=[0:mesh_step:1];Y1=[-0.5:mesh_step:1];
[xL,yL]=meshgrid(X1,Y1);
[p_Gaus,mu_x,mu_y,sigma_x,sigma_y,rho] = PvalueXY(x0_deg,y0_deg,xL,yL);
% 计算BCEA区域椭圆位置
p_edge_BCEA1 = findby2(p_Gaus,0.682,mesh_step,e);
p_edge_BCEA2 = findby2(p_Gaus,0.95,mesh_step,e);
index.BCEA_Cal = numel(find(p_Gaus>p_edge_BCEA1))*mesh_step*mesh_step
%% 画图_BCEA
%mesh(p_Gaus);shading interp; 
figure;
subplot(1,2,1);
pcolor(X1,Y1,p_Gaus);
shading interp; 
colorbar; colormap(flipud(hot));alpha(0.8);hold on
xlabel('X');ylabel('Y');title('BCEA 方法')
scatter(x0_deg',y0_deg',10,'k','MarkerEdgeAlpha',0.3);hold on
plot_BCEA1 = find(abs(p_Gaus-p_edge_BCEA1)<50*e);
raw=mod(plot_BCEA1,size(p_Gaus,1)); col=ceil(plot_BCEA1/size(p_Gaus,1));
scatter(X1(col),Y1(raw),'r.');hold on
plot_BCEA2 = find(abs(p_Gaus-p_edge_BCEA2)<30*e);
raw=mod(plot_BCEA2,size(p_Gaus,1)); col=ceil(plot_BCEA2/size(p_Gaus,1));
scatter(X1(col),Y1(raw),'r.');hold on
%%  方法3 等值线面积
% 进行核密度估计
gridx1 = X1;
gridx2 = Y1;
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
%ksdensity([x_bcea_deg;y_bcea_deg]',xi);
[f,xr,bw] = ksdensity([x0_deg;y0_deg]',xi);
p_ksd = reshape(f,[numel(gridx2),numel(gridx1)]);
% 计算Isoline区域
p_edge_IsA1 = findby2(p_ksd,0.682,mesh_step,e);
p_edge_IsA2 = findby2(p_ksd,0.95,mesh_step,e);
index.IsA_Cal = numel(find(p_ksd>p_edge_IsA1))*mesh_step*mesh_step;
%% 画图_Isoline Area
subplot(1,2,2);
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
scatter(X1(col),Y1(raw),'r.');hold on


