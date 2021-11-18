%% 2021.11.8 ��FEMs�������ݴ���
% �˳���� BCEA �� ISA ���������˸��ֺͷ���
% ���Ե������ݷ���
% ��Ҫ�ʼǣ�
% 1�� BCEA��ISA �����϶��Ǹ����ܶ���68%�����,һ���Ƕ�Ԫ��̬�ֲ�,һ���ǻ��ں˺����ĸ����ܶȹ���
% 2�� BCEA����ֱ���ù�ʽ����,��ʽ���Ƶ��ڡ�Retinal Location of the Preferred Retinal Locus...��������ϸ���Ƶ�����
%     ��matlab��meshgrid��������,���ַ�ȷ���ٽ�����ܶ�,������ֵ���ֵõ��Ľ����һ���ġ��������е�  index.bcea2 �� index.BCEA_Cal��
% 3�� ��Ҫ������meshgrid���־�ϸ�Ⱥ�������mesh_step��e��,��ֵ����Ծ�ϸ��Ҫ�󲻸�,��0.01/0.0001����ǧ��֮һ����
%     �����Ҫ���ڶ��ַ��Ľ����жϣ�ͬʱ�ڻ�ͼʱ��Ҳ����68%������50(��ֵӰ�컭ͼЧ��)*e�ĺ�㹹�ɱ߽����
%     ��ͼ��Ҫ�ϸ߾��ȣ�������ʱ�䳤�����ж��ַ�����ռ��ʱ�䲻������Ҫ��meshgrid��ĺ˺��������ܶȹ���ռ�úܳ�ʱ��
%     �Ƽ�������BCEA��0.0005-0.0001
%�����Ӧ����·����
%load('F:\����\����SSVEP���ȶ��Է���\ʵ�鷶ʽ\2021 FEMs\DataSave\DCH_20211108\DCH_1108155215stimevent.mat')
%load('F:\����\����SSVEP���ȶ��Է���\ʵ�鷶ʽ\2021 FEMs\DataSave\DCH_20211108\DCH_1108155215GazeData.mat')
%%
clc;clear;
EYE_EEG_loadData;
% ���۶��źŷ�5��ͨ������GazeDataTrail��
GazeData_All=collected_gaze_data;
GazeDataLength = size(GazeData_All,1);
%TimeStampStar = GazeData_All(1, 1).SystemTimeStamp;
GazeDataTrail = zeros(5,GazeDataLength);%���� ���� ʱ��� 
for i=1:GazeDataLength
GazeDataTrail(1:2,i) = GazeData_All(i,1).LeftEye.GazePoint.OnDisplayArea';
GazeDataTrail(3:4,i) = GazeData_All(i,1).RightEye.GazePoint.OnDisplayArea';
GazeDataTrail(5,i) = GazeData_All(i, 1).SystemTimeStamp;%-TimeStampStar;
end
EEG_stamp = squeeze(stimevent.stamp); % EEG�ź��Դ�ʱ��

for trail = 1:9
    % ����ʱ��  �۶����ݷֶ�                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
    [~,trail_start] = min(abs(GazeDataTrail(5,:)*10^-6-EEG_stamp(trail,1)));
    [~,trail_end] = min(abs(GazeDataTrail(5,:)*10^-6-EEG_stamp(trail,2)));
    GazeData{trail} = GazeDataTrail(:,trail_start:trail_start+719);
    disp(trail_end-trail_start);
    EEGData{trail} = DATA_online{1, trail}{1, 1}(1:7200,:);
end
%% ���ݴ���
% ����
trail =4;
x1 = GazeData{trail}(1,:)*3840;y1 = GazeData{trail}(2,:)*2160;
x2 = GazeData{trail}(3,:)*3840;y2 = GazeData{trail}(4,:)*2160;% ת��������ֵ
x = mean([x1;x2]);y = mean([y1;y2]);
x0 = x-stimevent.StimLocations(1,trail);y0 = y-stimevent.StimLocations(2,trail);%ת�����������ֵ
x0_deg = atan(x0*0.16/2/600)/pi*180*2;
y0_deg = atan(y0*0.16/2/600)/pi*180*2; %ת������ԽǶ�
distance = (x0.^2+y0.^2).^0.5;  %ת�������ĵ����ֵ
figure;
plot([x;y]');hold on;
plot(repmat(stimevent.StimLocations(:,trail),1,720)','LineWidth',1);

% ����1 P1 P2
pixel_1 = round(tan(1/180*pi)*600/0.16);
pixel_2 = round(tan(2/180*pi)*600/0.16);%��ת��������ֵ
index.p1 = numel(find(distance<pixel_1))/720
index.p2 = numel(find(distance<pixel_2))/720
figure;
plot(distance);hold on;
plot(repmat([pixel_1;pixel_2],1,720)','LineWidth',1);

%% ����2 BCEA �����Ǹ�˹�ֲ���
% �����ֲ� ����ֵ����� 
a=0.5:0.001:0.9;
X2=chi2inv(a,2);%�����ֲ��Ŀ���ֵ��pֵ�Ĺ�ϵ
%figure;plot(X2);
%X2_bcea=chi2inv(0.682,2); %��֪ȡ68%�ĵ��ʱ��  X2 = 2.2914

% �������Ǵ���ģ�������ȡ68%���������ݵ㣬����ȡ��̬�ֲ��ϵ�68%�����
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
%% BCEA�����ⷨ
bcea_deg=[x0_deg;y0_deg];
index.rH = std(x0_deg);
index.rV = std(y0_deg);
pHV =corrcoef(x0_deg,y0_deg);
pHV=pHV(1,2);
index.bcea2 = 2.291*pi*index.rH*index.rV*(1-pHV^2)^0.5;
% ������ܶȺ���
mesh_step=0.0005;
e=0.0001;
X1=[0:mesh_step:1];Y1=[-0.5:mesh_step:1];
[xL,yL]=meshgrid(X1,Y1);
[p_Gaus,mu_x,mu_y,sigma_x,sigma_y,rho] = PvalueXY(x0_deg,y0_deg,xL,yL);
% ����BCEA������Բλ��
p_edge_BCEA1 = findby2(p_Gaus,0.682,mesh_step,e);
p_edge_BCEA2 = findby2(p_Gaus,0.95,mesh_step,e);
index.BCEA_Cal = numel(find(p_Gaus>p_edge_BCEA1))*mesh_step*mesh_step
%% ��ͼ_BCEA
%mesh(p_Gaus);shading interp; 
figure;
subplot(1,2,1);
pcolor(X1,Y1,p_Gaus);
shading interp; 
colorbar; colormap(flipud(hot));alpha(0.8);hold on
xlabel('X');ylabel('Y');title('BCEA ����')
scatter(x0_deg',y0_deg',10,'k','MarkerEdgeAlpha',0.3);hold on
plot_BCEA1 = find(abs(p_Gaus-p_edge_BCEA1)<50*e);
raw=mod(plot_BCEA1,size(p_Gaus,1)); col=ceil(plot_BCEA1/size(p_Gaus,1));
scatter(X1(col),Y1(raw),'r.');hold on
plot_BCEA2 = find(abs(p_Gaus-p_edge_BCEA2)<30*e);
raw=mod(plot_BCEA2,size(p_Gaus,1)); col=ceil(plot_BCEA2/size(p_Gaus,1));
scatter(X1(col),Y1(raw),'r.');hold on
%%  ����3 ��ֵ�����
% ���к��ܶȹ���
gridx1 = X1;
gridx2 = Y1;
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
%ksdensity([x_bcea_deg;y_bcea_deg]',xi);
[f,xr,bw] = ksdensity([x0_deg;y0_deg]',xi);
p_ksd = reshape(f,[numel(gridx2),numel(gridx1)]);
% ����Isoline����
p_edge_IsA1 = findby2(p_ksd,0.682,mesh_step,e);
p_edge_IsA2 = findby2(p_ksd,0.95,mesh_step,e);
index.IsA_Cal = numel(find(p_ksd>p_edge_IsA1))*mesh_step*mesh_step;
%% ��ͼ_Isoline Area
subplot(1,2,2);
pcolor(X1,Y1,p_ksd);
shading interp; 
colorbar; colormap(flipud(hot));alpha(0.7);hold on
xlabel('X');ylabel('Y');title('IsolineArea ����')
scatter(x0_deg',y0_deg',10,'k','MarkerEdgeAlpha',0.3);hold on
plot_IsA1 = find(abs(p_ksd-p_edge_IsA1)<50*e);
raw=mod(plot_IsA1,size(p_ksd,1)); col=ceil(plot_IsA1/size(p_ksd,1));
scatter(X1(col),Y1(raw),'r.');hold on
plot_IsA2 = find(abs(p_ksd-p_edge_IsA2)<30*e);
raw=mod(plot_IsA2,size(p_ksd,1)); col=ceil(plot_IsA2/size(p_ksd,1));
scatter(X1(col),Y1(raw),'r.');hold on


