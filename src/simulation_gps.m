% 函数功能：考察gps定位误差椭球分布图 V
tic
clc;clear;
format long g
addpath(genpath("../../coordinateTransformation"));
satECEFPositionList = ...
[1.615202609797593E7, 3347999.121093184, 2.0654217845651966E7;...
 -9907042.824254453, 2.287805023963695E7, 8874684.143130684;...
 -1.9061270789226927E7, 3058189.769791643, 1.7922650491999816E7;...
 -1.542310695696636E7, 2.1054512181760143E7, -4117013.5301650744;...
 1.8525465084249824E7, 1.0357349901848707E7, 1.696803855309986E7;...
 -1.6719978820409523E7, 2.0520602009612914E7, -1927801.3735653749;...
 -1.829436750996808E7, -8079275.716146787, 1.742108925752629E7;...
-1.9302924837284807E7, 1.6167012674498076E7, 7735366.371604006;...
6131269.763786005, 2.5163298085870665E7, 4552437.166840457;...
-2093985.8407936988, 1.5185238960106703E7, 2.178821106732648E7];
satECEFVelocityList = ...
[931.638082885227, 2496.0533034633218, -1086.9526899935393;...
-884.1355175695837, 730.7897865848088, -2933.0537381917597;...
1439.8959741356891, -1836.5020286658437, 1834.5112596649858;...
-67.03857050317686, -638.6669534987794, -3080.36321658192;...
763.1332678124891, 1855.6928724616496, -2012.0498613051257;...
-79.89910497697765, -359.3605968666734, -3170.8818885708333;...
-821.4624063204889, -2107.2121764003555, -1775.7184609880733;...
398.06578559418415, -987.4290752337673, 2985.452930906086;...
-521.32077374667, -439.5794612784424, 3123.73711669872;...
-2720.703582581369, -358.64747298162075, -28.14697920702789;...
];
userECEFPosition = [-2195015.633774167, 4389738.432204695, 4059886.720907325];
userLLHPosition  = ecef2llh(userECEFPosition);
B0               = 7.3;
B1               = 0.03;
prError          = 4.5;
pdError          = 1;
c                = 299792458;
f0               = 1575.42E6;
iterrT           = 1000;
N                = size(satECEFPositionList,1);
withLEOFlag      = 1;
withLEOPrFlag    = 1;
withLEOPdFlag    = 1;
weightFlag       = 1;
scale            = 4;
lineWidth        = 2;
maxHeadSize      = 1;

 
fakeRangeList    = [];
fakeDopplerList  = [];
for i=1:N
    realRange        = norm(userECEFPosition-satECEFPositionList(i,:));
    fakeRange        = realRange + B0;
    fakeRangeList    = [fakeRangeList, fakeRange];
    fakedoppler      = (satECEFVelocityList(i,:)*(userECEFPosition - satECEFPositionList(i,:))'/realRange + B1)*(f0/c);       
    fakeDopplerList  = [fakeDopplerList, fakedoppler];            
end
% 星空图

[azimuthList, elevationList] = calAziAndEle(userLLHPosition, satECEFPositionList);
starsky(azimuthList, elevationList, '*');

errorByIterr      = [];
errorHoriByIterr  = [];
errorEastByIterr  = [];
errorNorthByIterr = [];
errorVetiByIterr  = [];
errorXByIterr     = [];
errorYByIterr     = [];
errorZByIterr     = [];
PDOPByIterr       = [];
HDOPByIterr       = [];
VDOPByIterr       = [];
B0ByIterr         = [];
for iterr = 1:iterrT
    prList = [];
    pdList = [];
    for i = 1:N
        pr     = fakeRangeList(i) + randn*prError;
        pd     = fakeDopplerList(i) + rand*pdError;
        prList = [prList, pr];
        pdList = [pdList, pd];
    end
    [x0List, y0List, z0List, B00List, B10List,...
    PDOP, HDOP, VDOP, iterTime] = ...
    LSM60(f0,...
    prList, pdList, satECEFPositionList, satECEFVelocityList,... % 卫星观测量
    [], [], [],...                 % VOR/DME观测量
    0,...                                                         % 高度观测量
    prError, pdError, 0, 0, 0,...% 观测误差
    withLEOFlag, withLEOPrFlag, withLEOPdFlag,... % LEO pr pd
    0, 0, 0, weightFlag,... % vor dme height    加权
    [0,0,0], 0, 0, 1);  % startPos b0 b1 clcOrder

    %%%%%%%%%%%%%%%%%数据统计%%%%%%%%%%%%%%%%%%
    calPos            = [x0List(end),y0List(end),z0List(end)];
    error             = norm(calPos-userECEFPosition);
    [errorEast, errorNorth, errorHorizontal, errorVertical]=...
    errorECEF2ENU(calPos, userECEFPosition);
    errorByIterr      = [errorByIterr,      error];
    errorHoriByIterr  = [errorHoriByIterr,  errorHorizontal];
    errorEastByIterr  = [errorEastByIterr,  errorEast];
    errorNorthByIterr = [errorNorthByIterr, errorNorth];
    errorVetiByIterr  = [errorVetiByIterr,  errorVertical];
    errorXByIterr     = [errorXByIterr,     calPos(1)-userECEFPosition(1)];
    errorYByIterr     = [errorYByIterr,     calPos(2)-userECEFPosition(2)];
    errorZByIterr     = [errorZByIterr,     calPos(3)-userECEFPosition(3)];
    B0ByIterr         = [B0ByIterr,         B00List(end)];
    PDOPByIterr       = [PDOPByIterr,       PDOP];
    HDOPByIterr       = [HDOPByIterr,       HDOP];
    VDOPByIterr       = [VDOPByIterr,       VDOP];
    
    if(rem(iterr,50) == 0)
        disp(['第', num2str(iterr), '次仿真结束']);
    end
end
%%%%%%%%    SEP     %%%%%%%%
[Vg_s, DOP1, DOP2, DOP3] = calSEP(...
f0,...
userECEFPosition,...
satECEFPositionList, satECEFVelocityList,[],...                                          
prError, pdError, 0, 0, 0,...% 观测误差
withLEOFlag, withLEOPrFlag, withLEOPdFlag,...
0, 0, 0,...
1, 1); % withclcshif / draft
disp('Vg_s:');
matrixSout(Vg_s)
disp(['DOP1 =           ',num2str(DOP1)]);
disp(['DOP2 =           ',num2str(DOP2)]);
disp(['DOP3 =           ',num2str(DOP3)]);
if(withLEOPrFlag == 1 && withLEOPdFlag == 0)
    v1 = Vg_s(:,1)*DOP1*prError;
    v2 = Vg_s(:,2)*DOP2*prError;
    v3 = Vg_s(:,3)*DOP3*prError;
elseif(withLEOPrFlag == 0 && withLEOPdFlag == 1)
    v1 = Vg_s(:,1)*DOP1*pdError;
    v2 = Vg_s(:,2)*DOP2*pdError;
    v3 = Vg_s(:,3)*DOP3*pdError;
elseif(withLEOPrFlag == 1 && withLEOPdFlag == 1)
    v1 = Vg_s(:,1)*DOP1;
    v2 = Vg_s(:,2)*DOP2;
    v3 = Vg_s(:,3)*DOP3;
else
end
lon0rad = userLLHPosition(1)*pi/180;
lat0rad = userLLHPosition(2)*pi/180;
R = [-sin(lon0rad), cos(lon0rad), 0;
    -sin(lat0rad)*cos(lon0rad), -sin(lat0rad)*sin(lon0rad), cos(lat0rad);
    cos(lat0rad)*cos(lon0rad), cos(lat0rad)*sin(lon0rad), sin(lat0rad)];
v1_ENU = R * v1;
v2_ENU = R * v2;
v3_ENU = R * v3;
disp('三轴转到东北天方向上的结果')
matrixSout(R*Vg_s);
vec = v1_ENU + v2_ENU + v3_ENU;
disp(['|vec| = ', num2str(norm(vec))]);
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    CRLB     %%%%%%%%
CRLBMatrix = calCRLB(...
f0,...
userECEFPosition,...
satECEFPositionList, satECEFVelocityList,[],...                                                      
prError, pdError, 0, 0, 0,...% 观测误差
withLEOFlag, withLEOPrFlag, withLEOPdFlag,...
0, 0, 0,...
1,1);

posCRLB_2 = 2*sqrt(CRLBMatrix(1,1)+CRLBMatrix(2,2)+CRLBMatrix(3,3));
disp(['posCRLB_2 =      ', num2str(posCRLB_2)]);

errorX_var = var(errorXByIterr);
errorY_var = var(errorYByIterr);
errorZ_var = var(errorZByIterr);
posCRLB = sqrt(CRLBMatrix(1,1)+CRLBMatrix(2,2)+CRLBMatrix(3,3));
disp(['error_std =      ', num2str(sqrt(errorX_var+errorY_var+errorZ_var))]);
disp(['posCRLB =        ', num2str(posCRLB)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 误差分布（ECEF + ENU + EN）%%%%%%%%
range1   = 25;
range2   = 25;
fontSize = 20;
numSize  = 20;
figure(2)
    subplot(1,3,1) 
    scatter3(errorEastByIterr, errorNorthByIterr, errorVetiByIterr, numSize, 'filled', 'r');
    hold on;
    quiver3(0,0,0,scale*v1_ENU(1,1),scale*v1_ENU(2,1),scale*v1_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
    hold on;
    quiver3(0,0,0,scale*v2_ENU(1,1),scale*v2_ENU(2,1),scale*v2_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
    hold on; 
    quiver3(0,0,0,scale*v3_ENU(1,1),scale*v3_ENU(2,1),scale*v3_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
    hold on;
    quiver3(0,0,0,-scale*v1_ENU(1,1),-scale*v1_ENU(2,1),-scale*v1_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
    hold on;
    quiver3(0,0,0,-scale*v2_ENU(1,1),-scale*v2_ENU(2,1),-scale*v2_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
    hold on; 
    quiver3(0,0,0,-scale*v3_ENU(1,1),-scale*v3_ENU(2,1),-scale*v3_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
    xlim([-range1 range1]);
    ylim([-range1 range1]);
    zlim([-range1 range1]);

    xlabel('东向误差/(m)', FontSize=fontSize);
    ylabel('北向误差/(m)', FontSize=fontSize);
    zlabel('天向误差/(m)', FontSize=fontSize);
    title('定位误差', FontSize=fontSize);
    
    subplot(1,3,2) 
    scatter(errorEastByIterr, errorNorthByIterr, numSize, 'filled', 'r');
    xlabel('东向误差/(m)', FontSize=fontSize);
    ylabel('北向误差/(m)', FontSize=fontSize);
    title('水平定位误差', FontSize=fontSize);
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-range1 range1])
    ylim([-range1 range1])
    title('水平误差');
    
    subplot(1,3,3)
    scatter(errorNorthByIterr, errorVetiByIterr, numSize, 'filled', 'r');
    xlabel('北向误差/(m)', FontSize=fontSize);
    ylabel('天向误差/(m)', FontSize=fontSize);
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-range2 range2])
    ylim([-range2 range2])
    title('北向天向误差', FontSize=fontSize);
%     subplot(1,3,3);
%     scatter3(errorXByIterr, errorYByIterr, errorZByIterr, 'filled', 'r');
%     hold on;
%     quiver3(0,0,0,scale*v1(1,1),scale*v1(2,1),scale*v1(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
%     hold on;
%     quiver3(0,0,0,scale*v2(1,1),scale*v2(2,1),scale*v2(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
%     hold on; 
%     quiver3(0,0,0,scale*v3(1,1),scale*v3(2,1),scale*v3(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
%     hold on;
%     quiver3(0,0,0,-scale*v1(1,1),-scale*v1(2,1),-scale*v1(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
%     hold on;
%     quiver3(0,0,0,-scale*v2(1,1),-scale*v2(2,1),-scale*v2(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
%     hold on; 
%     quiver3(0,0,0,-scale*v3(1,1),-scale*v3(2,1),-scale*v3(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
%     
%     xlabel('X误差/(m)');
%     ylabel('Y误差/(m)');
%     zlabel('Z误差/(m)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc