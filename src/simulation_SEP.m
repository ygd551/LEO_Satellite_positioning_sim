tic
clc;clear;
format long g

addpath(genpath("../../coordinateTransformation"));
paraAll = simulationPara();
%%%%%%%%%%%%%%%%%%%参数设置%%%%%%%%%%%%%%%%%%%%%%
% 用户位置
userLLHPosition  = paraAll.userLLHPosition;
userECEFPosition = llh2ecef(userLLHPosition);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%取出观测量%%%%%%%%%%%%%%%%%%%%
%%%%% LEO观测量%%%%%%
disp(paraAll.obsDataFilePath)
[ttList, satECEFPositionList, satECEFVelocityList,...
      rangeList, dopplerList] = getLEOObsData(paraAll.obsDataFilePath);
Ns = length(rangeList);
disp(['多历元卫星观测量的个数为：',num2str(Ns)]);

% [azimuthList, elevationList] = calAziAndEle(userLLHPosition, satECEFPositionList);
% starsky(azimuthList, elevationList, '.');
% %%%%% 所有理论上可见的VOR/DME得到的观测量%%%%%%
% VORDMEDataAll               = getVORDMEData(userECEFPosition);
% obsVORDMEECEFPostionListAll = VORDMEDataAll.obsVORDMEECEFPostionList;
% realDisListAll              = VORDMEDataAll.realDisList;
% realAngleListAll            = VORDMEDataAll.realAngleList;
% Nvd                         = length(realDisListAll);
% disp(['用户可见陆基台站个数：',num2str(Nvd)]);
% % 给出所有台站相对于用户的位置图
% enuList = [];
% for i = 1 : length(obsVORDMEECEFPostionListAll)
%     enu = ecef2enu(obsVORDMEECEFPostionListAll(i,:), userECEFPosition, userLLHPosition);
%     enuList = [enuList; enu];
% end
% index1 = 2;
% index2 = 5;
% obsVORDMEECEFPostionList = [obsVORDMEECEFPostionListAll(index1, :);
%     obsVORDMEECEFPostionListAll(index2, :)];
% realDisList = [realDisListAll(index1), realDisListAll(index2)];
% realAngleList = [realAngleListAll(index1), realAngleListAll(index2)];
% disp(['选取DME个数', num2str(length(realDisList))]);
% disp("==============================================")
%%%%%%%%%%%%%%%%%%得到初始值%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(paraAll.positionStartWay == 0)
    userECEFPositionStart = getInitialPos(rangeList, satECEFPositionList);
end
if(paraAll.positionStartWay == 1)
    userECEFPositionStart = [userECEFPosition(1) + paraAll.errorxyz(1), userECEFPosition(2) + paraAll.errorxyz(2), userECEFPosition(3) + paraAll.errorxyz(3)];
end
if(paraAll.positionStartWay == 2)
    userECEFPositionStart = [0, 0, 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errorByIterr       = [];
errorHoriByIterr   = [];
errorEastByIterr   = [];
errorEastByIterrJ  = [];
errorNorthByIterr  = [];
errorNorthByIterrJ = [];
errorVetiByIterr   = [];
errorVetiByIterrJ  = [];
errorXByIterr      = [];
errorYByIterr      = [];
errorZByIterr      = [];
PDOPByIterr        = [];
HDOPByIterr        = [];
VDOPByIterr        = [];
B0ByIterr          = [];
errorXIteraMatrix  = zeros(paraAll.iterrT, 10);
errorYIteraMatrix  = zeros(paraAll.iterrT, 10);
errorZIteraMatrix  = zeros(paraAll.iterrT, 10);
for iterr = 1:paraAll.iterrT
    %%%%%%%%%%%%%%%% 加观测误差%%%%%%%%%%%%%%%%
    % input:rangeList dopplerList
    prSatList = [];
    pdSatList = [];
    
    errorPrList = [];
    for ns = 1:Ns
        errorPr   = paraAll.pseudoRangeError*randn;
%         errorPrList = [errorPrList, errorPr];
        prSat     = rangeList(ns) + errorPr + paraAll.satError*randn + paraAll.ionoError*randn;
        prSatList = [prSatList, prSat];
        errorPd   = paraAll.dopplerError*randn;
        pdSat     = dopplerList(ns) + errorPd;
        pdSatList = [pdSatList, pdSat];
    end
    % output:prSatList pdSatList
    
    % input:realDisList realAngleList(deg)
%     prDMEList = [];
%     paVORList = [];
%     for nvd = 1:2
%         prDME     = realDisList(nvd) + paraAll.DMEError*randn;
%         prDMEList = [prDMEList, prDME];
%         paVOR     = realAngleList(nvd) * (1 + paraAll.VORError*randn) * pi/180;
%         paVORList = [paVORList, paVOR];
%     end
    % output:prDMEList paVORList(rad)
    
    % input:userLLHPosition(3)
    height = userLLHPosition(3) + (paraAll.hError + paraAll.hErrorSys)*randn;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     prSatList
%     pdSatList
%     prDMEList
%     paVORList
%     height
    
    %%%%%%%%%%%%%%%%LSM%%%%%%%%%%%%%%%%%%%%%
    [x0List, y0List, z0List, B00List, B10List,...
    PDOP, HDOP, VDOP, iterTime] = ...
    LSM60(paraAll.f0,...
    prSatList, pdSatList, satECEFPositionList, satECEFVelocityList,... % 卫星观测量
    [], [], [],...                 % VOR/DME观测量
    height,...                                                         % 高度观测量
    paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
    paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
    paraAll.withVORFlag, paraAll.withDMEFlag, paraAll.withHeightFlag, paraAll.weightFlag,...
    userECEFPositionStart, 0, 0, paraAll.clcOrder); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     errorList = [];
%     for ii = 1:length(x0List)
%         errorList = [errorList,norm(userECEFPosition - [x0List(ii),y0List(ii),z0List(ii)])];
%     end
%     errorList
    %%%%%%%%%%%%%%%%%数据统计%%%%%%%%%%%%%%%%%%
    calPos            = [x0List(end),y0List(end),z0List(end)];
    error             = norm(calPos-userECEFPosition);
    [errorEast, errorNorth, errorHorizontal, errorVertical]=...
    errorECEF2ENU(calPos, userECEFPosition);
    errorByIterr       = [errorByIterr,      error];
    errorHoriByIterr   = [errorHoriByIterr,  errorHorizontal];
    errorEastByIterr   = [errorEastByIterr,  errorEast];
    errorEastByIterrJ  = [errorEastByIterrJ,  abs(errorEast)];
    errorNorthByIterr  = [errorNorthByIterr, errorNorth];
    errorNorthByIterrJ = [errorNorthByIterrJ, abs(errorNorth)];
    errorVetiByIterr   = [errorVetiByIterr,  errorVertical];
    errorVetiByIterrJ  = [errorVetiByIterrJ,  abs(errorVertical)];
    errorXByIterr      = [errorXByIterr,     calPos(1)-userECEFPosition(1)];
    errorYByIterr      = [errorYByIterr,     calPos(2)-userECEFPosition(2)];
    errorZByIterr      = [errorZByIterr,     calPos(3)-userECEFPosition(3)];
    B0ByIterr          = [B0ByIterr,         B00List(end)];
    PDOPByIterr        = [PDOPByIterr,       PDOP];
    HDOPByIterr        = [HDOPByIterr,       HDOP];
    VDOPByIterr        = [VDOPByIterr,       VDOP];
    % x0List, y0List, z0List
    for itera = 1:length(x0List)
        errorXIteraMatrix(iterr,itera) = x0List(itera)-userECEFPosition(1);
        errorYIteraMatrix(iterr,itera) = y0List(itera)-userECEFPosition(2);
        errorZIteraMatrix(iterr,itera) = z0List(itera)-userECEFPosition(3);
    end
    if(rem(iterr,50) == 0)
        disp(['第', num2str(iterr), '次仿真结束']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
errorByIterrSort      = sort(errorByIterr);
errorHoriByIterrSort  = sort(errorHoriByIterr);
errorVetiByIterrJSort = sort(errorVetiByIterrJ);
errorEastByIterrJSort = sort(errorEastByIterrJ);
errorNorthByIterrJSort = sort(errorNorthByIterrJ);

%%%%%%%% 95%误差+DOP %%%%%%%%
disp(['error(95%)       = ', num2str(errorByIterrSort(paraAll.iterrT*0.95))]);
disp(['errorHori(95%)   = ', num2str(errorHoriByIterrSort(paraAll.iterrT*0.95))]);    
disp(['errorVetiJ(95%)  = ', num2str(errorVetiByIterrJSort(paraAll.iterrT*0.95))]); 
disp(['errorEastJ(95%)  = ', num2str(errorEastByIterrJSort(paraAll.iterrT*0.95))]); 
disp(['errorNorthJ(95%) = ', num2str(errorNorthByIterrJSort(paraAll.iterrT*0.95))]); 
disp(['error(65%)       = ', num2str(errorByIterrSort(paraAll.iterrT*0.65))]);
disp(['errorHori(65%)   = ', num2str(errorHoriByIterrSort(paraAll.iterrT*0.65))]);    
disp(['errorVetiJ(65%)  = ', num2str(errorVetiByIterrJSort(paraAll.iterrT*0.65))]); 
disp(['errorEastJ(65%)  = ', num2str(errorEastByIterrJSort(paraAll.iterrT*0.65))]); 
disp(['errorNorthJ(65%) = ', num2str(errorNorthByIterrJSort(paraAll.iterrT*0.65))]); 

disp(['PDOP             = ', num2str(mean(PDOPByIterr))]);
disp(['HDOP             = ', num2str(mean(HDOPByIterr))]);
disp(['VDOP             = ', num2str(mean(VDOPByIterr))]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%    SEP     %%%%%%%%
[Vg_s, DOP1, DOP2, DOP3] = calSEP(...
paraAll.f0,...
userECEFPosition,...
satECEFPositionList, satECEFVelocityList,[],...                                                           
paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
paraAll.withVORFlag, paraAll.withDMEFlag, paraAll.withHeightFlag,...
paraAll.withClcShift, paraAll.withClcDraft);
disp('Vg_s:');
matrixSout(Vg_s)
disp(['DOP1           = ', num2str(DOP1)]);
disp(['DOP2           = ', num2str(DOP2)]);
disp(['DOP3           = ', num2str(DOP3)]);
if(paraAll.withLEOPrFlag == 1 && paraAll.withLEOPdFlag == 0)
    v1 = Vg_s(:,1)*DOP1*paraAll.pseudoRangeError;
    v2 = Vg_s(:,2)*DOP2*paraAll.pseudoRangeError;
    v3 = Vg_s(:,3)*DOP3*paraAll.pseudoRangeError;
elseif(paraAll.withLEOPrFlag == 0 && paraAll.withLEOPdFlag == 1)
    v1 = Vg_s(:,1)*DOP1*paraAll.dopplerError;
    v2 = Vg_s(:,2)*DOP2*paraAll.dopplerError;
    v3 = Vg_s(:,3)*DOP3*paraAll.dopplerError;
elseif(paraAll.withLEOPrFlag == 1 && paraAll.withLEOPdFlag == 1)
    v1 = Vg_s(:,1)*DOP1;
    v2 = Vg_s(:,2)*DOP2;
    v3 = Vg_s(:,3)*DOP3;
else
end
lon0rad = userLLHPosition(1)*pi/180;
lat0rad = userLLHPosition(2)*pi/180;
R       = [-sin(lon0rad)         ,  cos(lon0rad)             , 0;
       -sin(lat0rad)*cos(lon0rad), -sin(lat0rad)*sin(lon0rad), cos(lat0rad);
        cos(lat0rad)*cos(lon0rad),  cos(lat0rad)*sin(lon0rad), sin(lat0rad)];
v1_ENU  = R * v1;
v2_ENU  = R * v2;
v3_ENU  = R * v3;
disp('三轴转到东北天方向上的结果')
R_Vg_s = R*Vg_s;
matrixSout(R_Vg_s);
disp(['sigmaEast  = ', num2str(sqrt((R_Vg_s(1,1)*DOP1)^2 + (R_Vg_s(1,2)*DOP2)^2 +(R_Vg_s(1,3)*DOP3)^2))])
disp(['sigmaNorth = ', num2str(sqrt((R_Vg_s(2,1)*DOP1)^2 + (R_Vg_s(2,2)*DOP2)^2 +(R_Vg_s(2,3)*DOP3)^2))])
disp(['sigmaVeti  = ', num2str(sqrt((R_Vg_s(3,1)*DOP1)^2 + (R_Vg_s(3,2)*DOP2)^2 +(R_Vg_s(3,3)*DOP3)^2))])
disp("验证：")
vec     = v1_ENU + v2_ENU + v3_ENU;
disp(['|vec| = ', num2str(norm(vec))]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%    CRLB     %%%%%%%%
CRLBMatrix = calCRLB(...
paraAll.f0,...
userECEFPosition,...
satECEFPositionList, satECEFVelocityList,[],...                                                                     
paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
paraAll.withVORFlag, paraAll.withDMEFlag, paraAll.withHeightFlag,...
paraAll.withClcShift, paraAll.withClcDraft);

posCRLB = sqrt(CRLBMatrix(1,1)+CRLBMatrix(2,2)+CRLBMatrix(3,3));
posCRLB_2 = 2*sqrt(CRLBMatrix(1,1)+CRLBMatrix(2,2)+CRLBMatrix(3,3));
disp(['posCRLB_2 =      ', num2str(posCRLB_2)]);

errorX_var = var(errorXByIterr);
errorY_var = var(errorYByIterr);
errorZ_var = var(errorZByIterr);
disp(['error_std =      ', num2str(sqrt(errorX_var+errorY_var+errorZ_var))]);
disp(['posCRLB =        ', num2str(posCRLB)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% 误差分布（ENU + EN）%%%%%%%%
scale       = 1;
lineWidth   = 3;
maxHeadSize = 1;
fontSize    = 20;
pointSize   = 6;

range1      = 50;
range2      = 50;
figure(2)
    subplot(1,3,1)
    scatter3(errorEastByIterr, errorNorthByIterr, errorVetiByIterr, pointSize, 'filled', 'r');
%     hold on;
%     quiver3(0,0,0,scale*v1_ENU(1,1),scale*v1_ENU(2,1),scale*v1_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
%     hold on;
%     quiver3(0,0,0,scale*v2_ENU(1,1),scale*v2_ENU(2,1),scale*v2_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
%     hold on; 
%     quiver3(0,0,0,scale*v3_ENU(1,1),scale*v3_ENU(2,1),scale*v3_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
%     hold on;
%     quiver3(0,0,0,-scale*v1_ENU(1,1),-scale*v1_ENU(2,1),-scale*v1_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
%     hold on;
%     quiver3(0,0,0,-scale*v2_ENU(1,1),-scale*v2_ENU(2,1),-scale*v2_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
%     hold on; 
%     quiver3(0,0,0,-scale*v3_ENU(1,1),-scale*v3_ENU(2,1),-scale*v3_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
    xlim([-range1 range1])
    ylim([-range1 range1])
    zlim([-range1 range1])

    xlabel('东向误差/(m)', FontSize=fontSize);
    ylabel('北向误差/(m)', FontSize=fontSize);
    zlabel('天向误差/(m)', FontSize=fontSize);
    grid on;
    title('定位误差', FontSize=fontSize);
    
	subplot(1,3,2)  
    scatter(errorEastByIterr, errorNorthByIterr, pointSize, 'filled', 'r');
    xlabel('东向误差/(m)', FontSize=fontSize);
    ylabel('北向误差/(m)', FontSize=fontSize);
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-range1 range1]);
    ylim([-range1 range1]);
    grid on;
    title('水平定位误差', FontSize=fontSize);

    subplot(1,3,3)
    scatter(errorNorthByIterr, errorVetiByIterr, pointSize, 'filled', 'r');
    xlabel('北向误差/(m)', FontSize=fontSize);
    ylabel('天向误差/(m)', FontSize=fontSize);
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-range2 range2]);
    ylim([-range2 range2]);
    grid on;
    title('北向天向误差', FontSize=fontSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%  误差回归过程  %%%%%%%%
% errorXIteraMatrix errorYIteraMatrix errorZIteraMatrix
% r Red g Green b Blue c Cyan m Magenta y Yellow k Black w White
% figure(3)  
%     subplot(1,2,1);
%     scatter3(errorXByIterr,errorYByIterr,errorZByIterr,'r');
%     hold on;
%     scatter3(errorXIteraMatrix(:,1),errorYIteraMatrix(:,1),errorZIteraMatrix(:,1),'g');
%     hold on;
%     scatter3(errorXIteraMatrix(:,2),errorYIteraMatrix(:,2),errorZIteraMatrix(:,2),'b');
%     hold on;
%     scatter3(errorXIteraMatrix(:,3),errorYIteraMatrix(:,3),errorZIteraMatrix(:,3),'c');
%     hold on;
%     scatter3(errorXIteraMatrix(:,4),errorYIteraMatrix(:,4),errorZIteraMatrix(:,4),'m');
%     hold on;
%     scatter3(errorXIteraMatrix(:,5),errorYIteraMatrix(:,5),errorZIteraMatrix(:,5),'y');
%     xlabel('errorX/(m)');
%     ylabel('errorY/(m)');
%     zlabel('errorZ/(m)');
%     ax1 = gca;
%     ax1.XAxisLocation = 'origin';
%     ax1.YAxisLocation = 'origin';
%     
%     subplot(1,2,2);
%     scatter(errorXByIterr,errorYByIterr,'r');
%     hold on;
%     scatter(errorXIteraMatrix(:,1),errorYIteraMatrix(:,1),'g');
%     hold on;
%     scatter(errorXIteraMatrix(:,2),errorYIteraMatrix(:,2),'b');
%     hold on;
%     scatter(errorXIteraMatrix(:,3),errorYIteraMatrix(:,3),'c');
%     hold on;
%     scatter(errorXIteraMatrix(:,4),errorYIteraMatrix(:,4),'m');
%     hold on;
%     scatter(errorXIteraMatrix(:,5),errorYIteraMatrix(:,5),'y');
%     xlabel('errorX/(m)');
%     ylabel('errorY/(m)');
%     ax1 = gca;
%     ax1.XAxisLocation = 'origin';
%     ax1.YAxisLocation = 'origin';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
         
toc