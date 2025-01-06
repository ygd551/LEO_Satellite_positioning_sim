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

[azimuthList, elevationList] = calAziAndEle(userLLHPosition, satECEFPositionList);
figure(1)
    starsky(azimuthList, elevationList, '.');
%% 取出观测量%%%%%%%%%%%%%%%%%%%%
%%%%% 所有理论上可见的VOR/DME得到的观测量%%%%%%
VORDMEDataAll               = getVORDMEData(userECEFPosition);
obsVORDMEECEFPostionListAll = VORDMEDataAll.obsVORDMEECEFPostionList;
realDisListAll              = VORDMEDataAll.realDisList;
realAngleListAll            = VORDMEDataAll.realAngleList;
Nvd                         = length(realDisListAll);
disp(['用户可见陆基台站个数：',num2str(Nvd)]);
% 给出所有台站相对于用户的位置图
enuList = [];
for i = 1 : length(obsVORDMEECEFPostionListAll)
    enu = ecef2enu(obsVORDMEECEFPostionListAll(i,:), userECEFPosition, userLLHPosition);
    enuList = [enuList; enu];
end
index1 = 2;
index2 = 5;
obsVORDMEECEFPostionList = [obsVORDMEECEFPostionListAll(index1, :);
    obsVORDMEECEFPostionListAll(index2, :)];
realDisList = [realDisListAll(index1), realDisListAll(index2)];
realAngleList = [realAngleListAll(index1), realAngleListAll(index2)];
disp(['选取DME个数', num2str(length(realDisList))]);
disp("==============================================")
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
errorXIteraMatrix = zeros(paraAll.iterrT, 10);
errorYIteraMatrix = zeros(paraAll.iterrT, 10);
errorZIteraMatrix = zeros(paraAll.iterrT, 10);
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
    prDMEList = [];
    paVORList = [];
    for nvd = 1:2
        prDME     = realDisList(nvd) + paraAll.DMEError*randn;
        prDMEList = [prDMEList, prDME];
        paVOR     = realAngleList(nvd) * (1 + paraAll.VORError*randn) * pi/180;
        paVORList = [paVORList, paVOR];
    end
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
    prDMEList, paVORList, obsVORDMEECEFPostionList,...                 % VOR/DME观测量
    height,...                                                         % 高度观测量
    paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
    paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
    paraAll.withVORFlag, paraAll.withDMEFlag, paraAll.withHeightFlag, paraAll.weightFlag,...
    userECEFPositionStart, 0, 0, paraAll.clcOrder); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    PDOPByIterr       = [PDOPByIterr,        PDOP];
    HDOPByIterr       = [HDOPByIterr,       HDOP];
    VDOPByIterr       = [VDOPByIterr,       VDOP];
    % x0List, y0List, z0List
    for itera = 1:length(x0List)
        errorXIteraMatrix(iterr,itera) = x0List(itera)-userECEFPosition(1);
        errorYIteraMatrix(iterr,itera) = y0List(itera)-userECEFPosition(2);
        errorZIteraMatrix(iterr,itera) = z0List(itera)-userECEFPosition(3);
    end
    disp(['第', num2str(iterr), '次仿真结束']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
errorByIterrSort     = sort(errorByIterr);
errorHoriByIterrSort = sort(errorHoriByIterr);
errorVetiByIterrSort = sort(errorVetiByIterr);

disp(['error(95%) = ',     num2str(errorByIterrSort(paraAll.iterrT*0.95))]);
disp(['errorHori(95%) = ', num2str(errorHoriByIterrSort(paraAll.iterrT*0.95))]);    
disp(['errorVeti(95%) = ', num2str(errorVetiByIterrSort(paraAll.iterrT*0.95))]); 
disp(['PDOP = ',           num2str(mean(PDOPByIterr))]);
disp(['HDOP = ',           num2str(mean(HDOPByIterr))]);
disp(['VDOP = ',           num2str(mean(VDOPByIterr))]);

CRLBMatrix = calCRLB(...
paraAll.f0,...
userECEFPosition,...
satECEFPositionList, satECEFVelocityList,... 
obsVORDMEECEFPostionList,...                                                                     
paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
paraAll.withVORFlag, paraAll.withDMEFlag, paraAll.withHeightFlag);

posCRLB_2 = 2*sqrt(CRLBMatrix(1,1)+CRLBMatrix(2,2)+CRLBMatrix(3,3));
disp(['posCRLB_2 = ', num2str(posCRLB_2)]);

errorX_var = var(errorXByIterr);
errorY_var = var(errorYByIterr);
errorZ_var = var(errorZByIterr);
posCRLB = sqrt(CRLBMatrix(1,1)+CRLBMatrix(2,2)+CRLBMatrix(3,3));
disp(['error_std = ', num2str(sqrt(errorX_var+errorY_var+errorZ_var))]);
disp(['posCRLB = ', num2str(posCRLB)]);


figure(2)
    subplot(1,2,1);
    scatter(errorEastByIterr, errorNorthByIterr);
    xlabel('东向误差/(m)');
    ylabel('北向误差/(m)');
    title('水平定位误差');
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    subplot(1,2,2);
    scatter3(errorXByIterr, errorYByIterr, errorZByIterr);
    xlabel('X误差/(m)');
    ylabel('Y误差/(m)');
    zlabel('Z误差/(m)');
    title('定位误差');
    ax2 = gca;
    ax2.XAxisLocation = 'origin';
    ax2.YAxisLocation = 'origin';

% errorXIteraMatrix errorYIteraMatrix errorZIteraMatrix
figure(4)
    scatter3(errorXIteraMatrix(1,:),errorYIteraMatrix(1,:),errorZIteraMatrix(1,:),'r');
    hold on;
    scatter3(errorXIteraMatrix(2,:),errorYIteraMatrix(2,:),errorZIteraMatrix(2,:),'g');
    hold on;
    scatter3(errorXIteraMatrix(3,:),errorYIteraMatrix(3,:),errorZIteraMatrix(3,:),'b');
    hold on;
    scatter3(errorXIteraMatrix(4,:),errorYIteraMatrix(4,:),errorZIteraMatrix(4,:),'c');
    hold on;
    scatter3(errorXIteraMatrix(5,:),errorYIteraMatrix(5,:),errorZIteraMatrix(5,:),'m');

% r Red g Green b Blue c Cyan m Magenta y Yellow k Black w White
figure(5)  
    subplot(1,2,1);
    scatter3(errorXByIterr,errorYByIterr,errorZByIterr,'r');
    hold on;
    scatter3(errorXIteraMatrix(:,1),errorYIteraMatrix(:,1),errorZIteraMatrix(:,1),'g');
    hold on;
    scatter3(errorXIteraMatrix(:,2),errorYIteraMatrix(:,2),errorZIteraMatrix(:,2),'b');
    hold on;
    scatter3(errorXIteraMatrix(:,3),errorYIteraMatrix(:,3),errorZIteraMatrix(:,3),'c');
    hold on;
    scatter3(errorXIteraMatrix(:,4),errorYIteraMatrix(:,4),errorZIteraMatrix(:,4),'m');
    hold on;
    scatter3(errorXIteraMatrix(:,5),errorYIteraMatrix(:,5),errorZIteraMatrix(:,5),'y');
    xlabel('errorX/(m)');
    ylabel('errorY/(m)');
    zlabel('errorZ/(m)');
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    
    subplot(1,2,2);
    scatter(errorXByIterr,errorYByIterr,'r');
    hold on;
    scatter(errorXIteraMatrix(:,1),errorYIteraMatrix(:,1),'g');
    hold on;
    scatter(errorXIteraMatrix(:,2),errorYIteraMatrix(:,2),'b');
    hold on;
    scatter(errorXIteraMatrix(:,3),errorYIteraMatrix(:,3),'c');
    hold on;
    scatter(errorXIteraMatrix(:,4),errorYIteraMatrix(:,4),'m');
    hold on;
    scatter(errorXIteraMatrix(:,5),errorYIteraMatrix(:,5),'y');
    xlabel('errorX/(m)');
    ylabel('errorY/(m)');
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    
errorXIteraMatrix(1:5,:)    
    
       
toc