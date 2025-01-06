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
starsky(azimuthList, elevationList, '.');
%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%得到初始值%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(paraAll.positionStartWay == 0)
    userECEFPositionStart = getInitialPos(rangeList, satECEFPositionList);
end
if(paraAll.positionStartWay == 1)
    userECEFPositionStart = [userECEFPosition(1) + paraAll.errorxyz, userECEFPosition(2) + paraAll.errorxyz, userECEFPosition(3) + paraAll.errorxyz];
end
if(paraAll.positionStartWay == 2)
    userECEFPositionStart = [0, 0, 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errorMatrix        = [];   % obsTime x iterrT
errorXMatrix       = [];
errorYMatrix       = [];
errorZMatrix       = [];

errorStdByObsTime  = [];
errorMeanByObsTime = [];
error65ByObsTime   = [];  % 65.3 
error95ByObsTime   = [];  % 95.3 
errorCRLBByObsTime = [];
PDOPByObsTime      = [];
HDOPByObsTime      = [];
VDOPByObsTime      = [];
for obsTime = paraAll.obsTimeList
%     ttList, satECEFPositionList, satECEFVelocityList,...
%       rangeList, dopplerList
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
    for iterr = 1:paraAll.iterrT
        %%%%%%%%%%%%%%%% 加观测误差%%%%%%%%%%%%%%%%
        % input:rangeList dopplerList
        prSatList = [];
        pdSatList = [];
        for ns = 1:obsTime
            errorPr   = paraAll.pseudoRangeError*randn;
            prSat     = rangeList(ns) + errorPr + paraAll.satError*randn + paraAll.ionoError*randn;
            prSatList = [prSatList, prSat];
            errorPd   = paraAll.dopplerError*randn;
            pdSat     = dopplerList(ns) + errorPd;
            pdSatList = [pdSatList, pdSat];
        end
        % output:prSatList pdSatList

        % input:userLLHPosition(3)
        height = userLLHPosition(3) + (paraAll.hError + paraAll.hErrorSys)*randn;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%LSM%%%%%%%%%%%%%%%%%%%%%
        [x0List, y0List, z0List, B00List, B10List,...
        PDOP, HDOP, VDOP, iterTime] = ...
        LSM60(paraAll.f0,...
            prSatList, pdSatList, satECEFPositionList(1:obsTime,:), satECEFVelocityList(1:obsTime,:),... % 卫星观测量
            [], [], [],...                 % VOR/DME观测量
            height,...                                                         % 高度观测量
            paraAll.pseudoRangeError, paraAll.dopplerError, 0, 0, paraAll.hError,...% 观测误差
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
        PDOPByIterr       = [PDOPByIterr,       PDOP];
        HDOPByIterr       = [HDOPByIterr,       HDOP];
        VDOPByIterr       = [VDOPByIterr,       VDOP];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    errorMatrix  = [errorMatrix;  errorByIterr];
    errorXMatrix = [errorXMatrix; errorXByIterr];
    errorYMatrix = [errorYMatrix; errorYByIterr];
    errorZMatrix = [errorZMatrix; errorZByIterr];
    sum = 0;
    for i = 1:paraAll.iterrT
        sum = sum + (errorByIterr(i))^2;
    end
    errorByIterrSort = sort(errorByIterr);
    errorStd         = sqrt(sum/paraAll.iterrT);
    
    CRLBMatrix = calCRLB(...
        paraAll.f0,...
        userECEFPosition,...
        satECEFPositionList(1:obsTime,:), satECEFVelocityList(1:obsTime,:),[],...                                                                      
        paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
        paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
        paraAll.withVORFlag, paraAll.withDMEFlag, paraAll.withHeightFlag);
    errorCRLB = sqrt(CRLBMatrix(1,1)+CRLBMatrix(2,2)+CRLBMatrix(3,3));
    
    errorStdByObsTime  = [errorStdByObsTime, errorStd];
    errorMeanByObsTime = [errorMeanByObsTime, mean(errorByIterr)];
    error65ByObsTime   = [error65ByObsTime, errorByIterrSort(paraAll.iterrT*0.653)];  % 65.3 
    error95ByObsTime   = [error95ByObsTime, errorByIterrSort(paraAll.iterrT*0.95)];  % 95.3 
    errorCRLBByObsTime = [errorCRLBByObsTime, errorCRLB];
   
    disp(['errorStd  = ', num2str(errorStd)]);
    disp(['errorMean = ', num2str(mean(errorByIterr))]);
    disp(['error65   = ', num2str(errorByIterrSort(paraAll.iterrT*0.653))]);
    disp(['error95   = ', num2str(errorByIterrSort(paraAll.iterrT*0.95))]);
    disp(['errorCRLB = ', num2str(errorCRLB)]);
    
    disp(['观测时长为', num2str(obsTime), '的仿真结束']);
    disp("--------------------");
end
% errorMatrix errorXMatrix errorYMatrix errorZMatrix
errorVarList  = [];
errorXVarList = [];
errorYVarList = [];
errorZVarList = [];
for i = 1:length(paraAll.obsTimeList)
    varX2 = 0;
    varY2 = 0;
    varZ2 = 0;
    for j = 1:paraAll.iterrT
        varX2 = varX2 + errorXMatrix(i,j)^2/paraAll.iterrT;
        varY2 = varY2 + errorYMatrix(i,j)^2/paraAll.iterrT;
        varZ2 = varZ2 + errorZMatrix(i,j)^2/paraAll.iterrT;
    end
    errorVarList  = [errorVarList,  sqrt(varX2+varY2+varZ2)];
    errorXVarList = [errorXVarList, var(errorXMatrix(i,:))];
    errorYVarList = [errorYVarList, var(errorYMatrix(i,:))];
    errorZVarList = [errorZVarList, var(errorZMatrix(i,:))];
end
figure
    plot(paraAll.obsTimeList, errorVarList, 'r');
    hold on;
    plot(paraAll.obsTimeList, errorCRLBByObsTime, 'g');
    xlabel('观测时长/(s)');
    ylabel('定位误差/(m)');
    legend('RMS', 'CRLB');
    title('误差随观测时长的变化');

figure(1)
    plot(paraAll.obsTimeList, error65ByObsTime, 'r');
    hold on;
    plot(paraAll.obsTimeList, errorCRLBByObsTime, 'g');
    xlabel('观测时长/(s)');
    ylabel('定位误差/(m)');
    legend('RMS', 'CRLB');
    title('误差随观测时长的变化');
figure(2)
    plot(paraAll.obsTimeList, error95ByObsTime/2, 'r');
    hold on;
    plot(paraAll.obsTimeList, errorCRLBByObsTime, 'g');
    xlabel('观测时长/(s)');
    ylabel('定位误差/(m)');
    legend('RMS', 'CRLB');
    title('误差随观测时长的变化');
figure(3)
    plot(paraAll.obsTimeList, errorStdByObsTime, 'r');
    hold on;
    plot(paraAll.obsTimeList, errorMeanByObsTime, 'g');
    hold on;
    plot(paraAll.obsTimeList, error65ByObsTime, 'c');
    hold on;
    plot(paraAll.obsTimeList, error95ByObsTime, 'm');
    hold on;
    plot(paraAll.obsTimeList, errorCRLBByObsTime, 'k');
    xlabel('观测时长/(s)');
    ylabel('定位误差/(m)');
    legend('RMS','Mean','65%','95%', 'CRLB');
    title('误差随观测时长的变化');
figure(4)
    plot(paraAll.obsTimeList, errorStdByObsTime, 'r', 'Linewidth', 2);
    hold on;
    plot(paraAll.obsTimeList, errorCRLBByObsTime, 'g', 'Linewidth', 2);
    xlabel('观测时长/(s)', FontSize=15);
    ylabel('定位误差/(m)', FontSize=15);
    legend('RMS', 'CRLB', FontSize=15);
    title('误差随观测时长的变化', FontSize=15);

% dlmwrite('./result/obsTime/errorMatrix.txt',errorMatrix);
% dlmwrite('./result/obsTime/errorXMatrix.txt',  errorXMatrix);
% dlmwrite('./result/obsTime/errorYMatrix.txt' ,  errorYMatrix);
% dlmwrite('./result/obsTime/errorZMatrix.txt' , errorZMatrix );
% 
% dlmwrite('./result/obsTime/errorStdByObsTime.txt' ,errorStdByObsTime);
% dlmwrite('./result/obsTime/errorMeanByObsTime.txt', errorMeanByObsTime); 
% dlmwrite('./result/obsTime/error65ByObsTime.txt' ,   error65ByObsTime);
% dlmwrite('./result/obsTime/error95ByObsTime.txt', error95ByObsTime  );
% dlmwrite('./result/obsTime/errorCRLBByObsTime.txt', errorCRLBByObsTime);
toc