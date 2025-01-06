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
[ttList_ori, satECEFPositionList, satECEFVelocityList,...
      rangeList_ori, dopplerList_ori] = getLEOObsData(paraAll.obsDataFilePath);
Ns = length(rangeList_ori);
disp(['多历元卫星观测量的个数为：',num2str(Ns)]);

[azimuthList, elevationList] = calAziAndEle(userLLHPosition, satECEFPositionList);
figure
    starsky(azimuthList, elevationList, '.');
%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%得到初始值%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(paraAll.positionStartWay == 0)
    userECEFPositionStart = getInitialPos(rangeList_ori, satECEFPositionList);
end
if(paraAll.positionStartWay == 1)
    userECEFPositionStart = [userECEFPosition(1) + paraAll.errorxyz, userECEFPosition(2) + paraAll.errorxyz, userECEFPosition(3) + paraAll.errorxyz];
end
if(paraAll.positionStartWay == 2)
    userECEFPositionStart = [0, 0, 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errorByClcDraft   = [];
for clcDraftIndex = 1:length(paraAll.clcDraftList)
    % ttList_ori, rangeList_ori, dopplerList_ori
    rangeList   = [];
    dopplerList = [];
    for obsIndex = 1:length(rangeList_ori)
        clcError    = paraAll.beta0+paraAll.clcDraftList(clcDraftIndex)*(ttList_ori(obsIndex)-paraAll.timeStart)+paraAll.beta2*(ttList_ori(obsIndex)-paraAll.timeStart)^2;
        rangeList   = [rangeList, rangeList_ori(obsIndex)+clcError];
        fError      = paraAll.clcDraftList(clcDraftIndex) + 2*paraAll.beta2*(ttList_ori(obsIndex)-paraAll.timeStart);
        dopplerList = [dopplerList, dopplerList_ori(obsIndex)+fError];
    end
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
        for ns = 1:Ns
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
            prSatList, pdSatList, satECEFPositionList, satECEFVelocityList,... % 卫星观测量
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
    sum = 0;
    for i = 1:paraAll.iterrT
        sum = sum + (errorByIterr(i))^2;
    end
    error_std       = sqrt(sum/paraAll.iterrT);
    errorByClcDraft = [errorByClcDraft, error_std];

%     figure(2)
%     scatter(errorEastByIterr, errorNorthByIterr);
%     xlabel('东向误差/(m)');
%     ylabel('北向误差/(m)');
%     title('水平定位误差');
%     ax1 = gca;
%     ax1.XAxisLocation = 'origin';
%     ax1.YAxisLocation = 'origin';
% 
%     figure(3)
%     scatter3(errorXByIterr, errorYByIterr, errorZByIterr);
%     xlabel('X误差/(m)');
%     ylabel('Y误差/(m)');
%     zlabel('Z误差/(m)');
%     title('定位误差');
%     ax2 = gca;
%     ax2.XAxisLocation = 'origin';
%     ax2.YAxisLocation = 'origin';
    disp(['error_std = ', num2str(error_std)]);
    disp(['第', num2str(clcDraftIndex), '组钟漂的仿真结束']);
    disp("--------------------")
end
CRLBMatrix = calCRLB(...
    paraAll.f0,...
    userECEFPosition,...
    satECEFPositionList, satECEFVelocityList,... 
    [],...                                                                     
    paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
    paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
    paraAll.withVORFlag, paraAll.withDMEFlag, paraAll.withHeightFlag);
errorCRLB = sqrt(CRLBMatrix(1,1)+CRLBMatrix(2,2)+CRLBMatrix(3,3));
disp(['errorCRLB = ', num2str(errorCRLB)]);
errorCRLBList = errorCRLB*ones(1, length(paraAll.clcDraftList));

figure
    plot(paraAll.clcDraftList, errorByClcDraft, 'r');
    hold on;
    plot(paraAll.clcDraftList, errorCRLBList, 'g');
    xlabel('clcDrift/(s)');
    ylabel('error/(m)');
    legend(['error', 'CRLB']);
    title('误差随钟漂变化');

toc