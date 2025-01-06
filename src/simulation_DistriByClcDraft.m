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
%%%%%%%%%%%%%%%%%%得到初始值%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(paraAll.positionStartWay == 0)
    userECEFPositionStart = getInitialPos(rangeList_ori, satECEFPositionList);
end
if(paraAll.positionStartWay == 1)
    userECEFPositionStart = [userECEFPosition(1) + paraAll.errorxyz(1), userECEFPosition(2) + paraAll.errorxyz(2), userECEFPosition(3) + paraAll.errorxyz(3)];
end
if(paraAll.positionStartWay == 2)
    userECEFPositionStart = [0, 0, 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errorMatrix      = [];
errorHoriMatrix  = [];
errorEastMatrix  = [];
errorNorthMatrix = [];
errorVetiMatrix  = [];
errorXMatrix     = [];
errorYMatrix     = [];
errorZMatrix     = [];
PDOPMatrix       = [];
HDOPMatrix       = [];
VDOPMatrix       = [];
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
    disp('==============================')
    disp(['钟漂为', num2str(paraAll.clcDraftList(clcDraftIndex)), '仿真开始']);
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
        PDOPByIterr       = [PDOPByIterr,       PDOP];
        HDOPByIterr       = [HDOPByIterr,       HDOP];
        VDOPByIterr       = [VDOPByIterr,       VDOP];
        if(rem(iterr,50) == 0)
            disp(['第', num2str(iterr), '次仿真结束']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    errorByIterrSort     = sort(errorByIterr);
    errorHoriByIterrSort = sort(errorHoriByIterr);
    errorVetiByIterrSort = sort(errorVetiByIterr);

    disp(['error(95%) =     ', num2str(errorByIterrSort(paraAll.iterrT*0.95))]);
    disp(['errorHori(95%) = ', num2str(errorHoriByIterrSort(paraAll.iterrT*0.95))]);    
    disp(['errorVeti(95%) = ', num2str(errorVetiByIterrSort(paraAll.iterrT*0.95))]); 
    disp(['PDOP =           ', num2str(mean(PDOPByIterr))]);
    disp(['HDOP =           ', num2str(mean(HDOPByIterr))]);
    disp(['VDOP =           ', num2str(mean(VDOPByIterr))]);
    %%% 结果统计
    errorMatrix      = [errorMatrix;      errorByIterr     ];
    errorHoriMatrix  = [errorHoriMatrix;  errorHoriByIterr ];
    errorEastMatrix  = [errorEastMatrix;  errorEastByIterr ];
    errorNorthMatrix = [errorNorthMatrix; errorNorthByIterr];
    errorVetiMatrix  = [errorVetiMatrix;  errorVetiByIterr ];
    errorXMatrix     = [errorXMatrix;     errorXByIterr    ];
    errorYMatrix     = [errorYMatrix;     errorYByIterr    ];
    errorZMatrix     = [errorZMatrix;     errorZByIterr    ];
    PDOPMatrix       = [PDOPMatrix;       PDOPByIterr      ];
    HDOPMatrix       = [HDOPMatrix;       HDOPByIterr      ];
    VDOPMatrix       = [VDOPMatrix;       VDOPByIterr      ];
end

figure(1)
    scatter3(errorXMatrix(1,:),errorYMatrix(1,:),errorZMatrix(1,:),'r');
    hold on;
    scatter3(errorXMatrix(2,:),errorYMatrix(2,:),errorZMatrix(2,:),'g');
    hold on;
    scatter3(errorXMatrix(3,:),errorYMatrix(3,:),errorZMatrix(3,:),'y');
    legend('B1 = 0.03(m/s)', 'B1 = 0.13(m/s)', 'B1 = 0.23(m/s)');
    xlabel('errorX/(m)');
    ylabel('errorY/(m)');
    zlabel('errorZ/(m)');
    title('钟漂大小对伪距+多普勒定位误差分布的影响');
figure(2)
    scatter(errorEastMatrix(1,:),errorEastMatrix(1,:),'r');
    hold on;
    scatter(errorEastMatrix(2,:),errorEastMatrix(2,:),'g');
    hold on;
    scatter(errorEastMatrix(3,:),errorEastMatrix(3,:),'y');
    xlabel('东向误差/(m)');
    ylabel('北向误差/(m)');
    legend('B1 = 0.03(m/s)', 'B1 = 0.13(m/s)', 'B1 = 0.23(m/s)');
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    title('钟漂大小对伪距+多普勒定位水平误差分布的影响');
figure(3)
    scatter(errorXMatrix(1,:),errorYMatrix(1,:),'r');
    hold on;
    scatter(errorXMatrix(2,:),errorYMatrix(2,:),'g');
    hold on;
    scatter(errorXMatrix(3,:),errorYMatrix(3,:),'y');
    legend('B1 = 0.03(m/s)', 'B1 = 0.13(m/s)', 'B1 = 0.23(m/s)');
    xlabel('errorX/(m)');
    ylabel('errorY/(m)');
    title('钟漂大小对伪距+多普勒对定位误差分布的影响(X/Y)');

toc