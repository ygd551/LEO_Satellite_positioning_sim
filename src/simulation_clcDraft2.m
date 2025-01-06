tic
clc;clear;
format long g

addpath(genpath("../../coordinateTransformation"));
paraAll = simulation_clcDraft2_para();
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


clcDraftList        = [0.03, 0.09, 0.15, 0.21,...
                       0.27, 0.33, 0.39, 0.45,...
                       0.51, 0.57, 0.63];

clcDraftLength      = length(clcDraftList);


errorByIterrM       = zeros(clcDraftLength, paraAll.iterrT);
errorHoriByIterrM   = zeros(clcDraftLength, paraAll.iterrT);
errorEastByIterrM   = zeros(clcDraftLength, paraAll.iterrT);
errorEastByIterrJM  = zeros(clcDraftLength, paraAll.iterrT);
errorNorthByIterrM  = zeros(clcDraftLength, paraAll.iterrT);
errorNorthByIterrJM = zeros(clcDraftLength, paraAll.iterrT);
errorVetiByIterrM   = zeros(clcDraftLength, paraAll.iterrT);
errorVetiByIterrJM  = zeros(clcDraftLength, paraAll.iterrT);
errorXByIterrM      = zeros(clcDraftLength, paraAll.iterrT);
errorYByIterrM      = zeros(clcDraftLength, paraAll.iterrT);
errorZByIterrM      = zeros(clcDraftLength, paraAll.iterrT);
PDOPByIterrM        = zeros(clcDraftLength, paraAll.iterrT);
HDOPByIterrM        = zeros(clcDraftLength, paraAll.iterrT);
VDOPByIterrM        = zeros(clcDraftLength, paraAll.iterrT);
B0ByIterrM          = zeros(clcDraftLength, paraAll.iterrT);

error95List         = zeros(1, clcDraftLength);
error65List         = zeros(1, clcDraftLength);
errorHori95List     = zeros(1, clcDraftLength);
errorHori65List     = zeros(1, clcDraftLength);
errorVeti95JList    = zeros(1, clcDraftLength);
errorVeti65JList    = zeros(1, clcDraftLength);
errorEastJ95List    = zeros(1, clcDraftLength);
errorEastJ65List    = zeros(1, clcDraftLength);
errorNorthJ95List   = zeros(1, clcDraftLength);
errorNorthJ65List   = zeros(1, clcDraftLength);
  
posCRLBList         = zeros(1, clcDraftLength);
errorStdList        = zeros(1, clcDraftLength);


figure;
for clcDraftIndex = 1:clcDraftLength
    clcDraft            = clcDraftList(clcDraftIndex);
    disp('===============================================')
    disp(['clcDraft为', num2str(clcDraft), '仿真开始']);
    %%%%%%%%%%%%%%%%%%得到初始值%%%%%%%%%%%%%%%%%%%%%
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
    for iterr = 1:paraAll.iterrT
        %%%%%%%%%%%%%%%% 加观测误差%%%%%%%%%%%%%%%%
        % input:rangeListAll dopplerListAll ttListAll
        prSatList = [];
        pdSatList = [];

        errorPrList = [];
        for ns = 1:length(rangeList)
            errorPr   = paraAll.pseudoRangeError*randn + paraAll.satError*randn + paraAll.ionoError*randn;
            prSat     = rangeList(ns) + ttList(ns)*clcDraft + errorPr;
            prSatList = [prSatList, prSat];
            errorPd   = paraAll.dopplerError*randn;
            pdSat     = dopplerList(ns) + clcDraft + errorPd;
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
        errorByIterrM(clcDraftIndex, iterr)       = error;
        errorHoriByIterrM(clcDraftIndex, iterr)   = errorHorizontal;
        errorEastByIterrM(clcDraftIndex, iterr)   = errorEast;
        errorEastByIterrJM(clcDraftIndex, iterr)  = abs(errorEast);
        errorNorthByIterrM(clcDraftIndex, iterr)  = errorNorth;
        errorNorthByIterrJM(clcDraftIndex, iterr) = abs(errorNorth);
        errorVetiByIterrM(clcDraftIndex, iterr)   = errorVertical;
        errorVetiByIterrJM(clcDraftIndex, iterr)  = abs(errorVertical);
        errorXByIterrM(clcDraftIndex, iterr)      = calPos(1)-userECEFPosition(1);
        errorYByIterrM(clcDraftIndex, iterr)      = calPos(2)-userECEFPosition(2);
        errorZByIterrM(clcDraftIndex, iterr)      = calPos(3)-userECEFPosition(3);
        B0ByIterrM(clcDraftIndex, iterr)          = B00List(end);
        PDOPByIterrM(clcDraftIndex, iterr)        = PDOP;
        HDOPByIterrM(clcDraftIndex, iterr)        = HDOP;
        VDOPByIterrM(clcDraftIndex, iterr)        = VDOP;

        if(rem(iterr,200) == 0)
            disp(['第', num2str(iterr), '次仿真结束']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end 

    errorByIterrSort        = sort(errorByIterrM(clcDraftIndex,:));
    errorHoriByIterrSort    = sort(errorHoriByIterrM(clcDraftIndex,:));
    errorVetiByIterrJSort   = sort(errorVetiByIterrJM(clcDraftIndex,:));
    errorEastByIterrJSort   = sort(errorEastByIterrJM(clcDraftIndex,:));
    errorNorthByIterrJSort  = sort(errorNorthByIterrJM(clcDraftIndex,:));

    %%%%%%%% 95%误差+DOP %%%%%%%%
    error95                 =      errorByIterrSort(paraAll.iterrT*0.95);
    disp(['error(95%)       =      ', num2str(error95)]);
    errorHori95             =      errorHoriByIterrSort(paraAll.iterrT*0.95);
    disp(['errorHori(95%)   =      ', num2str(errorHori95)]);    
    errorVetiJ95            =      errorVetiByIterrJSort(paraAll.iterrT*0.95);
    disp(['errorVetiJ(95%)  =      ', num2str(errorVetiJ95)]); 
    errorEastJ95            =      errorEastByIterrJSort(paraAll.iterrT*0.95);
    disp(['errorEastJ(95%)  =      ', num2str(errorEastJ95)]); 
    errorNorthJ95           =      errorNorthByIterrJSort(paraAll.iterrT*0.95);
    disp(['errorNorthJ(95%) =      ', num2str(errorNorthJ95)]); 
    error65                 =      errorByIterrSort(paraAll.iterrT*0.65);
    disp(['error(65%)       =      ', num2str(error65)]);
    errorHori65             =      errorHoriByIterrSort(paraAll.iterrT*0.65);
    disp(['errorHori(65%)   =      ', num2str(errorHori65)]); 
    errorVetiJ65            =      errorVetiByIterrJSort(paraAll.iterrT*0.65);
    disp(['errorVetiJ(65%)  =      ', num2str(errorVetiJ65)]); 
    errorEastJ65            =      errorEastByIterrJSort(paraAll.iterrT*0.65);
    disp(['errorEastJ(65%)  =      ', num2str(errorEastJ65)]); 
    errorNorthJ65           =      errorNorthByIterrJSort(paraAll.iterrT*0.65);
    disp(['errorNorthJ(65%) =      ', num2str(errorNorthJ65)]); 
    
    PDOP                    =      mean(PDOPByIterrM(clcDraftIndex,:));
    disp(['PDOP             =      ', num2str(PDOP)]);
    HDOP                    =      mean(HDOPByIterrM(clcDraftIndex,:));
    disp(['HDOP             =      ', num2str(PDOP)]);
    VDOP                    =      mean(VDOPByIterrM(clcDraftIndex,:));
    disp(['VDOP             =      ', num2str(VDOP)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%    SEP     %%%%%%%%
    [V, L1, L2, L3] = calSEP(...
    paraAll.f0,...
    userECEFPosition,...
    satECEFPositionList, satECEFVelocityList,[],...                                                           
    paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
    paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
    paraAll.withVORFlag, paraAll.withDMEFlag, paraAll.withHeightFlag,...
    paraAll.withClcShift, paraAll.withClcDraft);
    disp('V:');
    matrixSout(V)
    disp(['L1               =      ', num2str(L1)]);
    disp(['L2               =      ', num2str(L2)]);
    disp(['L3               =      ', num2str(L3)]);
    % 转到ENU
    lon0rad = userLLHPosition(1)*pi/180;
    lat0rad = userLLHPosition(2)*pi/180;
    R       = [-sin(lon0rad)         ,  cos(lon0rad)             , 0;
           -sin(lat0rad)*cos(lon0rad), -sin(lat0rad)*sin(lon0rad), cos(lat0rad);
            cos(lat0rad)*cos(lon0rad),  cos(lat0rad)*sin(lon0rad), sin(lat0rad)];
    if((paraAll.withClcShift + paraAll.withClcDraft) == 1)
       Rl       = [-sin(lon0rad)         ,  cos(lon0rad)             , 0           , 0;
               -sin(lat0rad)*cos(lon0rad), -sin(lat0rad)*sin(lon0rad), cos(lat0rad), 0;
                cos(lat0rad)*cos(lon0rad),  cos(lat0rad)*sin(lon0rad), sin(lat0rad), 0;
                                        0,                          0,            0, 1];
    else
       Rl       = [-sin(lon0rad)         ,  cos(lon0rad)             ,            0, 0, 0;
               -sin(lat0rad)*cos(lon0rad), -sin(lat0rad)*sin(lon0rad), cos(lat0rad), 0, 0;
                cos(lat0rad)*cos(lon0rad),  cos(lat0rad)*sin(lon0rad), sin(lat0rad), 0, 0;
                                        0,                          0,            0, 1, 0;
                                        0,                          0,            0, 0, 1];

    end
    V_ENU   = R * V;
    v1_ENU  = V_ENU(:,1);
    v2_ENU  = V_ENU(:,2);
    v3_ENU  = V_ENU(:,3);

    disp('三轴转到东北天方向上的结果')
    matrixSout(V_ENU);
    disp('东北天方向上的误差(1sigma):')
    disp(['sigmaEast        =      ', num2str(sqrt((V_ENU(1,1)*L1)^2 + (V_ENU(1,2)*L2)^2 +(V_ENU(1,3)*L3)^2))])
    disp(['sigmaNorth       =      ', num2str(sqrt((V_ENU(2,1)*L1)^2 + (V_ENU(2,2)*L2)^2 +(V_ENU(2,3)*L3)^2))])
    disp(['sigmaVeti        =      ', num2str(sqrt((V_ENU(3,1)*L1)^2 + (V_ENU(3,2)*L2)^2 +(V_ENU(3,3)*L3)^2))])
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

    CRLBMatrix_ENU = Rl*CRLBMatrix*Rl';

    errorStd                =      sqrt(var(errorXByIterrM(clcDraftIndex,:))+var(errorYByIterrM(clcDraftIndex,:))+var(errorZByIterrM(clcDraftIndex,:)));
    disp(['errorStd         =      ', num2str(errorStd)]);
    posCRLB                 =      sqrt(CRLBMatrix(1,1)+CRLBMatrix(2,2)+CRLBMatrix(3,3));
    disp(['posCRLB          =      ', num2str(posCRLB)]);
    errorXStd               =      sqrt(var(errorXByIterrM(clcDraftIndex,:)));
    disp(['errorXStd        =      ', num2str(errorXStd)]);
    errorXMean              =      mean(errorXByIterrM(clcDraftIndex,:));
    disp(['errorXMean       =      ', num2str(errorXMean)]);
    errorYStd               =      sqrt(var(errorYByIterrM(clcDraftIndex,:)));
    disp(['errorYStd        =      ', num2str(errorYStd)]);
    errorYMean              =      mean(errorYByIterrM(clcDraftIndex,:));
    disp(['errorYMean       =      ', num2str(errorYMean)]);
    errorZStd               =      sqrt(var(errorZByIterrM(clcDraftIndex,:)));
    disp(['errorZStd        =      ', num2str(errorZStd)]);
    errorZMean              =      mean(errorZByIterrM(clcDraftIndex,:));
    disp(['errorZMean       =      ', num2str(errorZMean)]);
    errorXCRLB              =      sqrt(CRLBMatrix(1,1));
    disp(['errorXCRLB       =      ', num2str(errorXCRLB)]);
    errorYCRLB              =      sqrt(CRLBMatrix(2,2));
    disp(['errorYCRLB       =      ', num2str(errorYCRLB)]);
    errorZCRLB              =      sqrt(CRLBMatrix(3,3));
    disp(['errorZCRLB       =      ', num2str(errorZCRLB)]);
    
    errorEastStd            =      sqrt(var(errorEastByIterrM(clcDraftIndex,:)));
    disp(['errorEastStd     =      ', num2str(errorEastStd)]);
    errorEastMean           =      mean(errorEastByIterrM(clcDraftIndex,:));
    disp(['errorEastMean    =      ', num2str(errorEastMean)]);
    errorNorthStd           =      sqrt(var(errorNorthByIterrM(clcDraftIndex,:)));
    disp(['errorNorthStd    =      ', num2str(errorNorthStd)]);
    errorNorthMean          =      mean(errorNorthByIterrM(clcDraftIndex,:));
    disp(['errorNorthMean   =      ', num2str(errorNorthMean)]);
    errorVetiStd            =      sqrt(var(errorVetiByIterrM(clcDraftIndex,:)));
    disp(['errorVetiStd     =      ', num2str(errorVetiStd)]);
    errorVetiMean           =      mean(errorVetiByIterrM(clcDraftIndex,:));
    disp(['errorVetiMean    =      ', num2str(errorVetiMean)]);
    errorEastCRLB           =      sqrt(CRLBMatrix_ENU(1,1));
    disp(['errorEastCRLB    =      ', num2str(errorEastCRLB)]);
    errorNorthCRLB          =      sqrt(CRLBMatrix_ENU(2,2));
    disp(['errorNorthCRLB   =      ', num2str(errorNorthCRLB)]);
    errorVetiCRLB           =      sqrt(CRLBMatrix_ENU(3,3));
    disp(['errorVetiCRLB    =      ', num2str(errorVetiCRLB)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error95List(1, clcDraftIndex)      = error95;
    error65List(1, clcDraftIndex)      = error65;
    errorHori95List(1, clcDraftIndex)  = errorHori95;
    errorVeti95JList(1, clcDraftIndex) = errorVetiJ95;

    posCRLBList(1, clcDraftIndex)      = posCRLB;
    errorStdList(1, clcDraftIndex)     = errorStd;
    
    %%%%%%%% 误差分布（ENU + EN）%%%%%%%%
    scale       = 1;
    lineWidth   = 3;
    maxHeadSize = 1;
    fontSize    = 20;
    pointSize   = 6;

    range1      = 800;
    range2      = 800;
    
    subplot(4,9, (clcDraftIndex-1)*3 + 1)
    scatter3(errorEastByIterrM(clcDraftIndex,:), errorNorthByIterrM(clcDraftIndex,:), errorVetiByIterrM(clcDraftIndex,:), pointSize, 'filled', 'r');
    xlim([-range1 range1])
    ylim([-range1 range1])
    zlim([-range1 range1])

    xlabel('东向误差/(m)', FontSize=fontSize);
    ylabel('北向误差/(m)', FontSize=fontSize);
    zlabel('天向误差/(m)', FontSize=fontSize);
    grid on;
    title('定位误差', FontSize=fontSize);

    subplot(4,9, (clcDraftIndex-1)*3 + 2)
    scatter(errorEastByIterrM(clcDraftIndex,:), errorNorthByIterrM(clcDraftIndex,:), pointSize, 'filled', 'r');
    xlabel('东向误差/(m)', FontSize=fontSize);
    ylabel('北向误差/(m)', FontSize=fontSize);
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-range1 range1]);
    ylim([-range1 range1]);
    grid on;
    title('水平定位误差', FontSize=fontSize);

    subplot(4,9, (clcDraftIndex-1)*3 + 3)
    scatter(errorNorthByIterrM(clcDraftIndex,:), errorVetiByIterrM(clcDraftIndex,:), pointSize, 'filled', 'r');
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
end
%%%%%%%%%%%%%%%%%%%%%%定位误差随观测间隔变化%%%%%%%%%%%%%%%%%%%%%%%
% clcDraftList = [0.03, 0.09, 0.15, 0.21,...
%                        0.27, 0.33, 0.39, 0.45,...
%                        0.51, 0.57, 0.63];
linewidth = 2;
figure;
    plot(clcDraftList, error95List,  'Linewidth', linewidth);
    hold on;
    plot(clcDraftList, error65List,  'Linewidth', linewidth);
    hold on;
    plot(clcDraftList, posCRLBList,  'Linewidth', linewidth);
    hold on;
    plot(clcDraftList, errorStdList, 'Linewidth', linewidth);
    grid on;
    set(gca, 'xTick', [0.03, 0.09, 0.15, 0.21,...
                       0.27, 0.33, 0.39, 0.45,...
                       0.51, 0.57, 0.63]);
    set(gca, 'XTickLabel',{'0.03', '0.09', '0.15', '0.21',...
                       '0.27', '0.33', '0.39', '0.45',...
                       '0.51', '0.57', '0.63'});
    xlabel('频偏/赫兹',FontSize=fontSize);
    ylabel('定位误差/米',FontSize=fontSize);
    title('定位误差随频偏变化图',FontSize=fontSize);
    legend('定位误差(95%)','定位误差(65%)','定位误差CRLB','定位误差均方差',FontSize=fontSize);

    
size(errorEastByIterrM)
colorList = ['r','b','g','y','c','m','k'];
figure;
    subplot(1,2,1)
    for clcDraftIndex = 1:length(colorList)
        scatter(errorEastByIterrM(clcDraftIndex,:), errorNorthByIterrM(clcDraftIndex,:), pointSize, 'filled', colorList(clcDraftIndex));
        hold on;
    end
    xlabel('东向误差/米', FontSize=fontSize);
    ylabel('北向误差/米', FontSize=fontSize);
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-800 200]);
    ylim([-300 100]);
    set(gca, 'xTick', -800:100:200);
    set(gca, 'yTick', -300:100:100);
    grid on;
    legend('频差0.03赫兹','频差0.09赫兹','频差0.15赫兹',...
        '频差0.21赫兹','频差0.27赫兹','频差0.33赫兹',...
        '频差0.39赫兹');
    title('水平定位误差', FontSize=fontSize);
    
    subplot(1,2,2)  
    for clcDraftIndex = 1:length(colorList)
        scatter(errorNorthByIterrM(clcDraftIndex,:), errorVetiByIterrM(clcDraftIndex,:), pointSize, 'filled', colorList(clcDraftIndex));
        hold on;
    end
    xlabel('北向误差/米', FontSize=fontSize);
    ylabel('天向误差/米', FontSize=fontSize);
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-400 200]);
    ylim([-200 1000]);
    set(gca, 'xTick', -400:100:200);
    set(gca, 'yTick', -200:100:1000);
    grid on;
    legend('频差0.03赫兹','频差0.09赫兹','频差0.15赫兹',...
        '频差0.21赫兹','频差0.27赫兹','频差0.33赫兹',...
        '频差0.39赫兹');
    title('北向天向定位误差', FontSize=fontSize);
toc