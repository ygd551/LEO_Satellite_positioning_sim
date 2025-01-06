% 函数功能：仿真分析采样间隔对低轨单星定位结果的影响   v
tic 
clc;clear;
format long g

addpath(genpath("../../coordinateTransformation"));
paraAll = simulation_sampling_Para();
%%%%%%%%%%%%%%%%%%%参数设置%%%%%%%%%%%%%%%%%%%%%%
% 用户位置
userLLHPosition  = paraAll.userLLHPosition;
userECEFPosition = llh2ecef(userLLHPosition);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%取出观测量%%%%%%%%%%%%%%%%%%%%
%%%%% LEO观测量%%%%%%
disp(paraAll.obsDataFilePath)
[ttListAll, satECEFPositionListAll, satECEFVelocityListAll,...
      rangeListAll, dopplerListAll] = getLEOObsData(paraAll.obsDataFilePath);
Ns = length(rangeListAll);
disp(['多历元卫星观测量的个数为：',num2str(Ns)]);


samplingList        = [1,2,4,8,16];
samplingLength      = length(samplingList);


errorByIterrM       = zeros(samplingLength, paraAll.iterrT);
errorHoriByIterrM   = zeros(samplingLength, paraAll.iterrT);
errorEastByIterrM   = zeros(samplingLength, paraAll.iterrT);
errorEastByIterrJM  = zeros(samplingLength, paraAll.iterrT);
errorNorthByIterrM  = zeros(samplingLength, paraAll.iterrT);
errorNorthByIterrJM = zeros(samplingLength, paraAll.iterrT);
errorVetiByIterrM   = zeros(samplingLength, paraAll.iterrT);
errorVetiByIterrJM  = zeros(samplingLength, paraAll.iterrT);
errorXByIterrM      = zeros(samplingLength, paraAll.iterrT);
errorYByIterrM      = zeros(samplingLength, paraAll.iterrT);
errorZByIterrM      = zeros(samplingLength, paraAll.iterrT);
PDOPByIterrM        = zeros(samplingLength, paraAll.iterrT);
HDOPByIterrM        = zeros(samplingLength, paraAll.iterrT);
VDOPByIterrM        = zeros(samplingLength, paraAll.iterrT);
B0ByIterrM          = zeros(samplingLength, paraAll.iterrT);

error95List         = zeros(1, samplingLength);
error65List         = zeros(1, samplingLength);
errorHori95List     = zeros(1, samplingLength);
errorHori65List     = zeros(1, samplingLength);
errorVeti95JList    = zeros(1, samplingLength);
errorVeti65JList    = zeros(1, samplingLength);
errorEastJ95List    = zeros(1, samplingLength);
errorEastJ65List    = zeros(1, samplingLength);
errorNorthJ95List   = zeros(1, samplingLength);
errorNorthJ65List   = zeros(1, samplingLength);
  
posCRLBList         = zeros(1, samplingLength);
errorStdList        = zeros(1, samplingLength);

figure;
for samplingIndex = 1:samplingLength
    sampling           = samplingList(samplingIndex);
    disp('===============================================')
    disp(['观测间隔为', num2str(sampling), '仿真开始']);
    %%%%%%%%%%%%%%%%%%得到观测值%%%%%%%%%%%%%%%%%%%%%
    satECEFPositionList = satECEFPositionListAll(1:sampling:Ns, :);
    satECEFVelocityList = satECEFVelocityListAll(1:sampling:Ns, :);
    rangeList           = rangeListAll(1:sampling:Ns);
    dopplerList         = dopplerListAll(1:sampling:Ns);
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
        % input:rangeList dopplerList
        prSatList = [];
        pdSatList = [];

        errorPrList = [];
        for ns = 1:length(rangeList)
            errorPr   = paraAll.pseudoRangeError*randn;
    %         errorPrList = [errorPrList, errorPr];
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
        errorByIterrM(samplingIndex, iterr)       = error;
        errorHoriByIterrM(samplingIndex, iterr)   = errorHorizontal;
        errorEastByIterrM(samplingIndex, iterr)   = errorEast;
        errorEastByIterrJM(samplingIndex, iterr)  = abs(errorEast);
        errorNorthByIterrM(samplingIndex, iterr)  = errorNorth;
        errorNorthByIterrJM(samplingIndex, iterr) = abs(errorNorth);
        errorVetiByIterrM(samplingIndex, iterr)   = errorVertical;
        errorVetiByIterrJM(samplingIndex, iterr)  = abs(errorVertical);
        errorXByIterrM(samplingIndex, iterr)      = calPos(1)-userECEFPosition(1);
        errorYByIterrM(samplingIndex, iterr)      = calPos(2)-userECEFPosition(2);
        errorZByIterrM(samplingIndex, iterr)      = calPos(3)-userECEFPosition(3);
        B0ByIterrM(samplingIndex, iterr)          = B00List(end);
        PDOPByIterrM(samplingIndex, iterr)        = PDOP;
        HDOPByIterrM(samplingIndex, iterr)        = HDOP;
        VDOPByIterrM(samplingIndex, iterr)        = VDOP;

        if(rem(iterr,200) == 0)
            disp(['第', num2str(iterr), '次仿真结束']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end 

    errorByIterrSort        = sort(errorByIterrM(samplingIndex,:));
    errorHoriByIterrSort    = sort(errorHoriByIterrM(samplingIndex,:));
    errorVetiByIterrJSort   = sort(errorVetiByIterrJM(samplingIndex,:));
    errorEastByIterrJSort   = sort(errorEastByIterrJM(samplingIndex,:));
    errorNorthByIterrJSort  = sort(errorNorthByIterrJM(samplingIndex,:));

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
    
    PDOP                    =      mean(PDOPByIterrM(samplingIndex,:));
    disp(['PDOP             =      ', num2str(PDOP)]);
    HDOP                    =      mean(HDOPByIterrM(samplingIndex,:));
    disp(['HDOP             =      ', num2str(PDOP)]);
    VDOP                    =      mean(VDOPByIterrM(samplingIndex,:));
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

    errorStd                =      sqrt(var(errorXByIterrM(samplingIndex,:))+var(errorYByIterrM(samplingIndex,:))+var(errorZByIterrM(samplingIndex,:)));
    disp(['errorStd         =      ', num2str(errorStd)]);
    posCRLB                 =      sqrt(CRLBMatrix(1,1)+CRLBMatrix(2,2)+CRLBMatrix(3,3));
    disp(['posCRLB          =      ', num2str(posCRLB)]);
    errorXStd               =      sqrt(var(errorXByIterrM(samplingIndex,:)));
    disp(['errorXStd        =      ', num2str(errorXStd)]);
    errorXMean              =      mean(errorXByIterrM(samplingIndex,:));
    disp(['errorXMean       =      ', num2str(errorXMean)]);
    errorYStd               =      sqrt(var(errorYByIterrM(samplingIndex,:)));
    disp(['errorYStd        =      ', num2str(errorYStd)]);
    errorYMean              =      mean(errorYByIterrM(samplingIndex,:));
    disp(['errorYMean       =      ', num2str(errorYMean)]);
    errorZStd               =      sqrt(var(errorZByIterrM(samplingIndex,:)));
    disp(['errorZStd        =      ', num2str(errorZStd)]);
    errorZMean              =      mean(errorZByIterrM(samplingIndex,:));
    disp(['errorZMean       =      ', num2str(errorZMean)]);
    errorXCRLB              =      sqrt(CRLBMatrix(1,1));
    disp(['errorXCRLB       =      ', num2str(errorXCRLB)]);
    errorYCRLB              =      sqrt(CRLBMatrix(2,2));
    disp(['errorYCRLB       =      ', num2str(errorYCRLB)]);
    errorZCRLB              =      sqrt(CRLBMatrix(3,3));
    disp(['errorZCRLB       =      ', num2str(errorZCRLB)]);
    
    errorEastStd            =      sqrt(var(errorEastByIterrM(samplingIndex,:)));
    disp(['errorEastStd     =      ', num2str(errorEastStd)]);
    errorEastMean           =      mean(errorEastByIterrM(samplingIndex,:));
    disp(['errorEastMean    =      ', num2str(errorEastMean)]);
    errorNorthStd           =      sqrt(var(errorNorthByIterrM(samplingIndex,:)));
    disp(['errorNorthStd    =      ', num2str(errorNorthStd)]);
    errorNorthMean          =      mean(errorNorthByIterrM(samplingIndex,:));
    disp(['errorNorthMean   =      ', num2str(errorNorthMean)]);
    errorVetiStd            =      sqrt(var(errorVetiByIterrM(samplingIndex,:)));
    disp(['errorVetiStd     =      ', num2str(errorVetiStd)]);
    errorVetiMean           =      mean(errorVetiByIterrM(samplingIndex,:));
    disp(['errorVetiMean    =      ', num2str(errorVetiMean)]);
    errorEastCRLB           =      sqrt(CRLBMatrix_ENU(1,1));
    disp(['errorEastCRLB    =      ', num2str(errorEastCRLB)]);
    errorNorthCRLB          =      sqrt(CRLBMatrix_ENU(2,2));
    disp(['errorNorthCRLB   =      ', num2str(errorNorthCRLB)]);
    errorVetiCRLB           =      sqrt(CRLBMatrix_ENU(3,3));
    disp(['errorVetiCRLB    =      ', num2str(errorVetiCRLB)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error95List(1, samplingIndex)      = error95;
    error65List(1, samplingIndex)      = error65;
    errorHori95List(1, samplingIndex)  = errorHori95;
    errorVeti95JList(1, samplingIndex) = errorVetiJ95;

    posCRLBList(1, samplingIndex)      = posCRLB;
    errorStdList(1, samplingIndex)     = errorStd;
    
    %%%%%%%% 误差分布（ENU + EN）%%%%%%%%
    scale       = 1;
    lineWidth   = 3;
    maxHeadSize = 1;
    fontSize    = 20;
    pointSize   = 6;

    range1      = 800;
    range2      = 800;
    
    subplot(3,6, (samplingIndex-1)*3 + 1)
    scatter3(errorEastByIterrM(samplingIndex,:), errorNorthByIterrM(samplingIndex,:), errorVetiByIterrM(samplingIndex,:), pointSize, 'filled', 'r');
    xlim([-range1 range1])
    ylim([-range1 range1])
    zlim([-range1 range1])

    xlabel('东向误差/(m)', FontSize=fontSize);
    ylabel('北向误差/(m)', FontSize=fontSize);
    zlabel('天向误差/(m)', FontSize=fontSize);
    grid on;
    title('定位误差', FontSize=fontSize);

    subplot(3,6, (samplingIndex-1)*3 + 2)
    scatter(errorEastByIterrM(samplingIndex,:), errorNorthByIterrM(samplingIndex,:), pointSize, 'filled', 'r');
    xlabel('东向误差/(m)', FontSize=fontSize);
    ylabel('北向误差/(m)', FontSize=fontSize);
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-range1 range1]);
    ylim([-range1 range1]);
    grid on;
    title('水平定位误差', FontSize=fontSize);

    subplot(3,6, (samplingIndex-1)*3 + 3)
    scatter(errorNorthByIterrM(samplingIndex,:), errorVetiByIterrM(samplingIndex,:), pointSize, 'filled', 'r');
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
linewidth = 2;
figure;
    plot(samplingList, error95List,  'Linewidth', linewidth);
    hold on;
    plot(samplingList, error65List,  'Linewidth', linewidth);
    hold on;
    plot(samplingList, posCRLBList,  'Linewidth', linewidth);
    hold on;
    plot(samplingList, errorStdList, 'Linewidth', linewidth);
    grid on;
    set(gca, 'xTick', [1,2,4,8,16]);
    set(gca, 'XTickLabel',{'1','2','4','8','16'});
    xlabel('观测间隔/秒',FontSize=fontSize);
    ylabel('定位误差/米',FontSize=fontSize);
    title('定位误差随观测间隔变化图',FontSize=fontSize);
    legend('定位误差(95%)','定位误差(65%)','定位误差CRLB','定位误差均方差',FontSize=fontSize);
    
colorList = ['r','b','g','y','c'];
figure;
    for samplingIndex = 1:samplingLength
        scatter(errorEastByIterrM(samplingIndex,:), errorNorthByIterrM(samplingIndex,:), pointSize, 'filled', colorList(samplingIndex));
        hold on;
        
        xlabel('东向误差/(m)', FontSize=fontSize);
        ylabel('北向误差/(m)', FontSize=fontSize);
        ax1 = gca;
        ax1.XAxisLocation = 'origin';
        ax1.YAxisLocation = 'origin';
        xlim([-range1 range1]);
        ylim([-range1 range1]);
        grid on;
        title('水平定位误差', FontSize=fontSize);
    end
figure;
    for samplingIndex = 1:samplingLength
        scatter(errorNorthByIterrM(samplingIndex,:), errorVetiByIterrM(samplingIndex,:), pointSize, 'filled', colorList(samplingIndex));
        hold on;
        
        xlabel('北向误差/(m)', FontSize=fontSize);
        ylabel('天向误差/(m)', FontSize=fontSize);
        ax1 = gca;
        ax1.XAxisLocation = 'origin';
        ax1.YAxisLocation = 'origin';
        xlim([-range1 range1]);
        ylim([-range1 range1]);
        grid on;
        title('北向天向定位误差', FontSize=fontSize);
    end
 figure;
    scatter(errorEastByIterrM(1,:), errorNorthByIterrM(1,:), pointSize, 'filled', colorList(samplingIndex));
    hold on;
    
    xlabel('东向误差/(m)', FontSize=fontSize);
    ylabel('北向误差/(m)', FontSize=fontSize);
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-range1 range1]);
    ylim([-range1 range1]);
    grid on;
    title('间隔1秒，水平定位误差', FontSize=fontSize);
  figure;
    scatter(errorEastByIterrM(3,:), errorNorthByIterrM(3,:), pointSize, 'filled', colorList(samplingIndex));
    hold on;
    
    xlabel('东向误差/(m)', FontSize=fontSize);
    ylabel('北向误差/(m)', FontSize=fontSize);
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-range1 range1]);
    ylim([-range1 range1]);
    grid on;
    title('间隔4秒，水平定位误差', FontSize=fontSize);
 figure;
    scatter(errorEastByIterrM(5,:), errorNorthByIterrM(5,:), pointSize, 'filled', colorList(samplingIndex));
    hold on;
    
    xlabel('东向误差/(m)', FontSize=fontSize);
    ylabel('北向误差/(m)', FontSize=fontSize);
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-range1 range1]);
    ylim([-range1 range1]);
    grid on;
    title('间隔16秒，水平定位误差', FontSize=fontSize);
  
toc