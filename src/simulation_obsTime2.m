% 函数功能：仿真观测时长对低轨单星多历元定位结果的影响  v
tic  
clc;clear;
format long g

addpath(genpath("../../coordinateTransformation"));
paraAll = simulation_obsTime2_para();
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

%%%%%%%%%%%%%%%%%%得到初始值%%%%%%%%%%%%%%%%%%%%%
if(paraAll.positionStartWay == 0)
    userECEFPositionStart = getInitialPos(rangeListAll, satECEFPositionListAll);
end
if(paraAll.positionStartWay == 1)
    userECEFPositionStart = [userECEFPosition(1) + paraAll.errorxyz(1), userECEFPosition(2) + paraAll.errorxyz(2), userECEFPosition(3) + paraAll.errorxyz(3)];
end
if(paraAll.positionStartWay == 2)
    userECEFPositionStart = [0, 0, 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obsTimeList         = 300:10:600;
obsTimeLength       = length(obsTimeList);

errorByIterrM       = zeros(obsTimeLength, paraAll.iterrT);
errorHoriByIterrM   = zeros(obsTimeLength, paraAll.iterrT);
errorEastByIterrM   = zeros(obsTimeLength, paraAll.iterrT);
errorEastByIterrJM  = zeros(obsTimeLength, paraAll.iterrT);
errorNorthByIterrM  = zeros(obsTimeLength, paraAll.iterrT);
errorNorthByIterrJM = zeros(obsTimeLength, paraAll.iterrT);
errorVetiByIterrM   = zeros(obsTimeLength, paraAll.iterrT);
errorVetiByIterrJM  = zeros(obsTimeLength, paraAll.iterrT);
errorXByIterrM      = zeros(obsTimeLength, paraAll.iterrT);
errorYByIterrM      = zeros(obsTimeLength, paraAll.iterrT);
errorZByIterrM      = zeros(obsTimeLength, paraAll.iterrT);
PDOPByIterrM        = zeros(obsTimeLength, paraAll.iterrT);
HDOPByIterrM        = zeros(obsTimeLength, paraAll.iterrT);
VDOPByIterrM        = zeros(obsTimeLength, paraAll.iterrT);
B0ByIterrM          = zeros(obsTimeLength, paraAll.iterrT);

error95List         = zeros(1, obsTimeLength);
error65List         = zeros(1, obsTimeLength);
errorHori95List     = zeros(1, obsTimeLength);
errorHori65List     = zeros(1, obsTimeLength);
errorVeti95JList    = zeros(1, obsTimeLength);
errorVeti65JList    = zeros(1, obsTimeLength);
errorEastJ95List    = zeros(1, obsTimeLength);
errorEastJ65List    = zeros(1, obsTimeLength);
errorNorthJ95List   = zeros(1, obsTimeLength);
errorNorthJ65List   = zeros(1, obsTimeLength);
  
posCRLBList         = zeros(1, obsTimeLength);
errorStdList        = zeros(1, obsTimeLength);

figure;
for obsTimeIndex = 1:obsTimeLength
    obsTime             = obsTimeList(obsTimeIndex);
    disp('===============================================')
    disp(['观测时长为', num2str(obsTime), '仿真开始']);
    %%%%%%%%%%%%%%%%%%得到观测值%%%%%%%%%%%%%%%%%%%%%
    satECEFPositionList = satECEFPositionListAll(1:obsTime, :);
    satECEFVelocityList = satECEFVelocityListAll(1:obsTime, :);
    rangeList           = rangeListAll(1:obsTime);
    dopplerList         = dopplerListAll(1:obsTime);
    
    for iterr = 1:paraAll.iterrT
        %%%%%%%%%%%%%%%% 加观测误差%%%%%%%%%%%%%%%%
        % input:rangeList dopplerList
        prSatList = [];
        pdSatList = [];
        for ns = 1:length(rangeList)
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
        errorByIterrM(obsTimeIndex, iterr)       = error;
        errorHoriByIterrM(obsTimeIndex, iterr)   = errorHorizontal;
        errorEastByIterrM(obsTimeIndex, iterr)   = errorEast;
        errorEastByIterrJM(obsTimeIndex, iterr)  = abs(errorEast);
        errorNorthByIterrM(obsTimeIndex, iterr)  = errorNorth;
        errorNorthByIterrJM(obsTimeIndex, iterr) = abs(errorNorth);
        errorVetiByIterrM(obsTimeIndex, iterr)   = errorVertical;
        errorVetiByIterrJM(obsTimeIndex, iterr)  = abs(errorVertical);
        errorXByIterrM(obsTimeIndex, iterr)      = calPos(1)-userECEFPosition(1);
        errorYByIterrM(obsTimeIndex, iterr)      = calPos(2)-userECEFPosition(2);
        errorZByIterrM(obsTimeIndex, iterr)      = calPos(3)-userECEFPosition(3);
        B0ByIterrM(obsTimeIndex, iterr)          = B00List(end);
        PDOPByIterrM(obsTimeIndex, iterr)        = PDOP;
        HDOPByIterrM(obsTimeIndex, iterr)        = HDOP;
        VDOPByIterrM(obsTimeIndex, iterr)        = VDOP;

        if(rem(iterr,200) == 0)
            disp(['第', num2str(iterr), '次仿真结束']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end 

    errorByIterrSort        = sort(errorByIterrM(obsTimeIndex,:));
    errorHoriByIterrSort    = sort(errorHoriByIterrM(obsTimeIndex,:));
    errorVetiByIterrJSort   = sort(errorVetiByIterrJM(obsTimeIndex,:));
    errorEastByIterrJSort   = sort(errorEastByIterrJM(obsTimeIndex,:));
    errorNorthByIterrJSort  = sort(errorNorthByIterrJM(obsTimeIndex,:));

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
    
    PDOP                    =      mean(PDOPByIterrM(obsTimeIndex,:));
    disp(['PDOP             =      ', num2str(PDOP)]);
    HDOP                    =      mean(HDOPByIterrM(obsTimeIndex,:));
    disp(['HDOP             =      ', num2str(PDOP)]);
    VDOP                    =      mean(VDOPByIterrM(obsTimeIndex,:));
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

    errorStd                =      sqrt(var(errorXByIterrM(obsTimeIndex,:))+var(errorYByIterrM(obsTimeIndex,:))+var(errorZByIterrM(obsTimeIndex,:)));
    disp(['errorStd         =      ', num2str(errorStd)]);
    posCRLB                 =      sqrt(CRLBMatrix(1,1)+CRLBMatrix(2,2)+CRLBMatrix(3,3));
    disp(['posCRLB          =      ', num2str(posCRLB)]);
    errorXStd               =      sqrt(var(errorXByIterrM(obsTimeIndex,:)));
    disp(['errorXStd        =      ', num2str(errorXStd)]);
    errorXMean              =      mean(errorXByIterrM(obsTimeIndex,:));
    disp(['errorXMean       =      ', num2str(errorXMean)]);
    errorYStd               =      sqrt(var(errorYByIterrM(obsTimeIndex,:)));
    disp(['errorYStd        =      ', num2str(errorYStd)]);
    errorYMean              =      mean(errorYByIterrM(obsTimeIndex,:));
    disp(['errorYMean       =      ', num2str(errorYMean)]);
    errorZStd               =      sqrt(var(errorZByIterrM(obsTimeIndex,:)));
    disp(['errorZStd        =      ', num2str(errorZStd)]);
    errorZMean              =      mean(errorZByIterrM(obsTimeIndex,:));
    disp(['errorZMean       =      ', num2str(errorZMean)]);
    errorXCRLB              =      sqrt(CRLBMatrix(1,1));
    disp(['errorXCRLB       =      ', num2str(errorXCRLB)]);
    errorYCRLB              =      sqrt(CRLBMatrix(2,2));
    disp(['errorYCRLB       =      ', num2str(errorYCRLB)]);
    errorZCRLB              =      sqrt(CRLBMatrix(3,3));
    disp(['errorZCRLB       =      ', num2str(errorZCRLB)]);
    
    errorEastStd            =      sqrt(var(errorEastByIterrM(obsTimeIndex,:)));
    disp(['errorEastStd     =      ', num2str(errorEastStd)]);
    errorEastMean           =      mean(errorEastByIterrM(obsTimeIndex,:));
    disp(['errorEastMean    =      ', num2str(errorEastMean)]);
    errorNorthStd           =      sqrt(var(errorNorthByIterrM(obsTimeIndex,:)));
    disp(['errorNorthStd    =      ', num2str(errorNorthStd)]);
    errorNorthMean          =      mean(errorNorthByIterrM(obsTimeIndex,:));
    disp(['errorNorthMean   =      ', num2str(errorNorthMean)]);
    errorVetiStd            =      sqrt(var(errorVetiByIterrM(obsTimeIndex,:)));
    disp(['errorVetiStd     =      ', num2str(errorVetiStd)]);
    errorVetiMean           =      mean(errorVetiByIterrM(obsTimeIndex,:));
    disp(['errorVetiMean    =      ', num2str(errorVetiMean)]);
    errorEastCRLB           =      sqrt(CRLBMatrix_ENU(1,1));
    disp(['errorEastCRLB    =      ', num2str(errorEastCRLB)]);
    errorNorthCRLB          =      sqrt(CRLBMatrix_ENU(2,2));
    disp(['errorNorthCRLB   =      ', num2str(errorNorthCRLB)]);
    errorVetiCRLB           =      sqrt(CRLBMatrix_ENU(3,3));
    disp(['errorVetiCRLB    =      ', num2str(errorVetiCRLB)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error95List(1, obsTimeIndex)      = error95;
    error65List(1, obsTimeIndex)      = error65;
    errorHori95List(1, obsTimeIndex)  = errorHori95;
    errorVeti95JList(1, obsTimeIndex) = errorVetiJ95;

    posCRLBList(1, obsTimeIndex)      = posCRLB;
    errorStdList(1, obsTimeIndex)     = errorStd;
    
    %%%%%%%% 误差分布（ENU + EN）%%%%%%%%
%     scale       = 1;
%     lineWidth   = 3;
%     maxHeadSize = 1;
%     fontSize    = 20;
%     pointSize   = 6;
% 
%     range1      = 800;
%     range2      = 800;
    
%     subplot(2,6, (obsTimeIndex-1)*3 + 1)
%     scatter3(errorEastByIterrM(obsTimeIndex,:), errorNorthByIterrM(obsTimeIndex,:), errorVetiByIterrM(obsTimeIndex,:), pointSize, 'filled', 'r');
%     xlim([-range1 range1])
%     ylim([-range1 range1])
%     zlim([-range1 range1])
% 
%     xlabel('东向误差/(m)', FontSize=fontSize);
%     ylabel('北向误差/(m)', FontSize=fontSize);
%     zlabel('天向误差/(m)', FontSize=fontSize);
%     grid on;
%     title('定位误差', FontSize=fontSize);
% 
%     subplot(2,6, (obsTimeIndex-1)*3 + 2)
%     scatter(errorEastByIterrM(obsTimeIndex,:), errorNorthByIterrM(obsTimeIndex,:), pointSize, 'filled', 'r');
%     xlabel('东向误差/(m)', FontSize=fontSize);
%     ylabel('北向误差/(m)', FontSize=fontSize);
%     ax1 = gca;
%     ax1.XAxisLocation = 'origin';
%     ax1.YAxisLocation = 'origin';
%     xlim([-range1 range1]);
%     ylim([-range1 range1]);
%     grid on;
%     title('水平定位误差', FontSize=fontSize);
% 
%     subplot(2,6, (obsTimeIndex-1)*3 + 3)
%     scatter(errorNorthByIterrM(obsTimeIndex,:), errorVetiByIterrM(obsTimeIndex,:), pointSize, 'filled', 'r');
%     xlabel('北向误差/(m)', FontSize=fontSize);
%     ylabel('天向误差/(m)', FontSize=fontSize);
%     ax1 = gca;
%     ax1.XAxisLocation = 'origin';
%     ax1.YAxisLocation = 'origin';
%     xlim([-range2 range2]);
%     ylim([-range2 range2]);
%     grid on;
%     title('北向天向误差', FontSize=fontSize);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

obsTimeList = 300:10:600;

% 有高度辅助的结果
% error95List = [24.1355  21.9031  21.0955  20.7067  19.3673  19.0904  18.4612  17.9244  18.2348  17.3032  16.4874  16.7504  16.7324  17.4136  15.9527  15.7317  16.6034  15.1543  15.4568  15.7336  16.0957  15.2828  15.5165  15.7353  14.853  15.0401  15.1469  14.8947  15.1251  15.024  15.4835  
% ];
% error65List = [12.9949  12.0709  11.3804  11.3547  10.1831  10.1211  9.7636  9.6177  9.3342  9.1956  8.8089  8.899  8.2817  8.4191  8.4995  8.0815  8.1675  8.078  7.9217  7.7528  8.0169  7.7435  7.8029  7.8368  7.635  7.7048  7.651  7.4254  7.8656  7.603  7.6893  
% ];
% posCRLBList = [12.72  12.0471  11.4694  10.9731  10.5463  10.1786  9.8616  9.5877  9.3508  9.1455  8.9672  8.8122  8.6771  8.5591  8.4561  8.3658  8.2867  8.2173  8.1564  8.1028  8.0558  8.0144  7.978  7.9461  7.9181  7.8936  7.8721  7.8533  7.8369  7.8226  7.8101  
% ];
% errorStdList = [12.7718  11.8359  11.4621  11.1078  10.4588  10.141  9.928  9.6687  9.4918  9.1196  8.8595  8.8732  8.625  8.8475  8.6154  8.4021  8.4049  8.0898  8.1877  8.2125  8.3325  7.9855  8.0719  8.1036  7.7497  7.8482  7.8974  7.7524  7.9117  7.7971  7.9356  
% ];

linewidth = 2;
figure;
    plot(obsTimeList, error95List, 'Linewidth', linewidth);
    hold on;
    plot(obsTimeList, error65List, 'Linewidth', linewidth);
    hold on;
    plot(obsTimeList, posCRLBList, 'Linewidth', linewidth);
    hold on;
    plot(obsTimeList, errorStdList, 'Linewidth', linewidth);
    grid on;
    xlabel('观测时长/秒');
    ylabel('定位误差/米');
    title('定位误差随观测时长变化图');
    legend('定位误差(95%)','定位误差(65%)','定位误差CRLB','定位误差均方差');
% colorList = ['r','b','g','y'];
% figure;
%     for obsTimeIndex = 1:obsTimeLength
%         scatter(errorEastByIterrM(obsTimeIndex,:), errorNorthByIterrM(obsTimeIndex,:), pointSize, 'filled', colorList(obsTimeIndex));
%         hold on;
%         
%         xlabel('东向误差/(m)', FontSize=fontSize);
%         ylabel('北向误差/(m)', FontSize=fontSize);
%         ax1 = gca;
%         ax1.XAxisLocation = 'origin';
%         ax1.YAxisLocation = 'origin';
%         xlim([-range1 range1]);
%         ylim([-range1 range1]);
%         grid on;
%         title('水平定位误差', FontSize=fontSize);
%     end
% figure;
%     for obsTimeIndex = 1:obsTimeLength
%         scatter(errorNorthByIterrM(obsTimeIndex,:), errorVetiByIterrM(obsTimeIndex,:), pointSize, 'filled', colorList(obsTimeIndex));
%         hold on;
%         
%         xlabel('北向误差/(m)', FontSize=fontSize);
%         ylabel('天向误差/(m)', FontSize=fontSize);
%         ax1 = gca;
%         ax1.XAxisLocation = 'origin';
%         ax1.YAxisLocation = 'origin';
%         xlim([-range1 range1]);
%         ylim([-range1 range1]);
%         grid on;
%         title('北向天向定位误差', FontSize=fontSize);
%     end
toc