% 函数功能：讨论接收机初始钟差对定位结果的影响 v
tic
clc;clear;
format long g

addpath(genpath("../../coordinateTransformation"));
paraAll = simulation_clc0_para();
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

clc0List            = [ 1, 2, 3,...
                        4, 5, 6, 7,...
                        8, 9, 10];
clcDraft            = 0.03;
           
clc0Length          = length(clc0List);

errorByIterrM       = zeros(clc0Length, paraAll.iterrT);
errorHoriByIterrM   = zeros(clc0Length, paraAll.iterrT);
errorEastByIterrM   = zeros(clc0Length, paraAll.iterrT);
errorEastByIterrJM  = zeros(clc0Length, paraAll.iterrT);
errorNorthByIterrM  = zeros(clc0Length, paraAll.iterrT);
errorNorthByIterrJM = zeros(clc0Length, paraAll.iterrT);
errorVetiByIterrM   = zeros(clc0Length, paraAll.iterrT);
errorVetiByIterrJM  = zeros(clc0Length, paraAll.iterrT);
errorXByIterrM      = zeros(clc0Length, paraAll.iterrT);
errorYByIterrM      = zeros(clc0Length, paraAll.iterrT);
errorZByIterrM      = zeros(clc0Length, paraAll.iterrT);
PDOPByIterrM        = zeros(clc0Length, paraAll.iterrT);
HDOPByIterrM        = zeros(clc0Length, paraAll.iterrT);
VDOPByIterrM        = zeros(clc0Length, paraAll.iterrT);
B0ByIterrM          = zeros(clc0Length, paraAll.iterrT);

error95List         = zeros(1, clc0Length);
error65List         = zeros(1, clc0Length);
errorHori95List     = zeros(1, clc0Length);
errorHori65List     = zeros(1, clc0Length);
errorVeti95JList    = zeros(1, clc0Length);
errorVeti65JList    = zeros(1, clc0Length);
errorEastJ95List    = zeros(1, clc0Length);
errorEastJ65List    = zeros(1, clc0Length);
errorNorthJ95List   = zeros(1, clc0Length);
errorNorthJ65List   = zeros(1, clc0Length);
  
posCRLBList         = zeros(1, clc0Length);
errorStdList        = zeros(1, clc0Length);


for clc0Index = 1:clc0Length
    clc0            = clc0List(clc0Index);
    disp('===============================================')
    disp(['clc0为', num2str(clc0), '仿真开始']);
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
            prSat     = rangeList(ns) + (clc0+ttList(ns)*clcDraft) + errorPr;
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
        errorByIterrM(clc0Index, iterr)       = error;
        errorHoriByIterrM(clc0Index, iterr)   = errorHorizontal;
        errorEastByIterrM(clc0Index, iterr)   = errorEast;
        errorEastByIterrJM(clc0Index, iterr)  = abs(errorEast);
        errorNorthByIterrM(clc0Index, iterr)  = errorNorth;
        errorNorthByIterrJM(clc0Index, iterr) = abs(errorNorth);
        errorVetiByIterrM(clc0Index, iterr)   = errorVertical;
        errorVetiByIterrJM(clc0Index, iterr)  = abs(errorVertical);
        errorXByIterrM(clc0Index, iterr)      = calPos(1)-userECEFPosition(1);
        errorYByIterrM(clc0Index, iterr)      = calPos(2)-userECEFPosition(2);
        errorZByIterrM(clc0Index, iterr)      = calPos(3)-userECEFPosition(3);
        B0ByIterrM(clc0Index, iterr)          = B00List(end);
        PDOPByIterrM(clc0Index, iterr)        = PDOP;
        HDOPByIterrM(clc0Index, iterr)        = HDOP;
        VDOPByIterrM(clc0Index, iterr)        = VDOP;

        if(rem(iterr,200) == 0)
            disp(['第', num2str(iterr), '次仿真结束']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end 

    errorByIterrSort        = sort(errorByIterrM(clc0Index,:));
    errorHoriByIterrSort    = sort(errorHoriByIterrM(clc0Index,:));
    errorVetiByIterrJSort   = sort(errorVetiByIterrJM(clc0Index,:));
    errorEastByIterrJSort   = sort(errorEastByIterrJM(clc0Index,:));
    errorNorthByIterrJSort  = sort(errorNorthByIterrJM(clc0Index,:));

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
    
    PDOP                    =      mean(PDOPByIterrM(clc0Index,:));
    disp(['PDOP             =      ', num2str(PDOP)]);
    HDOP                    =      mean(HDOPByIterrM(clc0Index,:));
    disp(['HDOP             =      ', num2str(PDOP)]);
    VDOP                    =      mean(VDOPByIterrM(clc0Index,:));
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

    errorStd                =      sqrt(var(errorXByIterrM(clc0Index,:))+var(errorYByIterrM(clc0Index,:))+var(errorZByIterrM(clc0Index,:)));
    disp(['errorStd         =      ', num2str(errorStd)]);
    posCRLB                 =      sqrt(CRLBMatrix(1,1)+CRLBMatrix(2,2)+CRLBMatrix(3,3));
    disp(['posCRLB          =      ', num2str(posCRLB)]);
    errorXStd               =      sqrt(var(errorXByIterrM(clc0Index,:)));
    disp(['errorXStd        =      ', num2str(errorXStd)]);
    errorXMean              =      mean(errorXByIterrM(clc0Index,:));
    disp(['errorXMean       =      ', num2str(errorXMean)]);
    errorYStd               =      sqrt(var(errorYByIterrM(clc0Index,:)));
    disp(['errorYStd        =      ', num2str(errorYStd)]);
    errorYMean              =      mean(errorYByIterrM(clc0Index,:));
    disp(['errorYMean       =      ', num2str(errorYMean)]);
    errorZStd               =      sqrt(var(errorZByIterrM(clc0Index,:)));
    disp(['errorZStd        =      ', num2str(errorZStd)]);
    errorZMean              =      mean(errorZByIterrM(clc0Index,:));
    disp(['errorZMean       =      ', num2str(errorZMean)]);
    errorXCRLB              =      sqrt(CRLBMatrix(1,1));
    disp(['errorXCRLB       =      ', num2str(errorXCRLB)]);
    errorYCRLB              =      sqrt(CRLBMatrix(2,2));
    disp(['errorYCRLB       =      ', num2str(errorYCRLB)]);
    errorZCRLB              =      sqrt(CRLBMatrix(3,3));
    disp(['errorZCRLB       =      ', num2str(errorZCRLB)]);
    
    errorEastStd            =      sqrt(var(errorEastByIterrM(clc0Index,:)));
    disp(['errorEastStd     =      ', num2str(errorEastStd)]);
    errorEastMean           =      mean(errorEastByIterrM(clc0Index,:));
    disp(['errorEastMean    =      ', num2str(errorEastMean)]);
    errorNorthStd           =      sqrt(var(errorNorthByIterrM(clc0Index,:)));
    disp(['errorNorthStd    =      ', num2str(errorNorthStd)]);
    errorNorthMean          =      mean(errorNorthByIterrM(clc0Index,:));
    disp(['errorNorthMean   =      ', num2str(errorNorthMean)]);
    errorVetiStd            =      sqrt(var(errorVetiByIterrM(clc0Index,:)));
    disp(['errorVetiStd     =      ', num2str(errorVetiStd)]);
    errorVetiMean           =      mean(errorVetiByIterrM(clc0Index,:));
    disp(['errorVetiMean    =      ', num2str(errorVetiMean)]);
    errorEastCRLB           =      sqrt(CRLBMatrix_ENU(1,1));
    disp(['errorEastCRLB    =      ', num2str(errorEastCRLB)]);
    errorNorthCRLB          =      sqrt(CRLBMatrix_ENU(2,2));
    disp(['errorNorthCRLB   =      ', num2str(errorNorthCRLB)]);
    errorVetiCRLB           =      sqrt(CRLBMatrix_ENU(3,3));
    disp(['errorVetiCRLB    =      ', num2str(errorVetiCRLB)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error95List(1, clc0Index)      = error95;
    error65List(1, clc0Index)      = error65;
    errorHori95List(1, clc0Index)  = errorHori95;
    errorVeti95JList(1, clc0Index) = errorVetiJ95;

    posCRLBList(1, clc0Index)      = posCRLB;
    errorStdList(1, clc0Index)     = errorStd;
end
%%%%%%%%%%%%%%%%%%%%%%定位误差随观测间隔变化%%%%%%%%%%%%%%%%%%%%%%%

figure;
    plot(clc0List, error95List);
    hold on;
    plot(clc0List, error65List);
    hold on;
    plot(clc0List, posCRLBList);
    hold on;
    plot(clc0List, errorStdList);
    grid on;
    set(gca, 'xTick', clc0List);
    set(gca, 'XTickLabel',{'1', '2', '3', '4',...
                       '5', '6', '7', '8',...
                       '9', '10'});
    xlabel('初始钟差/米');
    ylabel('定位误差/米');
    title('定位误差随初始钟差变化图');
    legend('定位误差(95%)','定位误差(65%)','定位误差CRLB','定位误差均方差');
 
pointSize = 6;
colorList = ['r','b','g','y'];
figure;
    subplot(1,2,1)
    for clc0Index = 1:length(colorList)
        scatter(errorEastByIterrM(clc0Index,:), errorNorthByIterrM(clc0Index,:), pointSize, 'filled', colorList(clc0Index));
        hold on;
    end
    xlabel('东向误差/米');
    ylabel('北向误差/米');
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-500 500]);
    ylim([-200 200]);
    set(gca, 'xTick', -500:100:500);
    set(gca, 'yTick', -200:100:200);
    grid on;
    legend('初始钟差1米','初始钟差2米','初始钟差3米',...
        '初始钟差4米');
    title('水平定位误差');
    
    subplot(1,2,2)  
    for clc0Index = 1:length(colorList)
        scatter(errorNorthByIterrM(clc0Index,:), errorVetiByIterrM(clc0Index,:), pointSize, 'filled', colorList(clc0Index));
        hold on;
    end
    xlabel('北向误差/米');
    ylabel('天向误差/米');
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-200 200]);
    ylim([-400 600]);
    set(gca, 'xTick', -200:100:200);
    set(gca, 'yTick', -400:100:600);
    grid on;
    title('北向天向定位误差');
toc