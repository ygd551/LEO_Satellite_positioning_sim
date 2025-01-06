% 函数功能：钟差信息辅助伪距定位 v
tic
clc;clear;
format long g

Time = clock; 
disp([num2str(Time(1)),'/',num2str(Time(2)),'/',num2str(Time(3)),' ',num2str(Time(4)),':',num2str(Time(5))])

addpath(genpath("../../coordinateTransformation"));
paraAll = simulation_clcAssist_para();
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
%%%%% 接收机钟差观测量%%%%%%
B0      = 7.5611;
B1      = paraAll.B1;
clcList = zeros(1, Ns);
t1      = ttList(1);
for i = 1:Ns
    clcListAll(i) = B0 + B1 * (ttList(i)-t1) + randn * paraAll.clcError;
end

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
lh               = 2;
errorByIterr       = zeros(lh, paraAll.iterrT);
errorHoriByIterr   = zeros(lh, paraAll.iterrT);
errorEastByIterr   = zeros(lh, paraAll.iterrT);
errorEastByIterrJ  = zeros(lh, paraAll.iterrT);
errorNorthByIterr  = zeros(lh, paraAll.iterrT);
errorNorthByIterrJ = zeros(lh, paraAll.iterrT);
errorVetiByIterr   = zeros(lh, paraAll.iterrT);
errorVetiByIterrJ  = zeros(lh, paraAll.iterrT);
errorXByIterr      = zeros(lh, paraAll.iterrT);
errorYByIterr      = zeros(lh, paraAll.iterrT);
errorZByIterr      = zeros(lh, paraAll.iterrT);
PDOPByIterr        = zeros(lh, paraAll.iterrT);
HDOPByIterr        = zeros(lh, paraAll.iterrT);
VDOPByIterr        = zeros(lh, paraAll.iterrT);
B0ByIterr          = zeros(lh, paraAll.iterrT);
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
    LSM601(paraAll.f0,...
    prSatList, pdSatList, satECEFPositionList, satECEFVelocityList,... % 卫星观测量
    [], [], [],...                 % VOR/DME观测量
    height,...  % 高度观测量
    clcList, ttList,...
    paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,paraAll.clcError,...% 观测误差
    paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
    paraAll.withVORFlag, paraAll.withDMEFlag, paraAll.withHeightFlag, paraAll.weightFlag,...
    0,...
    userECEFPositionStart, 0, 0, paraAll.clcOrder); 

    calPos            = [x0List(end),y0List(end),z0List(end)];
    error             = norm(calPos-userECEFPosition);
    [errorEast, errorNorth, errorHorizontal, errorVertical]=...
    errorECEF2ENU(calPos, userECEFPosition);
    errorByIterr(1, iterr)       = error;
    errorHoriByIterr(1, iterr)   = errorHorizontal;
    errorEastByIterr(1, iterr)   = errorEast;
    errorEastByIterrJ(1, iterr)  = abs(errorEast);
    errorNorthByIterr(1, iterr)  = errorNorth;
    errorNorthByIterrJ(1, iterr) = abs(errorNorth);
    errorVetiByIterr(1, iterr)   = errorVertical;
    errorVetiByIterrJ(1, iterr)  = abs(errorVertical);
    errorXByIterr(1, iterr)      = calPos(1)-userECEFPosition(1);
    errorYByIterr(1, iterr)      = calPos(2)-userECEFPosition(2);
    errorZByIterr(1, iterr)      = calPos(3)-userECEFPosition(3);
    B0ByIterr(1, iterr)          = B00List(end);


    %%%%%%%%%%%%%%%%LSM  clc%%%%%%%%%%%%%%%%%%%%%
%     [x0List, y0List, z0List, B00List, B10List,...
%     PDOP, HDOP, VDOP, iterTime] = ...
%     LSM601(paraAll.f0,...
%     prSatList, pdSatList, satECEFPositionList, satECEFVelocityList,... % 卫星观测量
%     [], [], [],...                 % VOR/DME观测量
%     height,...  % 高度观测量
%     clcList, ttList,...
%     paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,paraAll.clcError,...% 观测误差
%     paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
%     paraAll.withVORFlag, paraAll.withDMEFlag, paraAll.withHeightFlag, paraAll.weightFlag,...
%     1,...
%     userECEFPositionStart, 0, 0, paraAll.clcOrder); 


    [x0List, y0List, z0List, B00List, B10List,...
    PDOP, HDOP, VDOP, iterTime] = ...
    LSM602(paraAll.f0,...
    prSatList, pdSatList, satECEFPositionList, satECEFVelocityList,... % 卫星观测量
    [], [], [],...                 % VOR/DME观测量
    height,...  % 高度观测量
    clcList,...
    paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
    paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
    paraAll.withVORFlag, paraAll.withDMEFlag, paraAll.withHeightFlag, paraAll.weightFlag,...
    1,...
    userECEFPositionStart, 0, 0, paraAll.clcOrder); 

    calPos            = [x0List(end),y0List(end),z0List(end)];
    error             = norm(calPos-userECEFPosition);
    [errorEast, errorNorth, errorHorizontal, errorVertical]=...
    errorECEF2ENU(calPos, userECEFPosition);
    errorByIterr(2, iterr)       = error;
    errorHoriByIterr(2, iterr)   = errorHorizontal;
    errorEastByIterr(2, iterr)   = errorEast;
    errorEastByIterrJ(2, iterr)  = abs(errorEast);
    errorNorthByIterr(2, iterr)  = errorNorth;
    errorNorthByIterrJ(2, iterr) = abs(errorNorth);
    errorVetiByIterr(2, iterr)   = errorVertical;
    errorVetiByIterrJ(2, iterr)  = abs(errorVertical);
    errorXByIterr(2, iterr)      = calPos(1)-userECEFPosition(1);
    errorYByIterr(2, iterr)      = calPos(2)-userECEFPosition(2);
    errorZByIterr(2, iterr)      = calPos(3)-userECEFPosition(3);
    B0ByIterr(2, iterr)          = B00List(end);

    if(rem(iterr,50) == 0)
        disp(['第', num2str(iterr), '次仿真结束']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
errorByIterrSort1       = sort(errorByIterr(1,:));
errorHoriByIterrSort1   = sort(errorHoriByIterr(1,:));
errorVetiByIterrJSort1  = sort(errorVetiByIterrJ(1,:));
errorEastByIterrJSort1  = sort(errorEastByIterrJ(1,:));
errorNorthByIterrJSort1 = sort(errorNorthByIterrJ(1,:));

errorByIterrSort2       = sort(errorByIterr(2,:));
errorHoriByIterrSort2   = sort(errorHoriByIterr(2,:));
errorVetiByIterrJSort2  = sort(errorVetiByIterrJ(2,:));
errorEastByIterrJSort2  = sort(errorEastByIterrJ(2,:));
errorNorthByIterrJSort2 = sort(errorNorthByIterrJ(2,:));

%%%%%%%% 95%误差+DOP %%%%%%%%
disp('无钟差：')
disp(['error(95%)       = ', num2str(      errorByIterrSort1(paraAll.iterrT*0.95))]);
disp(['errorHori(95%)   = ', num2str(  errorHoriByIterrSort1(paraAll.iterrT*0.95))]);    
disp(['errorVetiJ(95%)  = ', num2str( errorVetiByIterrJSort1(paraAll.iterrT*0.95))]); 
disp(['errorEastJ(95%)  = ', num2str( errorEastByIterrJSort1(paraAll.iterrT*0.95))]); 
disp(['errorNorthJ(95%) = ', num2str(errorNorthByIterrJSort1(paraAll.iterrT*0.95))]); 
disp(['error(65%)       = ', num2str(      errorByIterrSort1(paraAll.iterrT*0.65))]);
disp(['errorHori(65%)   = ', num2str(  errorHoriByIterrSort1(paraAll.iterrT*0.65))]);    
disp(['errorVetiJ(65%)  = ', num2str( errorVetiByIterrJSort1(paraAll.iterrT*0.65))]); 
disp(['errorEastJ(65%)  = ', num2str( errorEastByIterrJSort1(paraAll.iterrT*0.65))]); 
disp(['errorNorthJ(65%) = ', num2str(errorNorthByIterrJSort1(paraAll.iterrT*0.65))]); 
disp('  ')
disp('有钟差：')
disp(['error(95%)       = ', num2str(      errorByIterrSort2(paraAll.iterrT*0.95))]);
disp(['errorHori(95%)   = ', num2str(  errorHoriByIterrSort2(paraAll.iterrT*0.95))]);    
disp(['errorVetiJ(95%)  = ', num2str( errorVetiByIterrJSort2(paraAll.iterrT*0.95))]); 
disp(['errorEastJ(95%)  = ', num2str( errorEastByIterrJSort2(paraAll.iterrT*0.95))]); 
disp(['errorNorthJ(95%) = ', num2str(errorNorthByIterrJSort2(paraAll.iterrT*0.95))]); 
disp(['error(65%)       = ', num2str(      errorByIterrSort2(paraAll.iterrT*0.65))]);
disp(['errorHori(65%)   = ', num2str(  errorHoriByIterrSort2(paraAll.iterrT*0.65))]);    
disp(['errorVetiJ(65%)  = ', num2str( errorVetiByIterrJSort2(paraAll.iterrT*0.65))]); 
disp(['errorEastJ(65%)  = ', num2str( errorEastByIterrJSort2(paraAll.iterrT*0.65))]); 
disp(['errorNorthJ(65%) = ', num2str(errorNorthByIterrJSort2(paraAll.iterrT*0.65))]); 

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
disp(['L1           = ', num2str(L1)]);
disp(['L2           = ', num2str(L2)]);
disp(['L3           = ', num2str(L3)]);
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
disp(['sigmaEast        = ', num2str(sqrt((V_ENU(1,1)*L1)^2 + (V_ENU(1,2)*L2)^2 +(V_ENU(1,3)*L3)^2))])
disp(['sigmaNorth       = ', num2str(sqrt((V_ENU(2,1)*L1)^2 + (V_ENU(2,2)*L2)^2 +(V_ENU(2,3)*L3)^2))])
disp(['sigmaVeti        = ', num2str(sqrt((V_ENU(3,1)*L1)^2 + (V_ENU(3,2)*L2)^2 +(V_ENU(3,3)*L3)^2))])
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
 
posCRLB        = sqrt(CRLBMatrix(1,1)+CRLBMatrix(2,2)+CRLBMatrix(3,3));

disp('无钟差：')
disp(['error_std        =      ', num2str(sqrt(var(errorXByIterr(1,:))+var(errorYByIterr(1,:))+var(errorZByIterr(1,:))))]);
disp(['posCRLB          =      ', num2str(posCRLB)]);
disp(['errorX_std       =      ', num2str(sqrt(var(errorXByIterr(1,:))))]);
disp(['errorX_mean      =      ', num2str(mean(errorXByIterr(1,:)))]);
disp(['errorY_std       =      ', num2str(sqrt(var(errorYByIterr(1,:))))]);
disp(['errorY_mean      =      ', num2str(mean(errorYByIterr(1,:)))]);
disp(['errorZ_std       =      ', num2str(sqrt(var(errorZByIterr(1,:))))]);
disp(['errorZ_mean      =      ', num2str(mean(errorZByIterr(1,:)))]);
disp(['errorX_CRLB      =      ', num2str(sqrt(CRLBMatrix(1,1)))]);
disp(['errorY_CRLB      =      ', num2str(sqrt(CRLBMatrix(2,2)))]);
disp(['errorZ_CRLB      =      ', num2str(sqrt(CRLBMatrix(3,3)))]);

disp(['errorEast_std    =      ', num2str(sqrt(var(errorEastByIterr(1,:)  )))]);
disp(['errorEast_mean   =      ', num2str(mean(errorEastByIterr(1,:)    ))]);
disp(['errorNorth_std   =      ', num2str(sqrt(var(errorNorthByIterr(1,:) )))]);
disp(['errorNorth_mean  =      ', num2str(mean(errorNorthByIterr(1,:)    ))]);
disp(['errorVeti_std    =      ', num2str(sqrt(var(errorVetiByIterr(1,:)  )))]);
disp(['errorVeti_mean   =      ', num2str(mean(errorVetiByIterr(1,:)  ))]);
disp(['errorEast_CRLB   =      ', num2str(sqrt(CRLBMatrix_ENU(1,1)))]);
disp(['errorNorthY_CRLB =      ', num2str(sqrt(CRLBMatrix_ENU(2,2)))]);
disp(['errorVeti_CRLB   =      ', num2str(sqrt(CRLBMatrix_ENU(3,3)))]);
disp(' ')
disp('有钟差：')
disp(['error_std        =      ', num2str(sqrt(var(errorXByIterr(2,:))+var(errorYByIterr(2,:))+var(errorZByIterr(2,:))))]);
disp(['posCRLB          =      ', num2str(posCRLB)]);
disp(['errorX_std       =      ', num2str(sqrt(var(errorXByIterr(2,:))))]);
disp(['errorX_mean      =      ', num2str(mean(errorXByIterr(2,:)))]);
disp(['errorY_std       =      ', num2str(sqrt(var(errorYByIterr(2,:))))]);
disp(['errorY_mean      =      ', num2str(mean(errorYByIterr(2,:)))]);
disp(['errorZ_std       =      ', num2str(sqrt(var(errorZByIterr(2,:))))]);
disp(['errorZ_mean      =      ', num2str(mean(errorZByIterr(2,:)))]);
disp(['errorX_CRLB      =      ', num2str(sqrt(CRLBMatrix(1,1)))]);
disp(['errorY_CRLB      =      ', num2str(sqrt(CRLBMatrix(2,2)))]);
disp(['errorZ_CRLB      =      ', num2str(sqrt(CRLBMatrix(3,3)))]);

disp(['errorEast_std    =      ', num2str(sqrt(var(errorEastByIterr(2,:)  )))]);
disp(['errorEast_mean   =      ', num2str(mean(errorEastByIterr(2,:)    ))]);
disp(['errorNorth_std   =      ', num2str(sqrt(var(errorNorthByIterr(2,:) )))]);
disp(['errorNorth_mean  =      ', num2str(mean(errorNorthByIterr(2,:)    ))]);
disp(['errorVeti_std    =      ', num2str(sqrt(var(errorVetiByIterr(2,:)  )))]);
disp(['errorVeti_mean   =      ', num2str(mean(errorVetiByIterr(2,:)  ))]);
disp(['errorEast_CRLB   =      ', num2str(sqrt(CRLBMatrix_ENU(1,1)))]);
disp(['errorNorthY_CRLB =      ', num2str(sqrt(CRLBMatrix_ENU(2,2)))]);
disp(['errorVeti_CRLB   =      ', num2str(sqrt(CRLBMatrix_ENU(3,3)))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% 误差分布（ENU + EN）%%%%%%%%
scale       = 1;
lineWidth   = 3;
maxHeadSize = 1;
fontSize    = 20;
pointSize   = 6;

range1      = 400;
range2      = 400;
figure(2)
    subplot(1,3,1)
    scatter3(errorEastByIterr(1,:), errorNorthByIterr(1,:), errorVetiByIterr(1,:), pointSize, 'filled', 'r');
    hold on;
    scatter3(errorEastByIterr(2,:), errorNorthByIterr(2,:), errorVetiByIterr(2,:), pointSize, 'filled', 'b');
    legend('无钟差辅助','有钟差辅助');
    
    xlim([-range1 range1])
    ylim([-range1 range1])
    zlim([-range1 range1])

    xlabel('东向误差/米', FontSize=fontSize);
    ylabel('北向误差/米', FontSize=fontSize);
    zlabel('天向误差/米', FontSize=fontSize);
    grid on;
    title('定位误差', FontSize=fontSize);
    %%%%%%%%%%%%%%%%%%%%%%
	subplot(1,3,2)  
    scatter(errorEastByIterr(1,:), errorNorthByIterr(1,:), pointSize, 'filled', 'r');
    hold on;
    scatter(errorEastByIterr(2,:), errorNorthByIterr(2,:), pointSize, 'filled', 'b');
    legend('无钟差辅助','有钟差辅助');
    
    xlabel('东向误差/米', FontSize=fontSize);
    ylabel('北向误差/米', FontSize=fontSize);
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-range1 range1]);
    ylim([-range1 range1]);
    grid on;
    title('水平定位误差', FontSize=fontSize);
    %%%%%%%%%%%%%%%%%%%%%%
    subplot(1,3,3)
    scatter(errorNorthByIterr(1,:), errorVetiByIterr(1,:), pointSize, 'filled', 'r');
    hold on;
    scatter(errorNorthByIterr(2,:), errorVetiByIterr(2,:), pointSize, 'filled', 'b');
    legend('无钟差辅助','有钟差辅助');
    
    xlabel('北向误差/米', FontSize=fontSize);
    ylabel('天向误差/米', FontSize=fontSize);
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-range2 range2]);
    ylim([-range2 range2]);
    grid on;
    title('北向天向误差', FontSize=fontSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
   
	subplot(1,2,1)  
    scatter(errorEastByIterr(1,:), errorNorthByIterr(1,:), pointSize, 'filled', 'r');
    hold on;
    scatter(errorEastByIterr(2,:), errorNorthByIterr(2,:), pointSize, 'filled', 'b');
    legend('无钟差辅助','有钟差辅助');
    
    xlabel('东向误差/(m)', FontSize=fontSize);
    ylabel('北向误差/(m)', FontSize=fontSize);
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-200 1000]);
    ylim([-200 200]);
    grid on;
    title('水平定位误差', FontSize=fontSize);
    %%%%%%%%%%%%%%%%%%%%%%
    subplot(1,2,2)
    scatter(errorNorthByIterr(1,:), errorVetiByIterr(1,:), pointSize, 'filled', 'r');
    hold on;
    scatter(errorNorthByIterr(2,:), errorVetiByIterr(2,:), pointSize, 'filled', 'b');
    legend('无钟差辅助','有钟差辅助');
    
    xlabel('北向误差/米', FontSize=fontSize);
    ylabel('天向误差/米', FontSize=fontSize);
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    xlim([-200 200]);
    ylim([-1000 200]);
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