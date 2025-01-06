% 对比有误高度信息辅助 v
% 函数功能：对比有无高度信息辅助下低轨单星定位结果
tic
clc;clear;
format long g

addpath(genpath("../../coordinateTransformation"));
paraAll = simulation_Hcompare_para();
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

%%%% 观测量--时间图 %%%%
% prSatList = [];
% pdSatList = [];
% for ns = 1:Ns
%     errorPr   = paraAll.pseudoRangeError*randn;
%     prSat     = rangeList(ns) + errorPr + paraAll.satError*randn + paraAll.ionoError*randn;
%     prSatList = [prSatList, prSat];
%     errorPd   = paraAll.dopplerError*randn;
%     pdSat     = dopplerList(ns) + errorPd;
%     pdSatList = [pdSatList, pdSat];
% end
% tList    = 1:Ns;
% fontSize = 20;
% figure
%     subplot(1,2,1);
%     plot(tList , prSatList, 'linewidth', 3);
%     xlabel('观测时间/秒', FontSize=fontSize);
%     ylabel('伪距/米', FontSize=fontSize);
%     ax1 = gca;
%     ax1.XAxisLocation = 'origin';
%     ax1.YAxisLocation = 'origin';
%     title('用户3伪距观测量随时间变化图', FontSize=fontSize);
%     subplot(1,2,2);
%     plot(tList , pdSatList, 'linewidth', 3);
%     xlabel('观测时间/秒', FontSize=fontSize);
%     ylabel('多普勒频率/赫兹', FontSize=fontSize);
%     ax1 = gca;
%     ax1.XAxisLocation = 'origin';
%     ax1.YAxisLocation = 'origin';
%     title('用户3多普勒频率观测量随时间变化图', FontSize=fontSize);

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
% 无高度
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

% 有高度
errorHByIterr       = [];
errorHHoriByIterr   = [];
errorHEastByIterr   = [];
errorHEastByIterrJ  = [];
errorHNorthByIterr  = [];
errorHNorthByIterrJ = [];
errorHVetiByIterr   = [];
errorHVetiByIterrJ  = [];
errorHXByIterr      = [];
errorHYByIterr      = [];
errorHZByIterr      = [];
PDOPHByIterr        = [];
HDOPHByIterr        = [];
VDOPHByIterr        = [];
B0HByIterr          = [];
errorHXIteraMatrix  = zeros(paraAll.iterrT, 10);
errorHYIteraMatrix  = zeros(paraAll.iterrT, 10);
errorHZIteraMatrix  = zeros(paraAll.iterrT, 10);
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
    
    % input: realDisList realAngleList(deg)
%     prDMEList = [];
%     paVORList = [];
%     for nvd = 1:2
%         prDME     = realDisList(nvd) + paraAll.DMEError*randn;
%         prDMEList = [prDMEList, prDME];
%         
%         paVOR     = realAngleList(nvd) * (1 + paraAll.VORError*randn) * pi/180;
%         paVORList = [paVORList, paVOR];
%     end
    prDMEList = []; paVORList = []; obsVORDMEECEFPostionList = [];
    % output: prDMEList paVORList(rad)
    
    % input:userLLHPosition(3)
    height = userLLHPosition(3) + (paraAll.hError + paraAll.hErrorSys)*randn;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%无高度LSM%%%%%%%%%%%%%%%%%%%%%
    [x0List, y0List, z0List, B00List, B10List,...
    PDOP, HDOP, VDOP, iterTime] = ...
    LSM60(paraAll.f0,...
    prSatList, pdSatList, satECEFPositionList, satECEFVelocityList,... % 卫星观测量
    prDMEList, paVORList, obsVORDMEECEFPostionList,...                 % VOR/DME观测量
    height,...                                                         % 高度观测量
    paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
    paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
    paraAll.withVORFlag, paraAll.withDMEFlag, 0, paraAll.weightFlag,...
    userECEFPositionStart, 0, 0, paraAll.clcOrder); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    %%%%%%%%%%%%%%%%有高度LSM%%%%%%%%%%%%%%%%%%%%%
    [x0List, y0List, z0List, B00List, B10List,...
    PDOP, HDOP, VDOP, iterTime] = ...
    LSM60(paraAll.f0,...
    prSatList, pdSatList, satECEFPositionList, satECEFVelocityList,... % 卫星观测量
    prDMEList, paVORList, obsVORDMEECEFPostionList,...                 % VOR/DME观测量
    height,...                                                         % 高度观测量
    paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
    paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
    paraAll.withVORFlag, paraAll.withDMEFlag, 1, paraAll.weightFlag,...
    userECEFPositionStart, 0, 0, paraAll.clcOrder); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%数据统计%%%%%%%%%%%%%%%%%%
    calPos            = [x0List(end),y0List(end),z0List(end)];
    error             = norm(calPos-userECEFPosition);
    [errorEast, errorNorth, errorHorizontal, errorVertical]=...
    errorECEF2ENU(calPos, userECEFPosition);
    errorHByIterr       = [errorHByIterr,      error];
    errorHHoriByIterr   = [errorHHoriByIterr,  errorHorizontal];
    errorHEastByIterr   = [errorHEastByIterr,  errorEast];
    errorHEastByIterrJ  = [errorHEastByIterrJ,  abs(errorEast)];
    errorHNorthByIterr  = [errorHNorthByIterr, errorNorth];
    errorHNorthByIterrJ = [errorHNorthByIterrJ, abs(errorNorth)];
    errorHVetiByIterr   = [errorHVetiByIterr,  errorVertical];
    errorHVetiByIterrJ  = [errorHVetiByIterrJ,  abs(errorVertical)];
    errorHXByIterr      = [errorHXByIterr,     calPos(1)-userECEFPosition(1)];
    errorHYByIterr      = [errorHYByIterr,     calPos(2)-userECEFPosition(2)];
    errorHZByIterr      = [errorHZByIterr,     calPos(3)-userECEFPosition(3)];
    B0HByIterr          = [B0HByIterr,         B00List(end)];
    PDOPHByIterr        = [PDOPHByIterr,       PDOP];
    HDOPHByIterr        = [HDOPHByIterr,       HDOP];
    VDOPHByIterr        = [VDOPHByIterr,       VDOP];
    % x0List, y0List, z0List
    for itera = 1:length(x0List)
        errorHXIteraMatrix(iterr,itera) = x0List(itera)-userECEFPosition(1);
        errorHYIteraMatrix(iterr,itera) = y0List(itera)-userECEFPosition(2);
        errorHZIteraMatrix(iterr,itera) = z0List(itera)-userECEFPosition(3);
    end
    if(rem(iterr,50) == 0)
        disp(['第', num2str(iterr), '次仿真结束']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
errorByIterrSort        = sort(errorByIterr);
errorHoriByIterrSort    = sort(errorHoriByIterr);
errorVetiByIterrJSort   = sort(errorVetiByIterrJ);
errorEastByIterrJSort   = sort(errorEastByIterrJ);
errorNorthByIterrJSort  = sort(errorNorthByIterrJ);

errorHByIterrSort       = sort(errorHByIterr);
errorHHoriByIterrSort   = sort(errorHHoriByIterr);
errorHVetiByIterrJSort  = sort(errorHVetiByIterrJ);
errorHEastByIterrJSort  = sort(errorHEastByIterrJ);
errorHNorthByIterrJSort = sort(errorHNorthByIterrJ);

%%%%%%%% 95%误差+DOP %%%%%%%%
disp('-------------------------------')
disp('无高度辅助：')
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


disp('-------------------------------')
disp('有高度辅助：')
disp(['error(95%)       = ', num2str(errorHByIterrSort(paraAll.iterrT*0.95))]);
disp(['errorHori(95%)   = ', num2str(errorHHoriByIterrSort(paraAll.iterrT*0.95))]);    
disp(['errorVetiJ(95%)  = ', num2str(errorHVetiByIterrJSort(paraAll.iterrT*0.95))]); 
disp(['errorEastJ(95%)  = ', num2str(errorHEastByIterrJSort(paraAll.iterrT*0.95))]); 
disp(['errorNorthJ(95%) = ', num2str(errorHNorthByIterrJSort(paraAll.iterrT*0.95))]); 
disp(['error(65%)       = ', num2str(errorHByIterrSort(paraAll.iterrT*0.65))]);
disp(['errorHori(65%)   = ', num2str(errorHHoriByIterrSort(paraAll.iterrT*0.65))]);    
disp(['errorVetiJ(65%)  = ', num2str(errorHVetiByIterrJSort(paraAll.iterrT*0.65))]); 
disp(['errorEastJ(65%)  = ', num2str(errorHEastByIterrJSort(paraAll.iterrT*0.65))]); 
disp(['errorNorthJ(65%) = ', num2str(errorHNorthByIterrJSort(paraAll.iterrT*0.65))]); 

disp(['PDOP             = ', num2str(mean(PDOPHByIterr))]);
disp(['HDOP             = ', num2str(mean(HDOPHByIterr))]);
disp(['VDOP             = ', num2str(mean(VDOPHByIterr))]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%   无高度 SEP     %%%%%%%%
disp('---------无高度---------')
[V, L1, L2, L3] = calSEP(...
paraAll.f0,...
userECEFPosition,...
satECEFPositionList, satECEFVelocityList,[],...                                                           
paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
paraAll.withVORFlag, paraAll.withDMEFlag, 0,...
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

disp('---------有高度---------')
[V, L1, L2, L3] = calSEP(...
paraAll.f0,...
userECEFPosition,...
satECEFPositionList, satECEFVelocityList,[],...                                                           
paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
paraAll.withVORFlag, paraAll.withDMEFlag, 1,...
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
disp('-----------无高度-----------')
CRLBMatrix = calCRLB(...
paraAll.f0,...
userECEFPosition,...
satECEFPositionList, satECEFVelocityList,[],...                                                                     
paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
paraAll.withVORFlag, paraAll.withDMEFlag, 0,...
paraAll.withClcShift, paraAll.withClcDraft);

CRLBMatrix_ENU = Rl*CRLBMatrix*Rl';
 
posCRLB        = sqrt(CRLBMatrix(1,1)+CRLBMatrix(2,2)+CRLBMatrix(3,3));

disp(['error_std        =      ', num2str(sqrt(var(errorXByIterr)+var(errorYByIterr)+var(errorZByIterr)))]);
disp(['posCRLB          =      ', num2str(posCRLB)]);

disp(['errorX_std       =      ', num2str(sqrt(var(errorXByIterr)))]);
disp(['errorX_mean      =      ', num2str(mean(errorXByIterr))]);
disp(['errorY_std       =      ', num2str(sqrt(var(errorYByIterr)))]);
disp(['errorY_mean      =      ', num2str(mean(errorYByIterr))]);
disp(['errorZ_std       =      ', num2str(sqrt(var(errorZByIterr)))]);
disp(['errorZ_mean      =      ', num2str(mean(errorZByIterr))]);
disp(['errorX_CRLB      =      ', num2str(sqrt(CRLBMatrix(1,1)))]);
disp(['errorY_CRLB      =      ', num2str(sqrt(CRLBMatrix(2,2)))]);
disp(['errorZ_CRLB      =      ', num2str(sqrt(CRLBMatrix(3,3)))]);

disp(['errorEast_std    =      ', num2str(sqrt(var(errorEastByIterr  )))]);
disp(['errorEast_mean   =      ', num2str(mean(errorEastByIterr    ))]);
disp(['errorNorth_std   =      ', num2str(sqrt(var(errorNorthByIterr )))]);
disp(['errorNorth_mean  =      ', num2str(mean(errorNorthByIterr    ))]);
disp(['errorVeti_std    =      ', num2str(sqrt(var(errorVetiByIterr  )))]);
disp(['errorVeti_mean   =      ', num2str(mean(errorVetiByIterr  ))]);
disp(['errorEast_CRLB   =      ', num2str(sqrt(CRLBMatrix_ENU(1,1)))]);
disp(['errorNorthY_CRLB =      ', num2str(sqrt(CRLBMatrix_ENU(2,2)))]);
disp(['errorVeti_CRLB   =      ', num2str(sqrt(CRLBMatrix_ENU(3,3)))]);

disp('-----------有高度-----------')
CRLBMatrix = calCRLB(...
paraAll.f0,...
userECEFPosition,...
satECEFPositionList, satECEFVelocityList,[],...                                                                     
paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
paraAll.withVORFlag, paraAll.withDMEFlag, 1,...
paraAll.withClcShift, paraAll.withClcDraft);

CRLBMatrix_ENU = Rl*CRLBMatrix*Rl';
 
posCRLB        = sqrt(CRLBMatrix(1,1)+CRLBMatrix(2,2)+CRLBMatrix(3,3));

disp(['error_std        =      ', num2str(sqrt(var(errorHXByIterr)+var(errorHYByIterr)+var(errorHZByIterr)))]);
disp(['posCRLB          =      ', num2str(posCRLB)]);

disp(['errorX_std       =      ', num2str(sqrt(var(errorHXByIterr)))]);
disp(['errorX_mean      =      ', num2str(mean(errorHXByIterr))]);
disp(['errorY_std       =      ', num2str(sqrt(var(errorHYByIterr)))]);
disp(['errorY_mean      =      ', num2str(mean(errorHYByIterr))]);
disp(['errorZ_std       =      ', num2str(sqrt(var(errorHZByIterr)))]);
disp(['errorZ_mean      =      ', num2str(mean(errorHZByIterr))]);
disp(['errorX_CRLB      =      ', num2str(sqrt(CRLBMatrix(1,1)))]);
disp(['errorY_CRLB      =      ', num2str(sqrt(CRLBMatrix(2,2)))]);
disp(['errorZ_CRLB      =      ', num2str(sqrt(CRLBMatrix(3,3)))]);

disp(['errorEast_std    =      ', num2str(sqrt(var(errorHEastByIterr  )))]);
disp(['errorEast_mean   =      ', num2str(mean(errorHEastByIterr    ))]);
disp(['errorNorth_std   =      ', num2str(sqrt(var(errorHNorthByIterr )))]);
disp(['errorNorth_mean  =      ', num2str(mean(errorHNorthByIterr    ))]);
disp(['errorVeti_std    =      ', num2str(sqrt(var(errorHVetiByIterr  )))]);
disp(['errorVeti_mean   =      ', num2str(mean(errorHVetiByIterr  ))]);
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
	subplot(1,2,1)  
    scatter(errorEastByIterr, errorNorthByIterr, pointSize, 'filled', 'r');
    hold on;
    scatter(errorHEastByIterr, errorHNorthByIterr, pointSize, 'filled', 'g');
    xlabel('东向误差/米', FontSize=fontSize);
    ylabel('北向误差/米', FontSize=fontSize);
    legend('无高度辅助', '有高度辅助');
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    set(gca,'XTick',-800:200:800);
    set(gca,'XTicklabel',{'-800','-600','-400','-200','0','200','400','600','800'});
    set(gca,'XTick',-800:200:800);
    set(gca,'XTicklabel',{'-800','-600','-400','-200','0','200','400','600','800'});
    set(gca,'FontSize', 20);
    xlim([-800 800]);
    ylim([-800 800]);
    grid on;
    title('水平定位误差', FontSize=fontSize);

    subplot(1,2,2)
    scatter(errorNorthByIterr, errorVetiByIterr, pointSize, 'filled', 'r');
    hold on; 
    scatter(errorHNorthByIterr, errorHVetiByIterr, pointSize, 'filled', 'g');
    xlabel('北向误差/米', FontSize=fontSize);
    ylabel('天向误差/米', FontSize=fontSize);
    legend('无高度辅助', '有高度辅助');
    ax1 = gca;
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    set(gca,'XTick',-800:200:800);
    set(gca,'XTicklabel',{'-800','-600','-400','-200','0','200','400','600','800'});
    set(gca,'YTick',-800:200:800);
    set(gca,'YTicklabel',{'-800','-600','-400','-200','0','200','400','600','800'});
    set(gca,'FontSize', 20);
    xlim([-800 800]);
    ylim([-800 800]);
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