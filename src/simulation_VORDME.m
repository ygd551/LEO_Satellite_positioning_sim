% VOR/DME定位结果   v
tic
clc;clear;
format long g

addpath(genpath("../../coordinateTransformation"));
paraAll = simulation_VORDME_para();
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
% starskyDraw(userLLHPosition, satECEFPositionList)

%%%%% 所有理论上可见的VOR/DME得到的观测量%%%%%%
VORDMEDataAll               = getVORDMEData(userECEFPosition);
obsVORDMEECEFPostionListAll = VORDMEDataAll.obsVORDMEECEFPostionList;
realDisListAll              = VORDMEDataAll.realDisList;
realAngleListAll            = VORDMEDataAll.realAngleList;
Nvd                         = length(realDisListAll);
disp(['用户可见陆基台站个数：',num2str(Nvd)]);

% % 给出所有台站相对于用户的位置图
enuList = [];
for i = 1 : size(obsVORDMEECEFPostionListAll, 1)
    enu = ecef2enu(obsVORDMEECEFPostionListAll(i,:), userECEFPosition, userLLHPosition);
    enuList = [enuList; enu];
end

figure
    sz = 100;
    c = 'r';
    scatter(enuList(:,1)', enuList(:,2)', sz, c, 'filled');
    set(gca,'XAxisLocation','origin');   
    set(gca,'YAxisLocation','origin'); % 将y轴的位置设置在x=0处。
    set(gca,'XDir','normal');  % 将x轴方向设置为普通(从左到右递增)              
    % set(gca,'XDir','reverse');  % 将x轴方向设置为反向(从右到左递增)。                  
    xlabel('东/米');
    ylabel('北/米');
    xlim([-3*10^5 3*10^5]);
    ylim([-3*10^5 3*10^5]);
    for i = 1 : length(enuList(:, 1)')
        if(i == 3 || i == 5)
            text(enuList(i, 1) + 300, enuList(i, 2) + 300, ['D',num2str(i)], 'Color', 'g');
        end
        text(enuList(i, 1) + 300, enuList(i, 2) + 300, ['D',num2str(i)]);
    end
    title('台站相对于用户的位置（水平）');

% 在地图上标出这些位置
obsVORDMELLHPostionListAll = [];
for i = 1 : length(obsVORDMEECEFPostionListAll)
    obsVORDMELLHPostionListAll = [obsVORDMELLHPostionListAll; ecef2llh(obsVORDMEECEFPostionListAll(i, :))];
end
% figure(2)
%     names = {'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7'};
%     mapOfDME(userLLHPosition, obsVORDMELLHPostionListAll, names);

% 选择台站(2,5) (3,5)  
index1 = 3;
index2 = 5;
obsVORDMEECEFPostionList = [obsVORDMEECEFPostionListAll(index1, :);
                            obsVORDMEECEFPostionListAll(index2, :)];
realDisList              = [realDisListAll(index1), realDisListAll(index2)];
realAngleList            = [realAngleListAll(index1), realAngleListAll(index2)];
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
    prDMEList = [];
    paVORList = [];
    for nvd = 1:2
        prDME     = realDisList(nvd) + paraAll.DMEError*randn;
        prDMEList = [prDMEList, prDME];

        paVOR     = realAngleList(nvd) * (1 + paraAll.VORError*randn) * pi/180;
        paVORList = [paVORList, paVOR];
    end
     % output: prDMEList paVORList(rad)
    
    % input:userLLHPosition(3)
    height = userLLHPosition(3) + (paraAll.hError + paraAll.hErrorSys)*randn;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%LSM%%%%%%%%%%%%%%%%%%%%%
    [x0List, y0List, z0List, B00List, B10List,...
    PDOP, HDOP, VDOP, iterTime] = ...
    LSM60(paraAll.f0,...
    prSatList, pdSatList, satECEFPositionList, satECEFVelocityList,... % 卫星观测量
    prDMEList, paVORList, obsVORDMEECEFPostionList,...                 % VOR/DME观测量
    height,...                                                         % 高度观测量
    paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
    paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
    0, paraAll.withDMEFlag, paraAll.withHeightFlag, paraAll.weightFlag,...
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
    if(rem(iterr,50) == 0)
        disp(['第', num2str(iterr), '次仿真结束']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
errorByIterrSort      = sort(errorByIterr);
errorHoriByIterrSort  = sort(errorHoriByIterr);
errorVetiByIterrJSort = sort(errorVetiByIterrJ);
errorEastByIterrJSort = sort(errorEastByIterrJ);
errorNorthByIterrJSort = sort(errorNorthByIterrJ);

%%%%%%%% 95%误差+DOP %%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%    SEP     %%%%%%%%
[V, L1, L2, L3] = calSEP(...
paraAll.f0,...
userECEFPosition,...
satECEFPositionList, satECEFVelocityList,[],...                                                           
paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
0, paraAll.withDMEFlag, paraAll.withHeightFlag,...
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
0, paraAll.withDMEFlag, paraAll.withHeightFlag,...
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% 误差分布（ENU + EN）%%%%%%%%
scale       = 1;
lineWidth   = 3;
maxHeadSize = 1;
fontSize    = 20;
pointSize   = 6;

range1      = 400;
range2      = 400;
figure;
	subplot(1,2,1)  
    scatter(errorEastByIterr, errorNorthByIterr, pointSize, 'filled', 'r');
    xlabel('东向误差/米', FontSize=fontSize);
    ylabel('北向误差/米', FontSize=fontSize);
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
    xlabel('北向误差/米', FontSize=fontSize);
    ylabel('天向误差/米', FontSize=fontSize);
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

         
toc