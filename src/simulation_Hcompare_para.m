function para = simulation_Hcompare_para()

    % 用户位置
    tjLLHPosition         = [117.6, 38.5, 1000];
    csLLHPosition         = [112.97, 28.18, 1000];
    user5LLHPosition      = [108.247, 30, 1000];
    user85LLHPosition     = [120.4, 84.9, 1000];
    userLLHPositionList   = [93.503, 0, 1000;...
                             108.503, 0, 1000;...
                             123.503, 0, 1000;...
                             93.247, 30, 1000;...
                             108.247, 30, 1000;...
                             123.247, 30, 1000;...
                             95.001, 60, 1000;...
                             110.001, 60, 1000;...
                             125.001, 60, 1000;...
                             130.733, 85, 1000;...
                             145.733, 85, 1000;...
                             160.733, 85, 1000];
    para.userIndex        = 4;  % 4、5、12
    para.userLLHPosition  = userLLHPositionList(para.userIndex, :);
    % para.userLLHPosition  = tjLLHPosition;
    para.obsSat           = "sat5_01";
    % 初始确定方式:0 最小残差法  1 自定义误差法  2 [0,0,0]
    para.positionStartWay = 0;
    para.errorxyz         = [20000,20000,10000];  % positionStartWay = 1 时有效;
    % 仿真次数
    para.iterrT           = 1000;
    % 计算SEP时是否考虑钟差,钟漂 0:不考虑  1:考虑
    para1                 = 1;
    para2                 = 0;
    
    para.withClcShift     = para1;   % = withLEOPrFlag
    para.withClcDraft     = para2;   % = withLEOPdFlag
    % 是否用LEO
    para.withLEOFlag      = 1;
    para.withLEOPrFlag    = para1;
    para.withLEOPdFlag    = para2;
    % 是否使用高度计
    para.withHeightFlag   = 0;
    % 是否钟差辅助
    para.withClcFlag   = 0;
    % 是否启用VOR/DME 
    para.withVORFlag      = 0;
    para.withDMEFlag      = 0;
    % 观测文件路径
    para.obsTime          = "5685_6236";
    para.obsTimeList      = 300:5:550;
%     para.obsTime          = "168977_169551";
    para.timeStart        = 5685;
    % 观测文件
%     para.obsDataFilePath  = "./data/user4_noclc.txt";
%     para.obsDataFilePath  = "./data/user4.txt";
    % para.obsDataFilePath  = "./data/tj_5685_6236.txt";
%     para.B1               = 0.09;
    para.obsDataFilePath  = "./data/user" + num2str(para.userIndex) + ".txt";
    % para.obsDataFilePath  = "./data/tj_signal1_168977_169551_unif.txt";
%     para.obsDataFilePath  = "obsData/tj_signal1_5685_6236_noclc_unif.txt";
%     para.obsDataFilePath  = "obsData/tj_signal1_5685_6236_onlyshift_unif.txt";
%     para.obsDataFilePath  = "/Users/yaoguodong/aDesk/myLEOExp/60sats/01代码/satObsData/substar5/user5_1s.txt";
%     para.obsDataFilePath  = "obsData/user85_1_12945_13736_unif.txt";
    % 接收机钟差模型阶数 0:BO  1:B0 B1
    if(para.withLEOPdFlag == 1)
        para.clcOrder     = 1;
    else
        para.clcOrder     = 0;
    end
        
    para.beta0            = 7.5611;
    para.clcDraftList     = [0.03, 0.09, 0.15, 0.21, 0.27, 0.33, 0.39, 0.45, 0.51, 0.57, 0.63];
    para.beta2            = 2.7738E-08;
%     para.clcError         = 0.05;
    para.clcError         = 0;
    % LEO 伪距多普勒观测误差
    para.pseudoRangeError = 4.5;
    para.dopplerError     = 1;
    % 载波频率
    para.f0               = 1520000000;
%     para.f0               = 1521000000;
    % 星历星钟误差
    para.satError         = 0;
    % 电离层误差
    para.ionoError        = 0;
    % DME误差设置
    para.DMEError         = 185;  % (m 1sigma)
    % VOR误差设置
    para.VORError         = 1*pi/180;  % (rad 1sigma)
    % barometer随机误差设置
    para.hError           = sqrt(30);  % (m 1sigma)
%     para.hErrorSys        = 0.1*userLLHPosition(3);
    para.hErrorSys        = 0;
    % 是否加权
    para.weightFlag       = 1;
    
    %% 参数输出提示
    disp("================参数设置================");
    if(para.userLLHPosition == tjLLHPosition)
        disp("用户位置在天津");
    elseif(para.userLLHPosition == csLLHPosition)
        disp("用户位置在长沙");
    else
        disp("未知用户位置！！！");
    end
    disp("观测卫星:" + para.obsSat);
    disp("观测时间段:" + para.obsTime);
%     disp("信号类型:" + para.signalType);
    disp("接收机钟差模型阶数:" + num2str(para.clcOrder) + "阶");
    if(para.weightFlag == 1)
        disp("加权方式");
    else
        disp("不加权");    
    end
    
    if(para.withLEOFlag == 1)
        if(para.withLEOPrFlag == 1 && para.withLEOPdFlag == 1)
            disp("使用LEO 伪距+多普勒");
        elseif(para.withLEOPrFlag == 1)
            disp("使用LEO 伪距");
        elseif(para.withLEOPdFlag == 1)
            disp("使用LEO 多普勒");
        else
        end    
    else
        disp("不使用LEO");
    end
    if(para.withHeightFlag == 1)
        disp("使用高度计");
    else
        disp("不使用高度计");
    end
    if(para.withVORFlag == 1)
        disp("使用VOR");
    else
        disp("不使用VOR");
    end
    if(para.withDMEFlag == 1)
        disp("使用DME");
    else
        disp("不使用DME");
    end
    if(para.positionStartWay == 0)
        disp("最小残差法确定初始值");
    elseif(para.positionStartWay == 1)
        disp("自定义误差法确定初始值");
    elseif(para.positionStartWay == 2)
        disp("初始值为[0,0,0]");
    else
    end
    disp("仿真重复次数:" + num2str(para.iterrT));
    disp("========================================");
end