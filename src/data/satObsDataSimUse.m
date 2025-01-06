% 函数功能：给定用户位置等参数
%          利用函数satObsDataSim函数，仿真得到低轨卫星观测数据
clc; clear all;
addpath(genpath("../../../coordinateTransformation"));
    
% user 1\2\3\4\5\6\7\8\9
userLLHPositionList = [93.503, 0, 1000;...
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
%%%% 参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 卫星星座配置
for l = 1:6
    OMEGA0(l) = (l-1)*pi/6;
end
for l = 1:6
    for m = 1:10
        omega0(l, m) = (m-1)*pi/5+(l-1)*pi/10;
    end
end
% 1 用户位置
% 1.1选择星下点轨迹附近用户1-12
% userIndex           = 4;
% userLLHPosition     = userLLHPositionList(userIndex, :);
% 1.2选择天津位置
% userLLHPosition     = [117.6, 38.5, 1000];
% userECEFPosition    = llh2ecef(userLLHPosition);
% 2 para
para.c              = 299792458;
para.beta0          = 7.5611;
% para.beta0          = 0;
para.beta1          = -0.005504;
% para.beta1          = 0;
% para.beta1          = 0.09;
para.beta2          = 2.7738E-08;
% para.beta2          = 0;

para.beamAngle      = 55;
para.shieldingAngle = 10;
para.inclination    = 86.5;
para.hs             = 1175000;

para.f0             = 1520*10^6;
para.tDelta         = 1;

% 4 选择观测的卫星
lm = [5, 1];
% 5 观测时间段
% sat1Intianjin_visbleTimeSectionSet =...  % sat(5,1)卫星观测时间段
% [[  5687,   6254]; [ 40078,  40793]; [ 46657,  47274];...  % 567, 715, 617
%  [ 84002,  84674]; [ 90545,  91218]; [125031, 125667];...  % 672, 673, 636
%  [131551, 132253]; [168992, 169561]; [175429, 176165];...  % 702, 569, 736
%  [210006, 210511]; [216459, 217209]; [254021, 254417]];    % 505, 750, 396
userTimeSection     = ...
    [[174816, 175351];[174700, 175460];[174809, 175344];...
     [175334, 175939];[175230, 176006];[175332, 175948];...
     [175809, 176537];[175778, 176560];[175825, 176569];...
     [176239, 177027];[176252, 177042];[176271, 177061]];
% obsTimeStart        = userTimeSection(userIndex, 1);
% obsTimeEnd          = userTimeSection(userIndex, 2);

obsTimeStart        = 5685;
obsTimeEnd          = 6236;

% 6 观测文件路径
dataFileName        = "./result.txt";

satObsDataSim(userLLHPosition,...
              OMEGA0, omega0, lm,...
              obsTimeStart, obsTimeEnd,...
              dataFileName, para);