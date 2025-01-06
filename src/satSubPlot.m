% 函数功能：计算星下点轨迹周期
clc;clear all;
format long g;
addpath(genpath("../../coordinateTransformation"));
tic
%%%%%%%%%%%% 星座参数 %%%%%%%%%%%%%%%%%%
%%%% 角度参数
beamAngle      = 55; 
shieldingAngle = 10; 
inclination    = 86.5 * pi / 180; 
%%%% 轨道高度
RE     = 6378137;
hs     = 1175000;
Rs     = RE + hs;
omegaS = sqrt(398601.2)/(Rs/1000)^(3/2);  % Satellite angular velocity[rad/s]
omegaE = 7.292115*10^-5;  % 地球自转角速度[rad/s]
%%%% 星座配置
for l = 1:6
    OMEGA0(l) = (l-1)*pi/6;
end
for l = 1:6
    for m = 1:10
        omega0(l, m) = (m-1)*pi/5+(l-1)*pi/10;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tsat    = 2*pi/omegaS;
Tearth  = 2*pi/omegaE;
disp(['卫星绕地球一周的时间', num2str(Tsat), '(秒):', second_change(Tsat)]);
ttEnd   = 24*24*60*60;
ttDelta = 0.05;
N       = ttEnd/ttDelta;
% latList = zeros(1, N);
% lonList = zeros(1, N);
l = 1;
m = 1;
for n = 1:N
    tt    = n*ttDelta;
    omega = omega0(l, m) + omegaS * tt;
    OMEGA = OMEGA0(l) - omegaE * tt;
    % tt时刻卫星的位置和速度
    satECEFPositionX = Rs*(cos(OMEGA)*cos(omega)-sin(OMEGA)*sin(omega)*cos(inclination));
    satECEFPositionY = Rs*(sin(OMEGA)*cos(omega)+cos(OMEGA)*sin(omega)*cos(inclination));
    satECEFPositionZ = Rs*sin(omega)*sin(inclination);
    satECEFPosition  = [satECEFPositionX, satECEFPositionY, satECEFPositionZ];
    satLLHPosition   = ecef2llh(satECEFPosition);
%     latList(n)       = satLLHPosition(2);
%     lonList(n)       = satLLHPosition(1);
    if((tt > 6532) && (abs(satLLHPosition(2)) <= 0.4) && (abs(satLLHPosition(1)) <= 0.4))
        disp(['星下点轨迹周期:', num2str(tt), '秒']);
    end
end
second_change(904796)
% latList(1)
% lonList(1)
% subStar(latList,lonList,'b');



% ERROR!!!
% Tsat    = 2*pi/omegaS;
% Tearth  = 2*pi/omegaE;
% for N = 1:1000
%     for M = 1:1000
%         if(abs(N * Tsat - M * Tearth ) < 10)
%             disp('N = ',num2str(N));
%             disp('M = ',num2str(M));
%             disp(num2str(N * Tsat));
%         end
%     end
% end

toc





