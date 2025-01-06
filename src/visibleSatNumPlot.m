clc;clear all;
format long g;
addpath(genpath("../../coordinateTransformation"));

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
stepLon           = 1;
stepLat           = 1;
userLonList       = -179:stepLon:179;
userLatList       = -89:stepLat:89;
tt                = 0;
mycolorList = [240/255 255/255 255/255
                252/255 230/255 201/255
                127/255 255/255 212/255
                1 1 0
                1 0 0];
%                 1 0.3804 0];
ttEnd             = 60;

for i = 1:length(userLonList)
    for j = 1:length(userLatList)
        userLLHPosition  = [userLonList(i)+180, userLatList(j), 0];
        userECEFPosition = llh2ecef(userLLHPosition);
        visibleNumTime   = zeros(1, ttEnd);
        for tt = 1:ttEnd
            visibleNum = 0;
            for l = 1:6
                for m = 1:10
                    omega = omega0(l, m) + omegaS * tt;
                    OMEGA = OMEGA0(l) - omegaE * tt;
                    % tt时刻卫星的位置和速度
                    satECEFPositionX = Rs*(cos(OMEGA)*cos(omega)-sin(OMEGA)*sin(omega)*cos(inclination));
                    satECEFPositionY = Rs*(sin(OMEGA)*cos(omega)+cos(OMEGA)*sin(omega)*cos(inclination));
                    satECEFPositionZ = Rs*sin(omega)*sin(inclination);
                    satECEFPosition  = [satECEFPositionX, satECEFPositionY, satECEFPositionZ];
                    isVisible = satVisibleJudging(satECEFPosition, userECEFPosition, beamAngle, shieldingAngle);
                    if isVisible == 1
                        visibleNum = visibleNum + 1;
                    end    
                end
            end
            visibleNumTime(tt) = visibleNum;
        end
        visibleNumMatrix(i, j) = mean(visibleNumTime);
        disp(['point:', num2str((i-1)*length(userLonList)+j)]);
    end
end

% size(visibleNumMatrix)
% max(max(visibleNumMatrix))
% min(min(visibleNumMatrix))

figure(1)
    
    for i = 1:length(userLonList)
        for j = 1:length(userLatList)
            hold on;
            if((visibleNumMatrix(i,j) > 0) && (visibleNumMatrix(i,j) <= 2))
                color = mycolorList(1,:);
            elseif((visibleNumMatrix(i,j) > 2) && (visibleNumMatrix(i,j) <= 4))
                color = mycolorList(2,:);
            elseif((visibleNumMatrix(i,j) > 4) && (visibleNumMatrix(i,j) <= 6))
                color = mycolorList(3,:);
            elseif((visibleNumMatrix(i,j) > 6) && (visibleNumMatrix(i,j) <= 8))
                color = mycolorList(4,:);
            else
                color = mycolorList(5,:);
            end
                
            rectangle('Position',[userLonList(i)-stepLon/2, userLatList(j)-stepLat/2, stepLon, stepLat], 'FaceColor', color, 'EdgeColor' , color)
  
        end
    end
    
    geoshow('landareas.shp', 'FaceColor', [1 1 1]);
    xlim([-180 180]);
    ylim([-90 90]);
    alpha(0.3);
    grid on;
    
    
    colormap(mycolorList)
    cb = colorbar;
    caxis([0 10]); 
    set(cb,'YTick',0:2:10); %色标值范围及显示间隔
    set(cb,'YTickLabel',{'0','2','4','6','8','10'}) %具体刻度赋值
    
saveas(gcf, 'a.png');

    
    
    
    