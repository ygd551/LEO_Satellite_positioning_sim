% 函数功能：绘制BDS-IGSO卫星星下点轨迹图
clc; clear;
t = 0:1:6;
we = 360/24;
u = we*t;
i = 30;
fai = asind( sind(i)*sind(u) );
deltalmd = atand( cosd(i)*tand(u) );

if(i==90)

deltalmd(end) = 90;

end

lmd = deltalmd - we*t;

% use symetry to generate the other data

for j = 1:6

lmd(j+7) = -lmd(7-j);

fai(j+7) = fai(7-j);

end

for j = 1:12

lmd(j+13) = lmd(13-j);

fai(j+13) = -fai(13-j);

end

h = geoshow('landareas.shp', 'FaceColor', [1 1 1]);

grid on

hold on

plot(lmd, fai); title(['GEO', num2str(i)])