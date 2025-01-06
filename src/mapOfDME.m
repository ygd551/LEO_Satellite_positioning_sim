% close all;
% clear all;
% clc, dbstop if error;
function mapOfDME(userLLHPosition, VORDMELLHPostionList, names)
    h1 = worldmap('china');
%     h1 = worldmap([30 50], [110 130]);
    
    setm(h1, 'mapprojection', 'mercator');%圆柱投影
    % setm(h1, 'mapprojection', 'lambert');%正形圆锥投影
    setm(h1, 'FFaceColor', 'w'); %图廓


    ChinaP = shaperead('./ChinaMap/bou1_4p.shp', 'UseGeoCoords', true);  % 中国面文件
    ChinaL = shaperead('./ChinaMap/bou2_4l.shp', 'UseGeoCoords', true);  % 中国省界线文件（包含九段线）
    N = size(VORDMELLHPostionList, 1);
    % disp(['N = ', num2str(N)]);
    CapLon = VORDMELLHPostionList(:, 1)';
    CapLat = VORDMELLHPostionList(:, 2)';

    % CapLon = [117.000923,115.48333,125.35000,127.63333,123.38333,111.670801,87.68333,103.73333,106.26667,112.53333,108.95000,113.65000,117.283042,119.78333,120.20000,118.30000,113.23333,115.90000,110.35000,108.320004,106.71667,113.00000,114.298572,104.06667,102.73333,91.00000,96.75000,117.20000,121.55333,106.45000,116.41667,121.30,114.10000,113.50000];
    % CapLat = [36.675807,38.03333,43.88333,47.75000,41.80000, 41.818311,43.76667,36.03333,37.46667,37.86667,34.26667,34.76667, 31.86119,32.05,30.266,26.0833,23.16667,28.68333,20.01667, 22.8240,26.56667,28.2166, 30.58435,30.6666,25.05000,30.600,36.5666,39.133,31.2000, 29.566,39.9166, 25.03,22.20,22.20];
    % names = {'济南','石家庄','长春','哈尔滨','沈阳', '呼和浩特','乌鲁木齐','兰州','银川','太原','西安','郑州','合肥','南京','杭州','福州','广州','南昌','海口','南宁','贵阳','长沙','武汉','成都','昆明','拉萨','西宁','天津','上海','重庆', '北京','台北','香港','澳门'};

    % 'red'、'green'、'blue'、'cyan'、
    % 'magenta'、'yellow'、'black'、'white' 和 'none'。有效

    % geoshow(ChinaP, 'Facecolor', [1 1 0.5])  % 显示面
    geoshow(ChinaP, 'Facecolor', 'white')  % 显示面
    geoshow(ChinaL, 'LineStyle', '-.', 'Color', 'k', 'LineWidth', 1)  % 显示界线
    % VOR/DME台站点
    geoshow(CapLat, CapLon, 'DisplayType', 'point', 'Marker', '.', 'MarkerEdgeColor', 'red', 'MarkerSize', 10)
    % 用户位置
    geoshow(userLLHPosition(2), userLLHPosition(1), 'DisplayType', 'point', 'Marker', 'p', 'MarkerEdgeColor', 'green', 'MarkerFaceColor', 'green', 'MarkerSize', 10)
    % VOR/DME标注
    for i = 1:numel(names)
        textm(CapLat(i)+1, CapLon(i)+1, names(i), 'color', 'k', 'FontSize', 10)
    end
    %图名
    title('VOR/DME台站位置', 'FontSize', 20);

    %图例
    legend({'省界线', 'VOR/DME站', '用户位置'}, 'FontSize', 12, 'Location', 'southwest');

    %比例尺
    scaleruler('units','km')
    setm(handlem('scaleruler1'), 'RulerStyle', 'lines', 'MajorTick', 0:500:1000,'MinorTick',0:250:500,'TickDir','down')

%     %海南岛及南海诸岛
%     h2 = axes('pos', [0.5922 0.15 0.13 0.2]);
%     worldmap([5.559248066 20.549868679], [106.680363685 122.034461754])
%     setm(h2,'FFaceColor', 'w')
%     insert1 = shaperead('./ChinaMap/bou2_4l.shp', 'UseGeoCoords', true);
%     geoshow([insert1.Lat], [insert1.Lon], 'Color', 'k', 'LineWidth', 1);
%     mlabel
%     plabel
%     gridm
%     setm(h2,'FFacecolor', 'c');
%     title('海南岛及南海诸岛', 'FontSize', 6);

    % 指北针
%     northarrow('latitude', 50, 'longitude', 62);
%     h=handlem('NorthArrow');
%     set(h, 'FaceColor', 'k');
end