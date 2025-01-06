close all, clear all, clc, dbstop if error

h1=worldmap('china')
setm(h1,'mapprojection','lambert');%Բ��ͶӰ
setm(h1,'FFaceColor','w')%ͼ��

ChinaP=shaperead('bou1_4p.shp','UseGeoCoords',true)
ChinaL=shaperead('bou2_4l.shp','UseGeoCoords',true)
CapLon=[117.000923,115.48333,125.35000,127.63333,123.38333,111.670801,87.68333,103.73333,106.26667,112.53333,108.95000,113.65000,117.283042,119.78333,120.20000,118.30000,113.23333,115.90000,110.35000,108.320004,106.71667,113.00000,114.298572,104.06667,102.73333,91.00000,96.75000,117.20000,121.55333,106.45000,116.41667,121.30,114.10000,113.50000];
CapLat=[36.675807,38.03333,43.88333,47.75000,41.80000, 41.818311,43.76667,36.03333,37.46667,37.86667,34.26667,34.76667, 31.86119,32.05,30.266,26.0833,23.16667,28.68333,20.01667, 22.8240,26.56667,28.2166, 30.58435,30.6666,25.05000,30.600,36.5666,39.133,31.2000, 29.566,39.9166, 25.03,22.20,22.20];
names={'����','ʯ��ׯ','����','������','����', '���ͺ���','��³ľ��','����','����','̫ԭ','����','֣��','�Ϸ�','�Ͼ�','����','����','����','�ϲ�','����','����','����','��ɳ','�人','�ɶ�','����','����','����','���','�Ϻ�','����', '����','̨��','���','����'};

geoshow(ChinaP,'Facecolor',[1 1 0.5])%��ʾ��
geoshow(ChinaL,'LineStyle','-.','Color','k','LineWidth',1)%��ʾ����
geoshow(CapLat,CapLon,'DisplayType','point','Marker','.','MarkerEdgeColor','red')%ʡ���
geoshow(39.9166,116.41667,'DisplayType','point','Marker','p','MarkerEdgeColor','red')%�׶�
%ʡ���ע
for i=1:numel(names)
    textm(CapLat(i)+0.3,CapLon(i)+0.3,names(i),'color','k','FontSize',8)
end
%ͼ��
title('�й�������ͼ','FontSize',20);

%ͼ��
legend({'ʡ����','ʡ��','�׶�','������'},'FontSize',12,'Location','southwest')

%������
scaleruler('units','km')
setm(handlem('scaleruler1'),'RulerStyle','lines','MajorTick',0:500:1000,'MinorTick',0:250:500,'TickDir','down')

%���ϵ����Ϻ��
h2=axes('pos',[0.5922 0.15 0.13 0.2])
worldmap([5.559248066 20.549868679],[106.680363685 122.034461754])
setm(h2,'FFaceColor','w')
insert1=shaperead('bou2_4l.shp','UseGeoCoords',true)
geoshow([insert1.Lat],[insert1.Lon],'Color','k','LineWidth',1)
mlabel
plabel
gridm
setm(h2,'FFacecolor','c')
title('���ϵ����Ϻ��','FontSize',6)

% ָ����
northarrow('latitude',50,'longitude',62)
h=handlem('NorthArrow');
set(h,'FaceColor','k')





