% 函数功能：绘制user4、5、12星空图轨迹
addpath(genpath("../../coordinateTransformation"));
userLLHPositionList ...  
= [93.247, 30, 1000;...
 108.247, 30, 1000;...
160.733, 85, 1000];

obsDataFilePathList...
= ["./data/user4.txt",...
    "./data/user5.txt",...
    "./data/user12.txt"];

figure(1)
for i = 1:3
    subplot(1,3,i)
    [ttList, satECEFPositionList, satECEFVelocityList,...
      rangeList, dopplerList] = getLEOObsData(obsDataFilePathList(i));
    [azimuthList, elevationList] = calAziAndEle(userLLHPositionList(i,:), satECEFPositionList);
    starsky(azimuthList, elevationList, '.');
end






