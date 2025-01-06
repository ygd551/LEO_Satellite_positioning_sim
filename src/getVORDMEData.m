function VORDMEData = getVORDMEData(userECEFPosition)
    % VORDMEData   pos(1:3) realrange realangle(0-360)
    addpath("../../coordinateTransformation")
    userLLHPosition               = ecef2llh(userECEFPosition);

    filePath                      = "./VORDMEPos.txt";
    fidin                         = fopen(filePath); 
    VORDMELLHPostionList          = [];
    VORDMEECEFPostionList         = []; 
    while(~feof(fidin))
        strline                   = fgetl(fidin);
        if(strline(1) == "#")
            llhPosStr             = strsplit(strline(2:17),', ');
            llhPos                = [str2double(llhPosStr(1)),str2double(llhPosStr(2)),str2double(llhPosStr(3))];
            VORDMELLHPostionList  = [VORDMELLHPostionList; llhPos];
            VORDMEECEFPostionList = [VORDMEECEFPostionList; llh2ecef(llhPos)]; 
        end
    end
    
    N                             = size(VORDMELLHPostionList,1);
    
    obsVORDMEECEFPostionList      = [];
    realDisList                   = [];
    realAngleList                 = [];
    obsVORDMEENUPosList           = [];
    
    effectiveRange                = 4.11*sqrt(userLLHPosition(3))*1852;
%     disp(['用户在该高度的有效DME观测距离为：',num2str(effectiveRange)]);
    for i = 1:N
        VORDMEENUPos = ecef2enu(VORDMEECEFPostionList(i,:), userECEFPosition, userLLHPosition);
        realDis      = norm(userECEFPosition - VORDMEECEFPostionList(i,:));
        
        realAngle    = atan2(VORDMEENUPos(2), VORDMEENUPos(1))*180/pi;
        if(VORDMEENUPos(1) == 0 && VORDMEENUPos(2) == 0)
            disp("VOR/DME位置错误");
            return;
        elseif(VORDMEENUPos(1) >= 0 || VORDMEENUPos(2) <= 0)
            realAngle = 90 - realAngle;
        else
            realAngle = 450 - realAngle;
        end
        
        if (realDis < effectiveRange)  
            obsVORDMEECEFPostionList    = [obsVORDMEECEFPostionList;...
                                           VORDMEECEFPostionList(i,:)];
            realDisList                 = [realDisList, realDis];
            realAngleList               = [realAngleList, realAngle];
            obsVORDMEENUPosList         = [obsVORDMEENUPosList; VORDMEENUPos];
        end
    end
    VORDMEData.obsVORDMEECEFPostionList = obsVORDMEECEFPostionList;
    VORDMEData.realDisList              = realDisList;
    VORDMEData.realAngleList            = realAngleList;
    VORDMEData.obsVORDMEENUPosList      = obsVORDMEENUPosList;
end