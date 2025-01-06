% 函数功能：生成低轨卫星仿真观测数据文件
function satObsDataSim(userLLHPosition,...
    OMEGA0, omega0, lm,...
    obsTimeStart, obsTimeEnd,...
    dataFileName, para)

    addpath(genpath("../../../coordinateTransformation"));
    userECEFPosition = llh2ecef(userLLHPosition);
    
    %%% parameters from para              
    c              = para.c;
    B0             = para.beta0;
    B1             = para.beta1;
    B2             = para.beta2;
    beamAngle      = para.beamAngle;
    shieldingAngle = para.shieldingAngle; 
    inclination    = para.inclination*pi/180;
    RE             = 6378137;  
    hs             = para.hs;
    Rs             = RE + hs; 
    f0             = para.f0;
    hGPS           = 20200*10^3;
    omegaS         = sqrt(398601.2)/(Rs/1000)^(3/2); 
    omegaE         = 7.292115*10^-5; 
    %%% 新建观测文件，写文件头
    fid = fopen(dataFileName, 'at');
    strUserPosllh = "##接收机参考位置（LLH）"+num2str(userLLHPosition(1))+" "+num2str(userLLHPosition(2))+" "+num2str(userLLHPosition(3))+"\n";
    fprintf(fid, strUserPosllh);
    strUserPosecef = "##接收机参考位置（ECEF）"+num2str(userECEFPosition(1))+" "+num2str(userECEFPosition(2))+" "+num2str(userECEFPosition(3))+"\n";
    fprintf(fid, strUserPosecef);
    strInitialTime = "##initialTime     "+num2str(obsTimeStart)+"\n";
    strEndTime = "##endTime         "+num2str(obsTimeEnd)+"\n";
    fprintf(fid, strInitialTime);
    fprintf(fid, strEndTime);
    
    strObsTimeLength = "##obsTimeLength   "+num2str(obsTimeEnd-obsTimeStart+1)+"\n";
    fprintf(fid, strObsTimeLength);
    
    fprintf(fid, "##数据格式：satx saty satz satvx satvy satvz rangeObs dopplerObs\n");
    fprintf(fid, "##rangeObs = range + receiverClcError\n");
    fprintf(fid, "##dopplerObs = doppler + receiverClcDriftError\n");
    strSignalType = "##信号频率为:" + num2str(f0) + "赫兹\n";
    fprintf(fid, strSignalType);
    
    SN = size(lm,1);  % 卫星个数
    strSN = "##观测卫星个数："+num2str(SN)+"\n";
    fprintf(fid, strSN);
    
    fprintf(fid, '##end of head\n');
    num2strFormat = '%013.4f';
    
    tt = obsTimeStart;
    while(tt <= obsTimeEnd) 
        epochVisibleSatNum = 0;
        dataStrList = [];
        for sn = 1:SN
            l = lm(sn, 1);
            m = lm(sn, 2);
            omega = omega0(l,m) + omegaS*tt;
            OMEGA = OMEGA0(l) - omegaE*tt;
            prnNum = (l - 1)*10 + m;
            if(prnNum < 10)
                prnStr = "H0" + num2str(prnNum);
            else
                prnStr = "H" + num2str(prnNum);
            end
            % tt时刻卫星的位置和速度
            satECEFPositionX = Rs*(cos(OMEGA)*cos(omega) - sin(OMEGA)*sin(omega)*cos(inclination));
            satECEFPositionY = Rs*(sin(OMEGA)*cos(omega) + cos(OMEGA)*sin(omega)*cos(inclination));
            satECEFPositionZ = Rs*sin(omega)*sin(inclination);
            satECEFVelocityX = Rs*(-sin(OMEGA)*cos(omega)*(-omegaE) - cos(OMEGA)*sin(omega)*omegaS - cos(OMEGA)*sin(omega)*(-omegaE)*cos(inclination) - sin(OMEGA)*cos(omega)*omegaS*cos(inclination));                  
            satECEFVelocityY = Rs*(cos(OMEGA)*cos(omega)*(-omegaE) - sin(OMEGA)*sin(omega)*omegaS - sin(OMEGA)*sin(omega)*(-omegaE)*cos(inclination) - cos(OMEGA)*cos(omega)*omegaS*cos(inclination));
            satECEFVelocityZ = Rs*cos(omega)*sin(inclination)*omegaS;
            satECEFPosition  = [satECEFPositionX, satECEFPositionY, satECEFPositionZ];
            satECEFVelocity  = [satECEFVelocityX, satECEFVelocityY, satECEFVelocityZ];

            isVisible = satVisibleJudging(satECEFPosition, userECEFPosition, beamAngle, shieldingAngle);
            if isVisible
                epochVisibleSatNum = epochVisibleSatNum + 1;
                realRange = norm(satECEFPosition - userECEFPosition);
                pseudorange = realRange + B0 + B1*(tt - obsTimeStart) + B2*(tt - obsTimeStart)^2;
                dopplerFake = satECEFVelocity*(userECEFPosition - satECEFPosition)'/(c*realRange)*f0 + (B1 + 2*B2*(tt - obsTimeStart))*(f0/c);       
                dataStr = prnStr+" "+num2str(satECEFPositionX, num2strFormat) + " " + num2str(satECEFPositionY, num2strFormat) + " " + num2str(satECEFPositionZ, num2strFormat) + " " + ...
                                num2str(satECEFVelocityX, num2strFormat) + " " + num2str(satECEFVelocityY, num2strFormat) + " " + num2str(satECEFVelocityZ, num2strFormat) + " " + ...
                                num2str(pseudorange, num2strFormat) + " " + num2str(dopplerFake, num2strFormat) + "\n";
                dataStrList = [dataStrList, dataStr];
            end  
        end
        strEpochHead = "> " + num2str(tt) + " " + num2str(epochVisibleSatNum) + "\n";
        fprintf(fid, strEpochHead);
        for i = 1 : epochVisibleSatNum
            fprintf(fid, dataStrList(i));
        end

        tt = tt + para.tDelta  ;
    end
    fclose(fid);
end