function [x0List, y0List, z0List, B00List, B10List,...
          PDOP, HDOP, VDOP, iterTime] = LSM602(...
f0,...
prSatList, pdSatList, satECEFPositionList, satECEFVelocityList,... % 卫星观测量
prDMEList, paVORList, obsVORDMEECEFPostionList,...                 % VOR/DME观测量
height,...                                                         % 高度观测量
clcList,...
pseudoRangeError, dopplerError, DMEError, VORError, hError,...% 观测误差
withLEOFlag, withLEOPrFlag, withLEOPdFlag,...
withVORFlag, withDMEFlag, withHeightFlag, weightFlag,...
withClc,...
userECEFPositionStart, B00In, B10In, clcOrder)                              % 最小二乘初始输入值
    
    addpath(genpath("../../coordinateTransformation"));
    c             = 299792458;
    Re            = 6378137;
    f             = 1/298.257222101;
    Rp            = Re * (1 - f);
    
    x0List        = [];
    y0List        = [];
    z0List        = [];
    B00List       = []; 
    B10List       = []; 
    
    x0            = userECEFPositionStart(1);
    y0            = userECEFPositionStart(2);
    z0            = userECEFPositionStart(3); 
    B00           = B00In;
    B10           = B10In;
    userPosition0 = [x0, y0, z0];
    x0List        = [x0List, x0];
    y0List        = [y0List, y0];
    z0List        = [z0List, z0];

    Ns  = length(prSatList);
    Nvd = length(prDMEList);

    if(withClc)
        for ns = 1:Ns
            prSatList(ns) = prSatList(ns) - clcList(ns);
        end
    end
    
    
    iterTime = 0;
    while(1)
        iterTime = iterTime + 1;
        if(iterTime > 200)
            break;
        end
        Z   = [];
        H   = [];
        W_v = [];
        
        % LEO
        if(withLEOFlag == 1)
            if(withLEOPrFlag == 1)
                for ns = 1:Ns
                    % prSat
                    rs0    = norm(userPosition0-satECEFPositionList(ns, :));
                    deltaC = prSatList(ns) - rs0;
                    dCdx   = (x0-satECEFPositionList(ns, 1))/rs0;
                    dCdy   = (y0-satECEFPositionList(ns, 2))/rs0;
                    dCdz   = (z0-satECEFPositionList(ns, 3))/rs0;
                    dCdB0  = 0;
                    dCdB1  = 0;
                    Z      = [Z; deltaC];
                    if(clcOrder == 0)
                        H  = [H; dCdx, dCdy, dCdz, dCdB0];
                    elseif(clcOrder == 1)
                        H  = [H; dCdx, dCdy, dCdz, dCdB0, dCdB1];
                    else
                        disp('errorClcOrder!!!');
                    end
                    if(weightFlag)
                        W_v = [W_v, 1/pseudoRangeError^2];
                    end
                end
            end
            if(withLEOPdFlag == 1)
                for ns = 1:Ns
                    % pdSat
                    rs0    = norm(userPosition0-satECEFPositionList(ns, :));
                    ds0    = (satECEFVelocityList(ns,:)*(userPosition0 - satECEFPositionList(ns,:))'/(c*rs0))*f0;            
                    deltaD = pdSatList(ns) - (ds0+B10);
                    dDdx   = (f0/c)*(satECEFVelocityList(ns,1)/rs0-(x0-satECEFPositionList(ns,1))*(userPosition0-satECEFPositionList(ns,:))*satECEFVelocityList(ns,:)'/rs0^3);                  
                    dDdy   = (f0/c)*(satECEFVelocityList(ns,2)/rs0-(y0-satECEFPositionList(ns,2))*(userPosition0-satECEFPositionList(ns,:))*satECEFVelocityList(ns,:)'/rs0^3);
                    dDdz   = (f0/c)*(satECEFVelocityList(ns,3)/rs0-(z0-satECEFPositionList(ns,3))*(userPosition0-satECEFPositionList(ns,:))*satECEFVelocityList(ns,:)'/rs0^3);
                    dDdB0  = 0;
                    dDdB1  = 1;
                    Z      = [Z; deltaD];
                    if(clcOrder == 0)
                        H = [H; dDdx, dDdy, dDdz, dDdB0];
                    elseif(clcOrder == 1)
                        H = [H; dDdx, dDdy, dDdz, dDdB0, dDdB1];
                    else
                        disp('errorClcOrder!!!');
                    end

                    if(weightFlag)
                        W_v = [W_v, 1/dopplerError^2];
                    end 
                end 
            end
        end
        
        % VOR
        if(withVORFlag)
            for nvd = 1:Nvd
                % paVORList
                angle0     = atan2(x0-obsVORDMEECEFPostionList(nvd, 1),y0-obsVORDMEECEFPostionList(nvd, 2));
                deltaAngle = paVORList(nvd) - angle0;
                R2         = (x0-obsVORDMEECEFPostionList(nvd, 1))^2+(y0-obsVORDMEECEFPostionList(nvd, 2))^2;
                dAdx       = (y0-obsVORDMEECEFPostionList(nvd, 2))/R2;
                dAdy       = -(x0-obsVORDMEECEFPostionList(nvd, 1))/R2;
                Z          = [Z; deltaAngle];
                
                if(clcOrder == 0)
                    H = [H; dAdx, dAdy, 0, 0];
                elseif(clcOrder == 1)
                    H = [H; dAdx, dAdy, 0, 0, 0];
                else
                    disp('errorClcOrder!!!');
                end
                if(weightFlag)
                    W_v = [W_v, 1/VORError^2];
                end
            end
        end
        
        % DME
        if(withDMEFlag)
            for nvd = 1:Nvd
                % prDME
                rd0      = norm(userPosition0 - obsVORDMEECEFPostionList(nvd,:));
                deltaDME = prDMEList(nvd) - rd0;
                dDMEdx   = (x0-obsVORDMEECEFPostionList(nvd,1))/rd0;
                dDMEdy   = (y0-obsVORDMEECEFPostionList(nvd,2))/rd0;
                dDMEdz   = (z0-obsVORDMEECEFPostionList(nvd,3))/rd0;
                Z = [Z; deltaDME];
                if(clcOrder == 0)
                    H = [H; dDMEdx, dDMEdy, dDMEdz, 0];
                elseif(clcOrder == 1)
                    H = [H; dDMEdx, dDMEdy, dDMEdz, 0, 0];
                else
                    disp('errorClcOrder!!!');
                end
                if(weightFlag)
                    W_v = [W_v, 1/DMEError^2];
                end
            end
        end
        
        % Baro
        if(withHeightFlag)
            for i = 1
                llh0      = ecef2llh(userPosition0);
                lon0      = llh0(1)*pi/180;
                lat0      = llh0(2)*pi/180;
                h0        = llh0(3);
                q0        = sqrt(x0^2+y0^2)-z0*cot(lat0);
                xv0       = q0*cos(lon0);
                yv0       = q0*sin(lon0);
                zv0       = 0;
                rv0       = norm([xv0,yv0,zv0]-userPosition0);
                deltaBaro = height - h0;
                dHdx      = (x0-xv0)/rv0;
                dHdy      = (y0-yv0)/rv0;
                dHdz      = (z0-zv0)/rv0;
                Z         = [Z; deltaBaro];
                if(clcOrder == 0)
                    H = [H; dHdx, dHdy, dHdz, 0];
                elseif(clcOrder == 1)
                    H = [H; dHdx, dHdy, dHdz, 0, 0];
                else
                    disp('errorClcOrder!!!');
                end
                if(weightFlag)
                    W_v = [W_v, 1/hError^2];
                end
            end
        end
        % 更新状态量
        W = diag(W_v);
        if(weightFlag)
            dX = pinv(H'*W*H)*H'*W*Z;
        else
            dX = pinv(H'*H)*H'*Z;
        end 
        x0  = x0  + dX(1);
        y0  = y0  + dX(2);
        z0  = z0  + dX(3);
        B00 = B00 + dX(4);
        if(clcOrder == 1)
        	B10 = B10 + dX(5);
        end
        userPosition0 = [x0, y0, z0];
        
        x0List  = [x0List,  x0];
        y0List  = [y0List,  y0];
        z0List  = [z0List,  z0];
        B00List = [B00List, B00];
        if(clcOrder == 1)
            B10List = [B10List, B10];
        end
        
        % 收敛判断
        if(abs(dX(1)) < 0.001 && abs(dX(2)) < 0.001 && abs(dX(3)) < 0.001)
            HcH_piv          = pinv(H'*H);
            userLLHPosition0 = ecef2llh(userPosition0);
            lon0rad          = userLLHPosition0(1)*pi/180;
            lat0rad          = userLLHPosition0(2)*pi/180;
            if(clcOrder == 0) 
                Rl = [
                -sin(lon0rad), cos(lon0rad), 0, 0;
                -sin(lat0rad)*cos(lon0rad), -sin(lat0rad)*sin(lon0rad), cos(lat0rad), 0;
                cos(lat0rad)*cos(lon0rad), cos(lat0rad)*sin(lon0rad), sin(lat0rad), 0;
                0 ,0 ,1 ,0;
                0 ,0 ,0 ,1];
            elseif(clcOrder == 1)
                Rl = [
                -sin(lon0rad), cos(lon0rad), 0, 0, 0;
                -sin(lat0rad)*cos(lon0rad), -sin(lat0rad)*sin(lon0rad), cos(lat0rad), 0, 0;
                cos(lat0rad)*cos(lon0rad), cos(lat0rad)*sin(lon0rad), sin(lat0rad), 0, 0;
                0 ,0 ,0 ,1, 0;
                0 ,0 ,0 ,0, 1];
            else
                Rl = 0;
            end
            HcH_piv_enu = Rl*pinv(H'*H)*Rl';
            PDOP        = sqrt(HcH_piv_enu(1,1)+HcH_piv_enu(2,2)+HcH_piv_enu(3,3));
            HDOP        = sqrt(HcH_piv_enu(1,1)+HcH_piv_enu(2,2));
            VDOP        = sqrt(HcH_piv_enu(3,3));
            break;
        end     
    end 
end