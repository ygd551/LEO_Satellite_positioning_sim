function CRLBMatrix = calCRLB(...
f0,...
userECEFPosition,...
satECEFPositionList, satECEFVelocityList,... 
obsVORDMEECEFPostionList,...                                                                   
pseudoRangeError, dopplerError, DMEError, VORError, hError,...% 观测误差
withLEOFlag, withLEOPrFlag, withLEOPdFlag,...
withVORFlag, withDMEFlag, withHeightFlag,...
withClcShift, withClcDraft) 
    
    c = 299792458;
    
    Ns  = size(satECEFPositionList, 1);
    Nvd = size(obsVORDMEECEFPostionList, 1);
    
    G   = [];
    W_v = [];
    if(withLEOFlag == 1)
        if(withLEOPrFlag == 1)
            for ns = 1:Ns
                r     = norm(userECEFPosition - satECEFPositionList(ns, :));
                drdx  = (userECEFPosition(1) - satECEFPositionList(ns, 1))/r;
                drdy  = (userECEFPosition(2) - satECEFPositionList(ns, 2))/r;
                drdz  = (userECEFPosition(3) - satECEFPositionList(ns, 3))/r;
                if(withClcShift == 1 && withClcDraft == 1)
                    G = [G; drdx, drdy, drdz, 1, 0];
                elseif(withClcShift == 1 && withClcDraft == 0)
                    G = [G; drdx, drdy, drdz, 1];
                elseif(withClcShift == 0 && withClcDraft == 0)
                    G = [G; drdx, drdy, drdz];
                else 
                end
                W_v   = [W_v, 1/(pseudoRangeError)^2];
            end
        end
        if(withLEOPdFlag == 1)
            for ns = 1:Ns
                r    = norm(userECEFPosition - satECEFPositionList(ns, :));           
                dddx = (f0/c)*(satECEFVelocityList(ns, 1)/r - (userECEFPosition(1)-satECEFPositionList(ns,1))*(userECEFPosition-satECEFPositionList(ns,:))*satECEFVelocityList(ns,:)'/r^3);                  
                dddy = (f0/c)*(satECEFVelocityList(ns, 2)/r - (userECEFPosition(2)-satECEFPositionList(ns,2))*(userECEFPosition-satECEFPositionList(ns,:))*satECEFVelocityList(ns,:)'/r^3);
                dddz = (f0/c)*(satECEFVelocityList(ns, 3)/r - (userECEFPosition(3)-satECEFPositionList(ns,3))*(userECEFPosition-satECEFPositionList(ns,:))*satECEFVelocityList(ns,:)'/r^3);
                if(withClcShift == 1 && withClcDraft == 1)
                    G = [G; dddx, dddy, dddz, 0, 1];
                elseif(withClcShift == 1 && withClcDraft == 0)
                    G = [G; dddx, dddy, dddz, 0];
                elseif(withClcShift == 0 && withClcDraft == 1)
                    G = [G; dddx, dddy, dddz, 1];
                elseif(withClcShift == 0 && withClcDraft == 0)
                    G = [G; dddx, dddy, dddz];
                else 
                end
                W_v  = [W_v, 1/(dopplerError)^2];
            end
        end
    end
    if(withDMEFlag == 1)
        for nvd = 1:Nvd
            rd    = norm(userECEFPosition - obsVORDMEECEFPostionList(nvd,:));
            drddx = (userECEFPosition(1) - obsVORDMEECEFPostionList(nvd,1))/rd;
            drddy = (userECEFPosition(1) - obsVORDMEECEFPostionList(nvd,2))/rd;
            drddz = (userECEFPosition(1) - obsVORDMEECEFPostionList(nvd,3))/rd;
            G     = [G; drddx, drddy, drddz, 0, 0];
            W_v   = [W_v, 1/(DMEError)^2];
        end
    end
    if(withHeightFlag == 1)
        for i = 1
            userLLHPosition = ecef2llh(userECEFPosition);
            lon             = userLLHPosition(1)*pi/180;
            lat             = userLLHPosition(2)*pi/180;
            h               = userLLHPosition(3);
            q               = sqrt(userECEFPosition(1)^2+userECEFPosition(2)^2)-userECEFPosition(3)*cot(lat);
            xv              = q*cos(lon);
            yv              = q*sin(lon);
            zv              = 0;
            rv              = norm([xv,yv,zv]-userECEFPosition);
            dHdx            = (userECEFPosition(1)-xv)/rv;
            dHdy            = (userECEFPosition(2)-yv)/rv;
            dHdz            = (userECEFPosition(3)-zv)/rv;
            if(withClcShift == 1 && withClcDraft == 1)
                G = [G; dHdx, dHdy, dHdz, 0, 0];
            elseif(withClcShift == 1 && withClcDraft == 0)           
                G = [G; dHdx, dHdy, dHdz, 0];
            elseif(withClcShift == 0 && withClcDraft == 1)
                G = [G; dHdx, dHdy, dHdz, 0];
            elseif(withClcShift == 0 && withClcDraft == 0)
                G = [G; dHdx, dHdy, dHdz];
            else 
            end
            W_v = [W_v, 1/(hError)^2];
        end
    end
    W = diag(W_v);
    CRLBMatrix = inv(G'*W*G);
end
