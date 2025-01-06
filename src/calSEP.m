function [Vg_s, L1, L2, L3] = calSEP(...
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
                r     = norm(userECEFPosition - satECEFPositionList(ns, :));           
                dddx  = (f0/c)*(satECEFVelocityList(ns, 1)/r - (userECEFPosition(1)-satECEFPositionList(ns,1))*(userECEFPosition-satECEFPositionList(ns,:))*satECEFVelocityList(ns,:)'/r^3);                  
                dddy  = (f0/c)*(satECEFVelocityList(ns, 2)/r - (userECEFPosition(2)-satECEFPositionList(ns,2))*(userECEFPosition-satECEFPositionList(ns,:))*satECEFVelocityList(ns,:)'/r^3);
                dddz  = (f0/c)*(satECEFVelocityList(ns, 3)/r - (userECEFPosition(3)-satECEFPositionList(ns,3))*(userECEFPosition-satECEFPositionList(ns,:))*satECEFVelocityList(ns,:)'/r^3);
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
            if(withClcCalSEP == 1) 
                G = [G; drddx, drddy, drddz, 0, 0];
            elseif(withClcCalSEP == 0)
                G = [G; drddx, drddy, drddz];
            end
            
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
%     if(withClcShift == 1 && withClcDraft == 1)
%         [V,Q,U]    = eig(inv(G'*W*G));
%         [Q,order]  = sort(diag(Q), 'descend');
%         Q = diag(Q');
%         V = V(:,order);
%     else
%         [V,Q,U]    = eig(inv(G'*G));
%         [Q,order] = sort(diag(Q), 'descend');
%         Q = diag(Q');
%         V = V(:,order);
%     end
    [V,Q,U]    = eig(inv(G'*W*G));
    [Q,order]  = sort(diag(Q), 'descend');
    Q = diag(Q');
    V = V(:,order);
%     disp("V：")
%     matrixSout(V)
%     disp("Q:")
%     matrixSout(Q)
%     disp(['(v1,v2) = ', num2str(V(:,1)' * V(:,2))]);
%     disp(['(v2,v3) = ', num2str(V(:,2)' * V(:,3))]);
%     disp(['(v3,v4) = ', num2str(V(:,3)' * V(:,4))]);
%     disp(['(v4,v5) = ', num2str(V(:,4)' * V(:,5))]);
%     disp(['(v5,v1) = ', num2str(V(:,5)' * V(:,1))]);
%     disp(['norm(v1) = ', num2str(norm(V(:,1)'))]);
%     disp(['norm(v2) = ', num2str(norm(V(:,2)'))]);
%     disp(['norm(v3) = ', num2str(norm(V(:,3)'))]);
%     disp(['norm(v4) = ', num2str(norm(V(:,4)'))]);
%     disp(['norm(v5) = ', num2str(norm(V(:,5)'))]);
    
    Vg = V(1:3,1:3);
    Qg = Q(1:3,1:3);
    [Vg_s, R] = qr_gramschmidt(Vg);
    
%     disp("V1':")
%     matrixSout(V1')
%     disp("inv(V1):")
%     matrixSout(inv(V1))
%     disp("Q1:")
%     matrixSout(Q1)
%     disp(['tr(Q1) = ', num2str(trace(Q1))])
%     disp("V:")
%     matrixSout(V)
%     disp("Q:")
%     matrixSout(Q)
%     disp("==================");
%     disp('Vg:');
%     matrixSout(Vg);
%     disp('Qg:');
%     matrixSout(Qg);
%     disp('Vg_s:');
%     matrixSout(Vg_s);
%     disp('R:');
%     matrixSout(R);

    
    L1 = sqrt(Qg(1,1))*R(1,1);
    L2 = sqrt(Qg(2,2))*R(2,2)*sqrt(1-R(1,2)^2);
    L3 = sqrt(Qg(3,3))*R(3,3)*sqrt(1-R(1,3)^2)*sqrt(1-R(2,3)^2);

%     invgg = inv(G'*G);
%     pdop  = sqrt(invgg(1,1)+invgg(2,2)+invgg(3,3));
%     disp(['pdop = ',    num2str(pdop)]);
%     disp(['dop_22 = ', num2str(sqrt(DOP1^2+DOP2^2+DOP3^2))]);
    
%     disp(['Q(4,4) = ', num2str(sqrt(Q(4,4)))]);
%     disp(['TDOP = ', num2str(sqrt(invgg(4,4)))]);
end
