function resultAnalysis(obsTimeList,errorByIterrM,...
    errorHoriByIterrM,...
    errorEastByIterrM,errorEastByIterrJM,...
    errorNorthByIterrM,errorNorthByIterrJM,...
    errorVetiByIterrM,errorVetiByIterrJM,...
    errorXByIterrM,errorYByIterrM,errorZByIterrM,...
    PDOPByIterrM,HDOPByIterrM,VDOPByIterrM,B0ByIterrM)

    figure;
    obsTimeLength = length(obsTimeList);
    for obsTimeIndex = 1:obsTimeLength
        obsTime               = obsTimeList(obsTimeIndex);
        disp('===============================================')
        disp(['观测时间为：', num2str(obsTime)]);
        errorByIterrSort      = sort(errorByIterrM(obsTimeIndex,:));
        errorHoriByIterrSort  = sort(errorHoriByIterrM(obsTimeIndex,:));
        errorVetiByIterrJSort = sort(errorVetiByIterrJM(obsTimeIndex,:));
        errorEastByIterrJSort = sort(errorEastByIterrJM(obsTimeIndex,:));
        errorNorthByIterrJSort = sort(errorNorthByIterrJM(obsTimeIndex,:));

        %%%%%%%% 95%误差+DOP %%%%%%%%
        disp(['error(95%)       = ', num2str(errorByIterrSort(paraAll.iterrT*0.95))]);
        disp(['errorHori(95%)   = ', num2str(errorHoriByIterrSort(paraAll.iterrT*0.95))]);    
        disp(['errorVetiJ(95%)  = ', num2str(errorVetiByIterrJSort(paraAll.iterrT*0.95))]); 
        disp(['errorEastJ(95%)  = ', num2str(errorEastByIterrJSort(paraAll.iterrT*0.95))]); 
        disp(['errorNorthJ(95%) = ', num2str(errorNorthByIterrJSort(paraAll.iterrT*0.95))]); 
        disp(['error(65%)       = ', num2str(errorByIterrSort(paraAll.iterrT*0.65))]);
        disp(['errorHori(65%)   = ', num2str(errorHoriByIterrSort(paraAll.iterrT*0.65))]);    
        disp(['errorVetiJ(65%)  = ', num2str(errorVetiByIterrJSort(paraAll.iterrT*0.65))]); 
        disp(['errorEastJ(65%)  = ', num2str(errorEastByIterrJSort(paraAll.iterrT*0.65))]); 
        disp(['errorNorthJ(65%) = ', num2str(errorNorthByIterrJSort(paraAll.iterrT*0.65))]); 

        disp(['PDOP             = ', num2str(mean(PDOPByIterrM(obsTimeIndex,:)))]);
        disp(['HDOP             = ', num2str(mean(HDOPByIterrM(obsTimeIndex,:)))]);
        disp(['VDOP             = ', num2str(mean(VDOPByIterrM(obsTimeIndex,:)))]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%    SEP     %%%%%%%%
        [V, L1, L2, L3] = calSEP(...
        paraAll.f0,...
        userECEFPosition,...
        satECEFPositionList, satECEFVelocityList,[],...                                                           
        paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
        paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
        paraAll.withVORFlag, paraAll.withDMEFlag, paraAll.withHeightFlag,...
        paraAll.withClcShift, paraAll.withClcDraft);
        disp('V:');
        matrixSout(V)
        disp(['L1           = ', num2str(L1)]);
        disp(['L2           = ', num2str(L2)]);
        disp(['L3           = ', num2str(L3)]);
        % 转到ENU
        lon0rad = userLLHPosition(1)*pi/180;
        lat0rad = userLLHPosition(2)*pi/180;
        R       = [-sin(lon0rad)         ,  cos(lon0rad)             , 0;
               -sin(lat0rad)*cos(lon0rad), -sin(lat0rad)*sin(lon0rad), cos(lat0rad);
                cos(lat0rad)*cos(lon0rad),  cos(lat0rad)*sin(lon0rad), sin(lat0rad)];
        if((paraAll.withClcShift + paraAll.withClcDraft) == 1)
           Rl       = [-sin(lon0rad)         ,  cos(lon0rad)             , 0           , 0;
                   -sin(lat0rad)*cos(lon0rad), -sin(lat0rad)*sin(lon0rad), cos(lat0rad), 0;
                    cos(lat0rad)*cos(lon0rad),  cos(lat0rad)*sin(lon0rad), sin(lat0rad), 0;
                                            0,                          0,            0, 1];
        else
           Rl       = [-sin(lon0rad)         ,  cos(lon0rad)             ,            0, 0, 0;
                   -sin(lat0rad)*cos(lon0rad), -sin(lat0rad)*sin(lon0rad), cos(lat0rad), 0, 0;
                    cos(lat0rad)*cos(lon0rad),  cos(lat0rad)*sin(lon0rad), sin(lat0rad), 0, 0;
                                            0,                          0,            0, 1, 0;
                                            0,                          0,            0, 0, 1];

        end
        V_ENU   = R * V;
        v1_ENU  = V_ENU(:,1);
        v2_ENU  = V_ENU(:,2);
        v3_ENU  = V_ENU(:,3);

        disp('三轴转到东北天方向上的结果')
        matrixSout(V_ENU);
        disp('东北天方向上的误差(1sigma):')
        disp(['sigmaEast        = ', num2str(sqrt((V_ENU(1,1)*L1)^2 + (V_ENU(1,2)*L2)^2 +(V_ENU(1,3)*L3)^2))])
        disp(['sigmaNorth       = ', num2str(sqrt((V_ENU(2,1)*L1)^2 + (V_ENU(2,2)*L2)^2 +(V_ENU(2,3)*L3)^2))])
        disp(['sigmaVeti        = ', num2str(sqrt((V_ENU(3,1)*L1)^2 + (V_ENU(3,2)*L2)^2 +(V_ENU(3,3)*L3)^2))])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%    CRLB     %%%%%%%%
        CRLBMatrix = calCRLB(...
        paraAll.f0,...
        userECEFPosition,...
        satECEFPositionList, satECEFVelocityList,[],...                                                                     
        paraAll.pseudoRangeError, paraAll.dopplerError, paraAll.DMEError, paraAll.VORError, paraAll.hError,...% 观测误差
        paraAll.withLEOFlag, paraAll.withLEOPrFlag, paraAll.withLEOPdFlag,...
        paraAll.withVORFlag, paraAll.withDMEFlag, paraAll.withHeightFlag,...
        paraAll.withClcShift, paraAll.withClcDraft);

        CRLBMatrix_ENU = Rl*CRLBMatrix*Rl';

        posCRLB        = sqrt(CRLBMatrix(1,1)+CRLBMatrix(2,2)+CRLBMatrix(3,3));

        disp(['error_std        =      ', num2str(sqrt(var(errorXByIterrM(obsTimeIndex,:))+var(errorYByIterrM(obsTimeIndex,:))+var(errorZByIterrM(obsTimeIndex,:))))]);
        disp(['posCRLB          =      ', num2str(posCRLB)]);

        disp(['errorX_std       =      ', num2str(sqrt(var(errorXByIterrM(obsTimeIndex,:))))]);
        disp(['errorX_mean      =      ', num2str(mean(errorXByIterrM(obsTimeIndex,:)))]);
        disp(['errorY_std       =      ', num2str(sqrt(var(errorYByIterrM(obsTimeIndex,:))))]);
        disp(['errorY_mean      =      ', num2str(mean(errorYByIterrM(obsTimeIndex,:)))]);
        disp(['errorZ_std       =      ', num2str(sqrt(var(errorZByIterrM(obsTimeIndex,:))))]);
        disp(['errorZ_mean      =      ', num2str(mean(errorZByIterrM(obsTimeIndex,:)))]);
        disp(['errorX_CRLB      =      ', num2str(sqrt(CRLBMatrix(1,1)))]);
        disp(['errorY_CRLB      =      ', num2str(sqrt(CRLBMatrix(2,2)))]);
        disp(['errorZ_CRLB      =      ', num2str(sqrt(CRLBMatrix(3,3)))]);

        disp(['errorEast_std    =      ', num2str(sqrt(var(errorEastByIterrM(obsTimeIndex,:)  )))]);
        disp(['errorEast_mean   =      ', num2str(mean(errorEastByIterrM(obsTimeIndex,:)    ))]);
        disp(['errorNorth_std   =      ', num2str(sqrt(var(errorNorthByIterrM(obsTimeIndex,:) )))]);
        disp(['errorNorth_mean  =      ', num2str(mean(errorNorthByIterrM(obsTimeIndex,:)    ))]);
        disp(['errorVeti_std    =      ', num2str(sqrt(var(errorVetiByIterrM(obsTimeIndex,:)  )))]);
        disp(['errorVeti_mean   =      ', num2str(mean(errorVetiByIterrM(obsTimeIndex,:)  ))]);
        disp(['errorEast_CRLB   =      ', num2str(sqrt(CRLBMatrix_ENU(1,1)))]);
        disp(['errorNorthY_CRLB =      ', num2str(sqrt(CRLBMatrix_ENU(2,2)))]);
        disp(['errorVeti_CRLB   =      ', num2str(sqrt(CRLBMatrix_ENU(3,3)))]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%% 误差分布（ENU + EN）%%%%%%%%
        scale       = 1;
        lineWidth   = 3;
        maxHeadSize = 1;
        fontSize    = 20;
        pointSize   = 6;

        range1      = 800;
        range2      = 800;

            subplot(2,6, (obsTimeIndex-1)*3 + 1)
            scatter3(errorEastByIterrM(obsTimeIndex,:), errorNorthByIterrM(obsTimeIndex,:), errorVetiByIterrM(obsTimeIndex,:), pointSize, 'filled', 'r');
        %     hold on;
        %     quiver3(0,0,0,scale*v1_ENU(1,1),scale*v1_ENU(2,1),scale*v1_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
        %     hold on;
        %     quiver3(0,0,0,scale*v2_ENU(1,1),scale*v2_ENU(2,1),scale*v2_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
        %     hold on; 
        %     quiver3(0,0,0,scale*v3_ENU(1,1),scale*v3_ENU(2,1),scale*v3_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
        %     hold on;
        %     quiver3(0,0,0,-scale*v1_ENU(1,1),-scale*v1_ENU(2,1),-scale*v1_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
        %     hold on;
        %     quiver3(0,0,0,-scale*v2_ENU(1,1),-scale*v2_ENU(2,1),-scale*v2_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
        %     hold on; 
        %     quiver3(0,0,0,-scale*v3_ENU(1,1),-scale*v3_ENU(2,1),-scale*v3_ENU(3,1), 0,Color = 'g', LineWidth = lineWidth, MaxHeadSize = maxHeadSize);
            xlim([-range1 range1])
            ylim([-range1 range1])
            zlim([-range1 range1])

            xlabel('东向误差/(m)', FontSize=fontSize);
            ylabel('北向误差/(m)', FontSize=fontSize);
            zlabel('天向误差/(m)', FontSize=fontSize);
            grid on;
            title('定位误差', FontSize=fontSize);

            subplot(2,6, (obsTimeIndex-1)*3 + 2)
            scatter(errorEastByIterrM(obsTimeIndex,:), errorNorthByIterrM(obsTimeIndex,:), pointSize, 'filled', 'r');
            xlabel('东向误差/(m)', FontSize=fontSize);
            ylabel('北向误差/(m)', FontSize=fontSize);
            ax1 = gca;
            ax1.XAxisLocation = 'origin';
            ax1.YAxisLocation = 'origin';
            xlim([-range1 range1]);
            ylim([-range1 range1]);
            grid on;
            title('水平定位误差', FontSize=fontSize);

            subplot(2,6, (obsTimeIndex-1)*3 + 3)
            scatter(errorNorthByIterrM(obsTimeIndex,:), errorVetiByIterrM(obsTimeIndex,:), pointSize, 'filled', 'r');
            xlabel('北向误差/(m)', FontSize=fontSize);
            ylabel('天向误差/(m)', FontSize=fontSize);
            ax1 = gca;
            ax1.XAxisLocation = 'origin';
            ax1.YAxisLocation = 'origin';
            xlim([-range2 range2]);
            ylim([-range2 range2]);
            grid on;
            title('北向天向误差', FontSize=fontSize);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
end