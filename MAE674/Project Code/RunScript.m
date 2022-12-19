%% Run Script for Single Run or Monte Carlo Run
clear
close all
clc

GenerateMeasurements
numruns             = 1;
filename            = 'AttitudeFilterCompare.xlsx';

if numruns == 1
    Bhat0           = [0 0 0]';
    qhat0           = (sqrt(2)/2)*[0, 1, 1, 0]'; 
    
    % Initial attitude covaraince - 1 deg^2 
    P0a             = 1 * (pi/180)^2 *ones(3,1); 

    % Initial bias covariance - .05 (deg/hr)^2
    P0b             = .05 *(pi/(180*3600))^2 *ones(3,1);
    
    tic
    MEKF
    EKF.RunTime = toc;
    
    pause(1);
    
    tic
    USQUE
    UKF.RunTime = toc;
    
    %% Plot Results EKF
    figure
    subplot(4,1,1)
    plot(t,EKF.xhat(1,:))
    ylabel('$\hat{q}_1$','Interpreter','latex')
    grid minor
    subplot(4,1,2)
    plot(t,EKF.xhat(2,:))
    ylabel('$\hat{q}_2$','Interpreter','latex')
    grid minor
    subplot(4,1,3)
    plot(t,EKF.xhat(3,:))
    ylabel('$\hat{q}_3$','Interpreter','latex')
    grid minor
    subplot(4,1,4)
    plot(t,EKF.xhat(4,:))
    grid minor
    ylabel('$\hat{q}_4$','Interpreter','latex')
    xlabel('Time [s]')
    sgtitle('EKF Quaternion Estimates')

    figure
    subplot(4,1,1)
    plot(t,EKF.delq(1,:))
    ylabel('$\hat{q}_1$','Interpreter','latex')
    grid minor
    subplot(4,1,2)
    plot(t,EKF.delq(2,:))
    ylabel('$\hat{q}_2$','Interpreter','latex')
    grid minor
    subplot(4,1,3)
    plot(t,EKF.delq(3,:))
    ylabel('$\hat{q}_3$','Interpreter','latex')
    grid minor
    subplot(4,1,4)
    plot(t,EKF.delq(4,:))
    grid minor
    ylabel('$\hat{q}_4$','Interpreter','latex')
    xlabel('Time [s]')
    sgtitle('EKF Quaternion Error')

    figure
    subplot(3,1,1)
    plot(t,rad2deg(EKF.xhat(5,:))*3600)
    ylabel('$\hat{\beta}_1$ [deg/hr]','Interpreter','latex')
    grid minor
    subplot(3,1,2)
    plot(t,rad2deg(EKF.xhat(6,:))*3600)
    ylabel('$\hat{\beta}_2$ [deg/hr]','Interpreter','latex')
    grid minor
    subplot(3,1,3)
    plot(t,rad2deg(EKF.xhat(7,:))*3600)
    ylabel('$\hat{\beta}_3$ [deg/hr]','Interpreter','latex')
    grid minor
    xlabel('Time [s]')
    sgtitle('EKF Gyro Bias/Drift Estimates')

    figure
    h1  = subplot(3,1,1);
    plot(t,rad2deg(Gyro.Bias(1,:) - EKF.xhat(5,:))*3600,t,...
    rad2deg(3*sqrt(EKF.Pdiag(4,:)))*3600,'--r',t,rad2deg(-3*sqrt(EKF.Pdiag(4,:)))*3600,'--r')
    ylabel('$\hat{\beta}_1$ Error [deg/hr]','Interpreter','latex')
    legend('Error','Upper 3\sigma Bound','Lower 3\sigma Bound')
    %ylim([-.1 .1])
    grid minor
    h2 = subplot(3,1,2);
    plot(t,rad2deg(Gyro.Bias(2,:) - EKF.xhat(6,:))*3600,t,...
    rad2deg(3*sqrt(EKF.Pdiag(5,:)))*3600,'--r',t,rad2deg(-3*sqrt(EKF.Pdiag(5,:)))*3600,'--r')
    ylabel('$\hat{\beta}_2$ [deg/hr]','Interpreter','latex')
    %ylim([-.1 .1])
    grid minor
    h3  = subplot(3,1,3);
    plot(t,rad2deg(Gyro.Bias(3,:) - EKF.xhat(7,:))*3600,t,...
    rad2deg(3*sqrt(EKF.Pdiag(6,:)))*3600,'--r',t,rad2deg(-3*sqrt(EKF.Pdiag(6,:)))*3600,'--r')
    ylabel('$\hat{\beta}_3$ [deg/hr]','Interpreter','latex')
    grid minor
    xlabel('Time [s]')
    %ylim([-.1 .1])
    sgtitle('EKF Gyro Drift Errors')
    linkaxes([h1,h2,h3])

    figure
    plot(t,vecnorm(EKF.xhat(1:4,:)))
    xlabel('Time [s]')
    ylabel(['Magnitude'])
    title('Magnitude of Quaternion Estimate')

    figure
    plot(t,q,t,EKF.xhat(1:4,:))
    xlabel('Time [s]')
    ylabel('Quaternion')
    legend('$q_1$','$\hat{q}_1$','$q_2$','$\hat{q}_2$','$q_3$','$\hat{q}_3$',...
        '$q_4$','$\hat{q}_4$','Interpreter','latex')
    title('EKF True vs Estimated Quaternion')
    grid minor

    figure
    subplot(3,1,1)
    plot(t,Gyro.Measurements(1,:),t,omegahat(1,:))
    legend('Measurement','Estimate')
    xlabel('Time [s]')
    ylabel('$\omega$ [rad/s]','Interpreter','latex')
    grid minor
    title('EKF Body Angular Velocity')
    subplot(3,1,2)
    plot(t,Gyro.Measurements(2,:),t,omegahat(2,:))
    xlabel('Time [s]')
    ylabel('$\omega$ [rad/s]','Interpreter','latex')
    grid minor
    subplot(3,1,3)
    plot(t,Gyro.Measurements(3,:),t,omegahat(3,:))
    xlabel('Time [s]')
    ylabel('$\omega$ [rad/s]','Interpreter','latex')
    grid minor

    % Zoom

    figure
    ax1 = subplot(3,1,1);
    plot(t,EKF.roll_error,t,rad2deg(3*sqrt(EKF.Pdiag(1,:))),'--r',t,rad2deg(-3*sqrt(EKF.Pdiag(1,:))),'--r')
    ylabel('Roll [deg]')
    ylim([-.01 .01])
    legend('Error','Upper 3\sigma Bound','Lower 3\sigma Bound')
    grid minor
    ax2 = subplot(3,1,2);
    plot(t,EKF.pitch_error,t,rad2deg(3*sqrt(EKF.Pdiag(2,:))),'--r',t,rad2deg(-3*sqrt(EKF.Pdiag(2,:))),'--r')
    ylabel('Pitch [deg]')
    ylim([-.01 .01])
    grid minor
    ax3 = subplot(3,1,3);
    plot(t,EKF.yaw_error,t,rad2deg(3*sqrt(EKF.Pdiag(3,:))),'--r',t,rad2deg(-3*sqrt(EKF.Pdiag(3,:))),'--r')
    ylim([-.01 .01])
    ylabel('Yaw [deg]')
    grid minor
    xlabel('Time [s]')
    sgtitle('EKF Attitude Errors')
    linkaxes([ax1,ax2,ax3]);


    
    %% Plot Results
    figure
    subplot(4,1,1)
    plot(t,UKF.qhat(1,:))
    ylabel('$\hat{q}_1$','Interpreter','latex')
    grid minor
    subplot(4,1,2)
    plot(t,UKF.qhat(2,:))
    ylabel('$\hat{q}_2$','Interpreter','latex')
    grid minor
    subplot(4,1,3)
    plot(t,UKF.qhat(3,:))
    ylabel('$\hat{q}_3$','Interpreter','latex')
    grid minor
    subplot(4,1,4)
    plot(t,UKF.qhat(4,:))
    grid minor
    ylabel('$\hat{q}_4$','Interpreter','latex')
    xlabel('Time [s]')
    sgtitle('UKF Quaternion Estimates')

    figure
    subplot(4,1,1)
    plot(t,delqtrue(1,:))
    ylabel('$\hat{q}_1$','Interpreter','latex')
    grid minor
    subplot(4,1,2)
    plot(t,delqtrue(2,:))
    ylabel('$\hat{q}_2$','Interpreter','latex')
    grid minor
    subplot(4,1,3)
    plot(t,delqtrue(3,:))
    ylabel('$\hat{q}_3$','Interpreter','latex')
    grid minor
    subplot(4,1,4)
    plot(t,delqtrue(4,:))
    grid minor
    ylabel('$\hat{q}_4$','Interpreter','latex')
    xlabel('Time [s]')
    sgtitle('UKF Quaternion Error')

    figure
    f1 =  subplot(3,1,1);
    plot(t,UKF.roll_error,t,rad2deg(3*sqrt(UKF.Pdiag(1,:))),'--r',t,rad2deg(-3*sqrt(UKF.Pdiag(1,:))),'--r')
    ylabel('Roll [deg]')
    legend('Error','Upper 3\sigma Bound','Lower 3\sigma Bound')
    grid minor
    f2  = subplot(3,1,2);
    plot(t,UKF.pitch_error,t,rad2deg(3*sqrt(UKF.Pdiag(2,:))),'--r',t,rad2deg(-3*sqrt(UKF.Pdiag(2,:))),'--r')
    ylabel('Pitch [deg]')
    grid minor
    f3 = subplot(3,1,3);
    plot(t,UKF.yaw_error,t,rad2deg(3*sqrt(UKF.Pdiag(3,:))),'--r',t,rad2deg(-3*sqrt(UKF.Pdiag(3,:))),'--r')
    ylabel('Yaw [deg]')
    grid minor
    xlabel('Time [s]')
    sgtitle('UKF Attitude Errors')
    linkaxes([f1,f2,f3])

    figure
    subplot(3,1,1)
    plot(t,rad2deg(UKF.xhat(4,:))*3600)
    ylabel('$\hat{\beta}_1$ [deg/hr]','Interpreter','latex')
    grid minor
    subplot(3,1,2)
    plot(t,rad2deg(UKF.xhat(5,:))*3600)
    ylabel('$\hat{\beta}_2$ [deg/hr]','Interpreter','latex')
    grid minor
    subplot(3,1,3)
    plot(t,rad2deg(UKF.xhat(6,:))*3600)
    ylabel('$\hat{\beta}_3$ [deg/hr]','Interpreter','latex')
    grid minor
    xlabel('Time [s]')
    sgtitle('UKF Gyro Bias/Drift Estimates')

    figure
    H1 = subplot(3,1,1);
    plot(t,rad2deg(Gyro.Bias(1,:) - UKF.xhat(4,:))*3600,t,...
    rad2deg(3*sqrt(UKF.Pdiag(4,:)))*3600,'--r',t,rad2deg(-3*sqrt(UKF.Pdiag(4,:)))*3600,'--r')
    ylabel('$\hat{\beta}_1$ Error [deg/hr]','Interpreter','latex')
    legend('Error','Upper 3\sigma Bound','Lower 3\sigma Bound')
    grid minor
    H2 = subplot(3,1,2);
    plot(t,rad2deg(Gyro.Bias(2,:) - UKF.xhat(5,:))*3600,t,...
    rad2deg(3*sqrt(UKF.Pdiag(5,:)))*3600,'--r',t,rad2deg(-3*sqrt(UKF.Pdiag(5,:)))*3600,'--r')
    ylabel('$\hat{\beta}_2$ [deg/hr]','Interpreter','latex')
    grid minor
    H3 = subplot(3,1,3);
    plot(t,rad2deg(Gyro.Bias(3,:) - UKF.xhat(6,:))*3600,t,...
    rad2deg(3*sqrt(UKF.Pdiag(6,:)))*3600,'--r',t,rad2deg(-3*sqrt(UKF.Pdiag(6,:)))*3600,'--r')
    ylabel('$\hat{\beta}_3$ [deg/hr]','Interpreter','latex')
    grid minor
    xlabel('Time [s]')
    sgtitle('UKF Gyro Drift Errors')
    linkaxes([H1,H2,H3]);

    figure
    plot(t,q,t,UKF.qhat)
    xlabel('Time [s]')
    ylabel('Quaternion')
    legend('$q_1$','$\hat{q}_1$','$q_2$','$\hat{q}_2$','$q_3$','$\hat{q}_3$',...
        '$q_4$','$\hat{q}_4$','Interpreter','latex')
    title('UKF True vs Estimated Quaternion')
    grid minor
    
    figure
    plot(t,vecnorm(UKF.qhat(1:4,:)))
    xlabel('Time [s]')
    ylabel(['Magnitude'])
    title('Magnitude of Quaternion Estimate')
    
    %% Calculate Mean Error Stats
    EKF.MeanYawError    = abs(mean(EKF.yaw_error));
    EKF.MeanRollError   = abs(mean(EKF.roll_error));
    EKF.MeanPitchError  = abs(mean(EKF.pitch_error));
    EKF.MeanBiasError   = abs(mean(rad2deg(Gyro.Bias - EKF.xhat(5:7,:))*3600,2));
    
    UKF.MeanYawError    = abs(mean(UKF.yaw_error));
    UKF.MeanRollError   = abs(mean(UKF.roll_error));
    UKF.MeanPitchError  = abs(mean(UKF.pitch_error));
    UKF.MeanBiasError   = abs(mean(rad2deg(Gyro.Bias - UKF.xhat(4:6,:))*3600,2));
    
    ParameterNames      = {'Run Time [s]','Mean Yaw Error [deg]','Mean Roll Error [deg]','Mean Pitch Error [deg]'...
                          ,'Gyro Bias 1 Error [deg/hr]','Gyro Bias 2 Error [deg/hr]','Gyro Bias 3 Error [deg/hr]'};
    
    TableData           = [EKF.RunTime, UKF.RunTime; EKF.MeanYawError UKF.MeanYawError; EKF.MeanPitchError UKF.MeanPitchError;...
                           EKF.MeanRollError, UKF.MeanRollError; EKF.MeanBiasError, UKF.MeanBiasError];

    Table                           = array2table(TableData);
    Table.Properties.VariableNames  = {'EKF','UKF'};
    Table.Properties.RowNames       = ParameterNames;
    writetable(Table,filename,'Sheet','SingleRun','WriteRowNames',true)
    
    
else
    
    %% Monte Carlo
    
    % Create Random Initial Conditions
    std         = .05;
    Bhat0_vec   = [0 0 0]' + std*abs(randn(3,numruns));
    qhat0_vec   = (sqrt(2)/2)*[0, 1, 1, 0]' + std*abs(randn(4,numruns));
    % Initial attitude covaraince - 1 deg^2 
    P0a             = 1 * (pi/180)^2 *ones(3,1); 

    % Initial bias covariance - .2 (deg/hr)^2
    P0b             = .2 *(pi/(180*3600))^2 *ones(3,1);
    
    for MC = 1:length(qhat0_vec)
        qhat0_vec(:,MC) = qhat0_vec(:,MC)/norm(qhat0_vec(:,MC));   
    end
    
    EKF.RunTime     = zeros(1,numruns);
    UKF.RunTime     = zeros(1,numruns);
    EKF.AttError    = zeros(1,numruns);
    EKF.BiasError   = zeros(1,numruns);

    UKF.AttError    = zeros(1,numruns);
    UKF.BiasError   = zeros(1,numruns);
    
    for MC = 1:length(qhat0_vec)
        qhat0   = qhat0_vec(:,MC);
        Bhat0   = Bhat0_vec(:,MC);
        
        tic
        MEKF
        EKF.RunTime(MC) = toc;
    
        pause(1);
    
        tic
        USQUE
        UKF.RunTime(MC) = toc;
        
        att_err                 = [EKF.yaw_error;EKF.roll_error;EKF.pitch_error];
        EKF.AttError(MC)        = sum(vecnorm(att_err));
        EKF.BiasError(MC)       = rad2deg(sum(vecnorm(Gyro.Bias - EKF.xhat(5:7,:))))*3600;

        att_err                 = [UKF.yaw_error;UKF.roll_error;UKF.pitch_error];
        UKF.AttError(MC)        = sum(vecnorm(att_err));
        UKF.BiasError(MC)       = rad2deg(sum(vecnorm(Gyro.Bias - UKF.xhat(4:6,:))))*3600;

         
    end

    figure
    plot((1:numruns),EKF.AttError,'*',(1:numruns),UKF.AttError,'o')
    ylabel('Sum of Norm of Attitude Errors [deg]');
    title('Attitude Errors')
    legend('MEKF','USQUE')
    grid minor
    xlabel('Num Runs')
   

    figure
    plot((1:numruns),EKF.BiasError,'*',(1:numruns),UKF.BiasError,'o')
    ylabel('Sum of Norm of Bias Errors [deg/hr]');
    title('Bias Errors')
    legend('MEKF','USQUE')
    grid minor
    xlabel('Num Runs')
    
    EKF.MeanAttError = mean(EKF.AttError);
    UKF.MeanAttError = mean(UKF.AttError);

end

