clear
clc
%% Data obtained from Chilbolton Site
% These tau values are obtained using the exponential fitting of the
% autocorrelation function exp(-t/tau). Using this model, we can obtain the
% relationshship between tau and the rain
Rs  = [0.6023 4.9970 8.9893 12.6261 16.2969 20.6991 28.2095 36.1591] ;      % Rain data binned as per chilbolton data [mm/hr]
TU  = [3.9706 4.2180 3.2370  6.7748  3.6766  5.0195  1.1731  0.0200] ;      % Tau correlated as per chilbolton data [seconds]
P   = polyfit(Rs,TU,1) ;                                                    % Linear fitting
%% Call model to obtain statistics (mu, sigma, theta)
% We can obtain the average and standard deviation of the PathLoss using
% the statistical model based on the measurements obtained
% Mu    : Model PathLoss mean of the rainrate and elevation angle
% Sigma : Model PathLoss std of the rainrate and elevation angle
[Mu,Sigma,~] = ModelFitting(0) ;
%% Satellite pass, elevation angle, and rain profile
% Load a satellite pass which has the elevation angle as a function of time
% and define the rain rate as a function of time for that pass as well as
% the sampling time
% R     : Expected rainrate as a function of time [mm/hr] (Forecast data)
% Ts    : Sampling Time
% Theta : Elevation angle as a function of time [Degrees]
load('SatellitePass.mat')
R       = 5*ones(size(Theta)) ;
Ts      = 0.1 ;
%% Use Correlation model and filter to obtain the pathloss
% y     : ModelingPathLoss as function of time
% tau   : Correlation time as a function of time
tau     = polyval(P,R) ;
for Ctr = 1 : length(tau)
    if isnan(Theta(Ctr))
        y(Ctr) = NaN ;
    else
        y(Ctr) = CorrelationModel(tau(Ctr),Mu(R(Ctr),Theta(Ctr)),Sigma(R(Ctr),Theta(Ctr)),Ts) ;
    end
end
%% Plot the results
% close all
Scale = 1 ;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[1000 150 800 800/1.618*Scale]);
set(gca,'Gridlinestyle','--') ;
ax = gca ;
hold on
yyaxis left
plot(Time,y,'.')
xlabel('Time, t [s]')
ylabel('Excess Path-loss, \eta [dB]')
ylim([-10 100])
box on
grid on
yyaxis right
plot(Time,Theta,'.')
xlabel('Time, t [s]')
ylabel('Elevation angle, \theta [dB]')
box on
grid on
ax.FontSize = 14 ;
ax.FontName = 'Times New Roman' ;
% 
% % Saving figure
% Filename='..\..\..\LatexSource\Figures\ExcessCMPrediction';
% print( '-depsc','-r600',Filename);
% print( '-dpng','-r600',Filename);
