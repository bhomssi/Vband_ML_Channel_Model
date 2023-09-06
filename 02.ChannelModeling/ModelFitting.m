function [Mu,Sigma,Theta] = ModelFitting(newData)
% This function provides the statistical model to drive the prediction
% You can use the current model, or can be updated using a new set of data
% The new data needs to have the PathLoss mean (AttMu), the PathLoss
% standard deviaton (AttSg), the parameters of the measurements and the
% RinRate in mm/hr for the mean and standard deviation
% newData:          Switch for new data (1) or old data (0)
% Mu:               The average obtained using the model [dB]
% Sigma :           The standard deviation obtained using the model [dB]
% ModelElevationAngle: The Elevation angle used in the modeling [Degrees]
if newData == 0
    load('BinnedData_chilbolton.mat')
elseif newData == 1
    load('NewData.mat')
end
%% ITU Model for Path-loss
aj      = [-5.33980,-0.35351,-0.23789,-0.94158] ;
bj      = [-0.10008,1.26970,0.86036,0.64552] ;
cj      = [1.13098,0.45400,0.15354,0.16817] ;
mk      = -0.18961 ;
ck      = 0.71147 ;
kH      = 10.^(sum(aj.*exp(-((log10(par.Freq'/1e9) - bj)./cj).^2),2) + mk.*log10(par.Freq'/1e9) + ck) ;
aj      = [-3.80595,-3.44965,-0.39902,0.50167] ;
bj      = [0.56934,-0.22911,0.73042,1.07319] ;
cj      = [0.81061,0.51059,0.11899,0.27195] ;
mk      = -0.16398 ;
ck      = 0.63297 ;
kV      = 10.^(sum(aj.*exp(-((log10(par.Freq'/1e9) - bj)./cj).^2),2) + mk.*log10(par.Freq'/1e9) + ck) ;
aj      = [-0.14318 0.29591 0.32177 -5.37610 16.1721] ;
bj      = [1.82442 0.77564 0.63773 -0.96230 -3.29980] ;
cj      = [-0.55187 0.19822 0.13164 1.47828 3.43990] ;
mk      = 0.67849 ;
ck      = -1.95537 ;
alphaH  = sum(aj.*exp(-((log10(par.Freq'/1e9) - bj)./cj).^2),2) + mk.*log10(par.Freq'/1e9) + ck ;
aj      = [-0.07771 0.56727 -0.20238 -48.2991 48.5833] ;
bj      = [2.33840 0.95545 1.14520 0.791669 0.791459] ;
cj      = [-0.76284 0.54039 0.26809 0.116226 0.116479] ;
mk      = -0.053739 ;
ck      = 0.83433 ;
alphaV  = sum(aj.*exp(-((log10(par.Freq'/1e9) - bj)./cj).^2),2) + mk.*log10(par.Freq'/1e9) + ck ;

k       = @(theta) (kH + kV + (kH - kV).*cosd(theta).^2.*cosd(2*45))/2 ;
alpha   = @(theta) (kH.*alphaH + kV.*alphaV + (kH.*alphaH - kV.*alphaV).*cosd(theta).^2.*cosd(2.*45))./(2.*k(theta)) ;
%% Developed Model for the Average and Stndard Deviation
Ls      = @(theta) (par.hR - 0.1)./sind(theta) ;
LG      = @(theta) (par.hR - 0.1)./sind(theta).*cosd(theta) ;
%% Model Fitting for Mu and Sigma
% fitfunM = fittype( @(a,x,y) k(y).*x.^alpha(y).*Ls(y).*a,'independent', {'x', 'y'},'dependent', 'z') ;
% Mu      = fit([Rbins' par.El*ones(size(Rbins'))],AttMu',fitfunM,'StartPoint',1)
% fitfunS = fittype( @(a,b,x,y) a.*x.^b.*LG(y),'independent', {'x', 'y'},'dependent', 'z') ;
% Sigma   = fit([Rbins' par.El*ones(size(Rbins'))],AttSg',fitfunS,'StartPoint',[1 1])
Theta   = par.El ;
Mu      = @(R,theta) k(theta).*R.^alpha(theta).*Ls(theta).*0.4354 ;
Sigma   = @(R,theta) 0.2962.*R.^0.6732.*LG(theta) ;
end