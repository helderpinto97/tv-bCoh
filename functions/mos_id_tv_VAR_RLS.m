%% Model Order Selection for identification time-varying VAR model
% The model order selection follow the approach introduced in [1].
% [1]- Costa, Antonio H., and Stephan Hengstler. "Adaptive time–frequency 
% analysis based on autoregressive modeling." Signal Processing (2011)"
%%% input:
% Y, R x M x N matrix of time series 
% pmax, maximum tested model order
% forgetting factor for RLS

%%% output:
% popt_RLS: model order optimized with  RLS
% 

function [popt_RLS] = mos_id_tv_VAR_RLS(Y,pmax,c)

for p=1:pmax
    % RLS
    [~,~,~,Z_RLS(p,:)]=tvID_VAR_RLS(Y,p,c);        
end

Z_RLS=Z_RLS(:,p+1:end);
Z_RLS1=mean(Z_RLS,2);
[popt_RLS] = knee_pt(Z_RLS1);
