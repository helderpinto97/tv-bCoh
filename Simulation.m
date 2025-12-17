%% test IIR and total interdependence measures on a 7-variate VAR model
clear; close all; clc;
addpath('functions');

c1max = 0.8;
fs    = 100;
nfft  = 1000;

f_osc = 0.25;
nobs  = 1000;
k     = 1:nobs;
t     = k/fs;
DC    = 50;

c1 = -c1max*square(2*pi*f_osc*t,DC);
c1(c1==-c1max)=0;
c2 = 1-c1; 

M = 7;
c_forget_factor = 0.98;

par.poles{1}=([0.9 0.35]);
par.poles{4}=([0.9 0.1]);
par.poles{6}=([0.9 0.1]);
par.poles{7}=[];

%% 1) TRUE time-varying VAR
for n=1:length(t)
    C = c1(n);
    par.coup=[1 2 1 0.8;
              1 3 1 C;
              4 5 2 C;
              6 5 1 0.8;
              5 7 2 C;
              2 7 2 0.8];

    par.Su = ones(1,M);
    [Am1,Su1] = theoreticalVAR(M,par);

    Am1 = Am1';
    Am1 = reshape(Am1,[M,M,size(Am1,2)/M]);

    AmT(:,:,:,n) = Am1;
    SuT(:,:,n)   = Su1;
end

p = size(AmT,3);

for n=1:nobs
    SuT(:,:,n) = (SuT(:,:,n)+SuT(:,:,n)')/2 + 1e-12*eye(M);
end

%% 2) Generate realization (VAR filter)
Xgen = var_to_tsdata_nonstat(AmT,SuT,1);
Xsim = Xgen(:,:,1);

X = permute(Xsim,[1 3 2]);
X = zscore(X,0,3); % Normalize data

%% 3) TV-VAR estimation
[A_i,Su_i] = idMVAR(squeeze(X),p,0);
[Am_n_rls,Su_n_rls] = tvID_VAR_RLS_IC(X,p,c_forget_factor,A_i,Su_i);

%% 4) Time-varying spectra
S_T = complex(zeros(M,M,nfft,nobs));
for n=1:nobs
    [S_T(:,:,:,n),~,f] = sir_VARspectra(Am_n_rls(:,:,n),Su_n_rls(:,:,n),nfft,fs);
end

%% 5) COMPLETE NETWORK
for tt=1:size(S_T,4)
    out=iispectral(S_T(:,:,:,tt),7,[1 2 3],[4 5 6]);
    iy_x1_x2(:,tt)=out.iy_x1_x2;
    fy_x1x2(:,tt)=out.fy_x1x2;
    fy_x1(:,tt)=out.fy_x1;
    fy_x2(:,tt)=out.fy_x2;
    Fy_x1x2(tt)=out.Fy_x1x2;
    Fy_x1(tt)=out.Fy_x1;
    Fy_x2(tt)=out.Fy_x2;
    Iy_x1_x2(tt)=out.Iy_x1_x2;
end

figure('WindowState','maximized');
subplot(2,4,4), plot(t,Iy_x1_x2), ylim([-1 1]), title('I_{Y;X_1;X_2}');
subplot(2,4,3), plot(t,Fy_x1x2), ylim([0 2]), title('F_{Y;X_1,X_2}');
subplot(2,4,1), plot(t,Fy_x1),   ylim([0 2]), title('F_{Y;X_1}');
subplot(2,4,2), plot(t,Fy_x2),   ylim([0 2]), title('F_{Y;X_2}');

t_in = p+5;
climF = [-5 5];
climI = [-1 1];

subplot_tight(2,4,5), imagesc(t(t_in:end),f,fy_x1(:,t_in:end),climF); axis tight; colorbar
subplot_tight(2,4,6), imagesc(t(t_in:end),f,fy_x2(:,t_in:end),climF); axis tight; colorbar
subplot_tight(2,4,7), imagesc(t(t_in:end),f,fy_x1x2(:,t_in:end),climF); axis tight; colorbar
subplot_tight(2,4,8), imagesc(t(t_in:end),f,iy_x1_x2(:,t_in:end),climI); axis tight; colorbar

%% FIRST BLOCK X1
for tt=1:size(S_T,4)
    out=iispectral(S_T(1:3,1:3,:,tt),1,2,3);
    iy_x1_x2(:,tt)=out.iy_x1_x2;
    fy_x1x2(:,tt)=out.fy_x1x2;
    fy_x1(:,tt)=out.fy_x1;
    fy_x2(:,tt)=out.fy_x2;
    Fy_x1x2(tt)=out.Fy_x1x2;
    Fy_x1(tt)=out.Fy_x1;
    Fy_x2(tt)=out.Fy_x2;
    Iy_x1_x2(tt)=out.Iy_x1_x2;
end

figure('WindowState','maximized');
subplot(2,4,4), plot(t,Iy_x1_x2), ylim([0 1.5])
subplot(2,4,3), plot(t,Fy_x1x2), ylim([0 1.5])
subplot(2,4,1), plot(t,Fy_x1),   ylim([0 1.5])
subplot(2,4,2), plot(t,Fy_x2),   ylim([0 1.5])

subplot(2,4,5), imagesc(t(t_in:end),f,fy_x1(:,t_in:end),climF); axis tight; colorbar
subplot(2,4,6), imagesc(t(t_in:end),f,fy_x2(:,t_in:end),climF); axis tight; colorbar
subplot(2,4,7), imagesc(t(t_in:end),f,fy_x1x2(:,t_in:end),climF); axis tight; colorbar
subplot(2,4,8), imagesc(t(t_in:end),f,iy_x1_x2(:,t_in:end),climI); axis tight; colorbar

%% SECOND BLOCK X2
for tt=1:size(S_T,4)
    out=iispectral(S_T(4:6,4:6,:,tt),2,1,3);
    iy_x1_x2(:,tt)=out.iy_x1_x2;
    fy_x1x2(:,tt)=out.fy_x1x2;
    fy_x1(:,tt)=out.fy_x1;
    fy_x2(:,tt)=out.fy_x2;
    Fy_x1x2(tt)=out.Fy_x1x2;
    Fy_x1(tt)=out.Fy_x1;
    Fy_x2(tt)=out.Fy_x2;
    Iy_x1_x2(tt)=out.Iy_x1_x2;
end

figure('WindowState','maximized');
subplot(2,4,4), plot(t,Iy_x1_x2), ylim([-1 1])
subplot(2,4,3), plot(t,Fy_x1x2), ylim([0 1.5])
subplot(2,4,1), plot(t,Fy_x1),   ylim([0 1.5])
subplot(2,4,2), plot(t,Fy_x2),   ylim([0 1.5])

subplot(2,4,5), imagesc(t(t_in:end),f,fy_x1(:,t_in:end),climF); axis tight; colorbar
subplot(2,4,6), imagesc(t(t_in:end),f,fy_x2(:,t_in:end),climF); axis tight; colorbar
subplot(2,4,7), imagesc(t(t_in:end),f,fy_x1x2(:,t_in:end),climF); axis tight; colorbar
subplot(2,4,8), imagesc(t(t_in:end),f,iy_x1_x2(:,t_in:end),climI); axis tight; colorbar
