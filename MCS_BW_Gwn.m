%% Monte-Carlo simulation of MDOF Bouc-Wen model enforced by Gaussian white noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %                                       
% Copyright (C) 2023                                                      %
% Jian-Bing Chen, Meng-Ze Lyu, Jia-Shu Yang, Ting-Ting Sun, and Yi Luo    %
% College of Civiling  Engineering, Tongji Unversity, Shanghai, China.    %
% All Rights Reserved.                                                    %
% You may use, distribute and modify this code under the terms of the     %
% license file "LICENSE.txt"                                              %
%                                                                         %                                                                       %
% Disclaimer:                                                             %
% The authors reserve all rights but do not guarantee that the code is    %
% free from errors. Furthermore, The authors shall not be not responsible %
% for the correctness of the results derived by the program.              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc,clear
%% Input of structural parameters
dof=10;             % Degrees of freedom
l_itrst=10;         % DOF of interest
str_colnum=3;                                               % Number of column of each story
str_m=ones(dof,1)*2.6e5;                                    % Lumped mass of each story (kg) [dof by 1 vector]
str_h=[4,3,3,3,3,3,3,3,3,3]';                               % Height of each story (m) [dof by 1 vector]
str_a=[0.5,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4]';           % Length of column section of each story (m) [dof by 1 vector]
str_b=[0.5,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4]';           % Width of column section of each story (m) [dof by 1 vector]
str_E=ones(dof,1)*3e10;                                     % Initial Young's modulus of each story (Pa) [dof by 1 vector]
zeta=[0.05,0.05]';                                          % Damping ratios of the first two modes [2 by 1 vector]
bwm_alpha=0.04;         % Ratio of linear to nonlinear response
bwm_A=1;                % Basic hysteresis shape controls of Bouc-Wen model
bwm_n=1;               	% Sharpness of yield of Bouc-Wen model
bwm_beta=15;        	% Basic hysteresis shape controls of Bouc-Wen model (/m)
bwm_gamma=150;        	% Basic hysteresis shape controls of Bouc-Wen model (/m)
bwm_d_nu=1000;       	% Strength degradation of Bouc-Wen model (/m2)
bwm_d_eta=1000;         % Stiffness degradation of Bouc-Wen model (/m2)
bwm_p=1000;           	% Pinching slope of Bouc-Wen model (/m2)
bwm_q=0.25;            	% Pinching initiation of Bouc-Wen model
bwm_d_psi=5;          	% Pinching rate of Bouc-Wen model (/m2)
bwm_lambda=0.5;       	% Pinching severity/rate interaction of Bouc-Wen model
bwm_zeta_s=0.99;        % Measure of total slip of Bouc-Wen model
bwm_psi=0.05;           % Pinching magnitude of Bouc-Wen model (m)
D=0.0072;                 	% Intensity of Gaussian white noise (m2/s3)
x0=zeros(4*dof,1);          % Initial values of displacements, velocitys, hysteretic displacements, and hysteretic-dissipated energies [4*dof by 1 vector]
%% Input of numerical parameters
T=10;                   % Time domain (s)
dt=0.005;               % Time step (s)
N=1e6;                  % Number of MCS
Nrec=10;                % Number for records
Ttyp=[4,7,10];          % Typical instants (s)
threshold=0.1:0.02:0.3;         % Thresholds fot Xl (m)
N_bar=50;           % Number of bars for histogram
%% Calculation of structural property
bwm=[bwm_alpha,bwm_A,bwm_n,bwm_beta,bwm_gamma,bwm_d_nu,bwm_d_eta,bwm_p,bwm_q,bwm_d_psi,bwm_lambda,bwm_zeta_s,bwm_psi]';         % Parameters of Bouc-Wen model [13 by 1 vector]
str_I=str_b.*str_a.^3/12;                                                           % Moment of inertia of column section of each story (m4) [dof by 1 vector]
str_k=str_colnum*12*str_E.*str_I./str_h.^3;                                         % Initial lateral inter-story stiffness of each story (N/m) [dof by 1 vector]
M=diag(str_m);                                                                      % Lumped mass (kg) [dof by dof matrix]
invM=diag(1./str_m);                                                                % Inverse of lumped mass (/kg) [dof by dof matrix]
K=diag(str_k+[str_k(2:end);0])-diag(str_k(2:end),1)-diag(str_k(2:end),-1);          % Initial lateral inter-story stiffness (N/m) [dof by dof matrix]
eig_val=eig(M\K);                                                                   % Eigenvalues (/s2) [dof by 1 vector]
omega=sort(sqrt(eig_val));                                                          % Circular frequency (/s) [dof by 1 vector]
dmp_b=2*(zeta(1)*omega(1)-zeta(2)*omega(2))/(omega(1)^2-omega(2)^2);                % Coefficient of Rayleigh damping (s)
dmp_a=2*zeta(1)*omega(1)-dmp_b*omega(1)^2;                                          % Coefficient of Rayleigh damping (/s)
C=dmp_a*M+dmp_b*K;                                                                  % Rayleigh damping (kg/s) [dof by dof matrix]
b=[zeros(dof,1);ones(dof,1);zeros(2*dof,1)];                                        % Load position [4*dof by 1 vector]
%% Treatment of numerical parameters
Nt=T/dt;                        % Number of time steps
t=0:dt:T;                       % Time series (s) [(Nt+1) by 1 vector]
%% MCS by 2-order stochastic Runge-Kutta algorithm
Dhist=zeros(Nrec,Nt+1);                         % Pre-allocated storage for histories of inter-story drift of the bottom story [Nrec by (Nt+1) matrix]
Fhist=zeros(Nrec,Nt+1);                         % Pre-allocated storage for histories of inter-story restoring force of the bottom story [Nrec by (Nt+1) matrix]
Xsample=zeros(N,length(Ttyp));                  % Pre-allocated storage for samples of Xl at typical instants [N by length(Ttyp) matrix]
Vsample=zeros(N,length(Ttyp));                  % Pre-allocated storage for samples of Vl at typical instants [N by length(Ttyp) matrix]
varX_mcs=zeros(1,Nt+1);                         % Pre-allocated storage for variance of Xl [1 by (Nt+1) vector]    
varV_mcs=zeros(1,Nt+1);                         % Pre-allocated storage for variance of Vl [1 by (Nt+1) vector] 
kurtX_mcs=3*ones(1,Nt+1);                       % Pre-allocated storage for kurtosis of Xl [1 by (Nt+1) vector] 
kurtV_mcs=3*ones(1,Nt+1);                       % Pre-allocated storage for kurtosis of Vl [1 by (Nt+1) vector] 
Pf_mcs=zeros(length(threshold),Nt+1);           % Pre-allocated storage for failure probability [length(threshold) by (Nt+1) matrix] 
Xtemp=x0*ones(1,N);                     % Initial value for all samples [4*dof by N matrix]
Ztemp=x0(l_itrst)*ones(N,1);            % Initial value for extreme value samples of Xl [N by 1 vector]
for j=1:Nt
    dW=sqrt(D*dt)*randn(1,N);                               % Wiener increments (m/s) [1 by N vector]
    F1=f_BoucWen(Xtemp,M,C,str_k,bwm);                      % 2-order stochastic Runge-Kutta algorithm [4*dof by N matrix] {call function 'f_BoucWen'}
    F2=f_BoucWen(Xtemp+F1*dt+b*dW,M,C,str_k,bwm);           % 2-order stochastic Runge-Kutta algorithm [4*dof by N matrix] {call function 'f_BoucWen'}
    dX=(F1+F2)/2*dt+b*dW;                                   % Responses increments [4*dof by N matrix]
    Xtemp=Xtemp+dX;                                         % Responses at the j-th step [4*dof by N matrix]
    Xitrst=Xtemp(l_itrst,:)';                               % Samples of Xl (m) [N by 1 vector]
    Vitrst=Xtemp(dof+l_itrst,:)';                           % Samples of Vl (m/s) [N by 1 vector]
    DX=Xtemp(1:dof,:)-[zeros(1,N);Xtemp(1:dof-1,:)];                        % Inter-story drift of each story at the j-th step (m) [dof by N matrix]
    DG=str_k.*(bwm_alpha*DX+(1-bwm_alpha)*Xtemp(2*dof+1:3*dof,:));          % Inter-story restoring force of each story at the j-th step (N) [dof by N matrix]
    Dhist(:,j+1)=Xtemp(1,1:Nrec)';                                          % Histories of inter-story drift of the bottom story (m)
    Fhist(:,j+1)=DG(1,1:Nrec)';                                             % Histories of inter-story restoring force of the bottom story (N)
    for i=1:length(Ttyp)
        if j==Ttyp(i)/dt
            Xsample(:,i)=Xitrst;            % Samples of Xl at typical instants (m)
            Vsample(:,i)=Vitrst;            % Samples of Vl at typical instants (m/s)
        end
    end
    varX_mcs(j+1)=sum(Xitrst.^2)/N;                             % Variance of Xl (m2)
    varV_mcs(j+1)=sum(Vitrst.^2)/N;                             % Variance of Vl (m2/s2)
    kurtX_mcs(j+1)=sum(Xitrst.^4)/N/varX_mcs(j+1)^2;            % Kurtosis of Xl
    kurtV_mcs(j+1)=sum(Vitrst.^4)/N/varV_mcs(j+1)^2;            % Kurtosis of Vl
    Ztemp=max([abs(Xitrst),Ztemp],[],2);                            % Extreme value samples of Xl (m) [N by 1 vector]
    for i=1:length(threshold)
        Pf_mcs(i,j+1)=length(find(Ztemp>threshold(i)))/N;           % Failure probability
    end
    fprintf('MCS at the %i-th step is in progress\n',j);
end
%% Post-processing for histogram
x_bar=zeros(N_bar+1,length(Ttyp));              % Pre-allocated storage for grid sequence of Xl [(Nbar+1) by length(Ttyp) matrix]
v_bar=zeros(N_bar+1,length(Ttyp));              % Pre-allocated storage for grid sequence of Vl [(Nbar+1) by length(Ttyp) matrix]
PDF_xl_mcs=zeros(N_bar,length(Ttyp));           % Pre-allocated storage for PDF of Xl [Nbar by length(Ttyp) matrix]
PDF_vl_mcs=zeros(N_bar,length(Ttyp));           % Pre-allocated storage for PDF of Vl [Nbar by length(Ttyp) matrix]
for i=1:length(Ttyp)
    bar_infX=min(Xsample(:,i))*(1+0.1/(N_bar-0.2))-max(Xsample(:,i))*0.1/(N_bar-0.2);           % Inf of grid sequence of Xl (m)
    bar_supX=max(Xsample(:,i))*(1+0.1/(N_bar-0.2))-min(Xsample(:,i))*0.1/(N_bar-0.2);           % Sup of grid sequence of Xl (m)
    d_barX=(bar_supX-bar_infX)/N_bar;                                                           % Grid size of Xl (m)
    x_bar(:,i)=(bar_infX:d_barX:bar_supX)';                                                     % Grid sequence of Xl (m)
    for k1=1:N_bar
        PDF_xl_mcs(k1,i)=length(find(Xsample(:,i)>x_bar(k1,i) & Xsample(:,i)<=x_bar(k1+1,i)))/(N*d_barX);          % PDF of Xl (/m)
    end
    bar_infV=min(Vsample(:,i))*(1+0.1/(N_bar-0.2))-max(Vsample(:,i))*0.1/(N_bar-0.2);           % Inf of grid sequence of Vl (m/s)
    bar_supV=max(Vsample(:,i))*(1+0.1/(N_bar-0.2))-min(Vsample(:,i))*0.1/(N_bar-0.2);           % Sup of grid sequence of Vl (m/s)
    d_barV=(bar_supV-bar_infV)/N_bar;                                                           % Grid size of Vl (m/s)
    v_bar(:,i)=(bar_infV:d_barV:bar_supV)';                                                     % Grid sequence of Vl (m/s)
    for k2=1:N_bar
        PDF_vl_mcs(k2,i)=length(find(Vsample(:,i)>v_bar(k2,i) & Vsample(:,i)<=v_bar(k2+1,i)))/(N*d_barV);          % PDF of Vl (s/m)
    end
end
save('MCS_BW_Gwn.mat','Xsample','Vsample','varX_mcs','varV_mcs','kurtX_mcs','kurtV_mcs','PDF_xl_mcs','PDF_vl_mcs','Pf_mcs','x_bar','v_bar','Ttyp','threshold','N','N_bar');
%% Figures
figure(1),hold on           % Inter-story restoring force versus drift of the bottom story for a typical sample
    i_rec=1;
    plot(Dhist(i_rec,:),Fhist(i_rec,:),'k','LineWidth',0.5)
    xlabel('Drift (m)'),ylabel('Restoring force (N)')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);
figure(2),hold on           % Histogram PDF of Xl at a typical instant
    i_typ=3;
    histogX=reshape([x_bar(1:end-1,i_typ)';x_bar(1:end-1,i_typ)';x_bar(2:end,i_typ)';x_bar(2:end,i_typ)'],4*N_bar,1);
    histogPx=reshape([ones(1,N_bar)*1e-10;PDF_xl_mcs(:,i_typ)';PDF_xl_mcs(:,i_typ)';ones(1,N_bar)*1e-10],4*N_bar,1);
    plot(histogX,histogPx,'b','LineWidth',0.5)
    xlabel('\itx\rm_{10} (m)'),ylabel('PDF (m^{-1})')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);
figure(3),hold on           % Histogram PDF of Vl at a typical instant
    i_typ=3;
    histogV=reshape([v_bar(1:end-1,i_typ)';v_bar(1:end-1,i_typ)';v_bar(2:end,i_typ)';v_bar(2:end,i_typ)'],4*N_bar,1);
    histogPv=reshape([ones(1,N_bar)*1e-10;PDF_vl_mcs(:,i_typ)';PDF_vl_mcs(:,i_typ)';ones(1,N_bar)*1e-10],4*N_bar,1);
    plot(histogV,histogPv,'b','LineWidth',0.5)
    xlabel('\itv\rm_{10} (m/s)'),ylabel('PDF (s/m)')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);
figure(4),hold on           % Empirical CDF of Xl at a typical instant
    i_typ=3;
    stairs(sort(Xsample(:,i_typ)),1/N:1/N:1)
    xlabel('\itx\rm_{10} (m)'),ylabel('CDF')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);
figure(5),hold on           % Empirical CDF of Vl at a typical instant
    i_typ=3;
    stairs(sort(Vsample(:,i_typ)),1/N:1/N:1)
    xlabel('\itv\rm_{10} (m/s)'),ylabel('CDF')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16); 
figure(6),hold on           % Standard deviation history of Xl
    plot(t,sqrt(varX_mcs),'k','LineWidth',0.5)
    xlabel('\itt \rm(s)'),ylabel('StD of \itx\rm_{10} (m)')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);
figure(7),hold on           % Standard deviation history of Vl
    plot(t,sqrt(varV_mcs),'k','LineWidth',0.5)
    xlabel('\itt \rm(s)'),ylabel('StD of \itv\rm_{10} (m/s)')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);
figure(8),hold on           % Kurtosis history of Xl
    plot(t,kurtX_mcs,'k','LineWidth',0.5)
    xlabel('\itt \rm(s)'),ylabel('Kurt of \itx\rm_{10}')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);
figure(9),hold on           % Kurtosis history of Vl
    plot(t,kurtV_mcs,'k','LineWidth',0.5)
    xlabel('\itt \rm(s)'),ylabel('Kurt of \itv\rm_{10}')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);
figure(10),hold on          % Time-variant failure probability under different thresholds
    for i=1:length(threshold)
        plot(t,Pf_mcs(i,:))
        xlabel('\itt \rm(s)'),ylabel('Failure probability')
        box(gca,'on');
        grid(gca,'on');
        set(gca,'FontName','Times New Roman','FontSize',16);
    end