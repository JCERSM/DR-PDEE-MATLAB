%% Stochastic response determination of MDOF Bouc-Wen model enforced by Gaussian white noise via physically-driven DR-PDEE
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
str_column=3;                                               % Number of column of each story
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
T=10;             	% Time domain (s)
dt=0.005;          	% Time step (s)
Neq=800;          	% Number of representative deterministic analyses
a_x=0.5;           	% Domain of Xl (m)
a_v=0.5;          	% Domain of Vl (m/s)
Nxfit=40;          	% Number of grids of Xl for identifying EDC
Nvfit=40;          	% Number of grids of Vl for identifying EDC
Nxl=200;          	% Number of grids of Xl for solving DR-PDEE
Nvl=200;          	% Number of grids of Vl for solving DR-PDEE
alpha=0.2;          % Local smooth scale in (0,1] for identifying EDC
%% Calculation of structural property
bwm=[bwm_alpha,bwm_A,bwm_n,bwm_beta,bwm_gamma,bwm_d_nu,bwm_d_eta,bwm_p,bwm_q,bwm_d_psi,bwm_lambda,bwm_zeta_s,bwm_psi]';         % Parameters of Bouc-Wen model [13 by 1 vector]
str_I=str_b.*str_a.^3/12;                                                           % Moment of inertia of column section of each story (m4) [dof by 1 vector]
str_k=str_column*12*str_E.*str_I./str_h.^3;                                         % Initial lateral inter-story stiffness of each story (N/m) [dof by 1 vector]
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
dxfit=a_x*2/Nxfit;              % Grid size of Xl for identifying EDC (m)
dvfit=a_v*2/Nvfit;              % Grid size of Vl for identifying EDC (m/s)
xfit=-a_x:dxfit:a_x;            % Grid sequence of Xl for identifying EDC (m) [(Nxfit+1) by 1 vector]
vfit=-a_v:dvfit:a_v;            % Grid sequence of Vl for identifying EDC (m/s) [(Nvfit+1) by 1 vector]
dxl=a_x*2/Nxl;                  % Grid size of Xl for solving DR-PDEE (m)
dvl=a_v*2/Nvl;                  % Grid size of Vl for solving DR-PDEE (m/s)
xl=-a_x:dxl:a_x;                % Grid sequence of Xl for solving DR-PDEE (m) [(Nxl+1) by 1 vector]
vl=-a_v:dvl:a_v;                % Grid sequence of Vl for solving DR-PDEE (m/s) [(Nvl+1) by 1 vector]
%% Representative deterministic analyses by 2-order stochastic Runge-Kutta algorithm
Xl=zeros(Neq,Nt);               % Pre-allocated storage for data of Xl [Neq by Nt matrix]
Vl=zeros(Neq,Nt);               % Pre-allocated storage for data of Vl [Neq by Nt matrix]
fl=zeros(Neq,Nt);               % Pre-allocated storage for data of drift force [Neq by Nt matrix]
Xtemp=x0*ones(1,Neq);           % Initial value for all data [4*dof by Neq matrix]
for j=1:Nt
    dW=sqrt(D*dt)*randn(1,Neq);                             % Wiener increments (m/s) [1 by Neq vector]
    F1=f_BoucWen(Xtemp,M,C,str_k,bwm);                      % 2-order stochastic Runge-Kutta algorithm [4*dof by Neq matrix] {call function 'f_BoucWen'}
    F2=f_BoucWen(Xtemp+F1*dt+b*dW,M,C,str_k,bwm);           % 2-order stochastic Runge-Kutta algorithm [4*dof by Neq matrix] {call function 'f_BoucWen'}
    dX=(F1+F2)/2*dt+b*dW;                                   % Responses increments [4*dof by Neq matrix]
    Xtemp=Xtemp+dX;                                         % Responses at the j-th step [4*dof by Neq matrix]
    Xl(:,j)=Xtemp(l_itrst,:)';                                                	% Data of Xl (m)
    Vl(:,j)=Xtemp(dof+l_itrst,:)';                                             	% Data of Vl (m/s)
    DX=Xtemp(1:dof,:)-[zeros(1,Neq);Xtemp(1:dof-1,:)];                         	% Inter-story drift of each story at the j-th step (m) [dof by Neq vector]
    DG=str_k.*(bwm_alpha*DX+(1-bwm_alpha)*Xtemp(2*dof+1:3*dof,:));              % Inter-story restoring force of each story at the j-th step (N) [dof by Neq vector]
    f_vct=-invM*(C*Xtemp(dof+1:2*dof,:)+DG-[DG(2:end,:);zeros(1,Neq)]);         % Drift force of each story at the j-th step (m/s2) [dof by Neq vector]
    fl(:,j)=f_vct(l_itrst,:)';                                               	% Data of drift force (m/s2)
end
%% Identification of effective drift coefficients
Aeff=zeros(Nvfit+1,Nxfit+1,Nt);         % Pre-allocated storage for effective drift coefficient [(Nvfit+1) by (Nxfit+1) by Nt array]
parfor j=1:Nt
    Xdata=[Xl(:,j);-Xl(:,j)];           % Used data of Xl (m)
    Vdata=[Vl(:,j);-Vl(:,j)];        	% Used data of Vl (m/s)
    Fdata=[fl(:,j);-fl(:,j)];          	% Used data of drift force (m/s2)
    Aeff(:,:,j)=f_EDC_lowess(xfit',vfit',Xdata,Vdata,Fdata,alpha);          % Effective drift coefficient (m/s2) {call function 'f_EDC_lowess'}
    fprintf('Identification of EDC at the %i-th step is in progress\n',j)
end
%% Solution of DR-PDEE
PDF_xl=zeros(Nxl+1,Nt+1);                   % Pre-allocated storage for PDF of Xl [(Nxl+1) by (Nt+1) matrix]
PDF_vl=zeros(Nvl+1,Nt+1);                   % Pre-allocated storage for PDF of Vl [(Nvl+1) by (Nt+1) matrix]
PDF=zeros(Nvl+1,Nxl+1);                     % Initial value of joint PDF of Xl and Vl (s/m2) [(Nvl+1) by (Nxl+1) matrix] PDF(i,j)=p[x(j),v(i)]
PDF(Nvl/2+1,Nxl/2+1)=1/dxl/dvl;             % Initial value of joint PDF of Xl and Vl (s/m2)
[XXfit,VVfit]=meshgrid(xfit,vfit);          % Mesh division for identifying EDC
[XXl,VVl]=meshgrid(xl,vl);                  % Mesh division for solving DR-PDEE
for j=1:Nt
    Aeff_used=interp2(XXfit,VVfit,Aeff(:,:,j),XXl,VVl);         % EDC on refined grid (m/s2) [(Nvl+1) by (Nxl+1) matrix]
    PDF=f_GEGDEE_Gwn(PDF,xl',vl',Aeff_used,D,dt);               % Transient joint PDF of Xl and Vl (s/m2) [(Nvl+1) by (Nxl+1) matrix]
    PDF_xl(:,j+1)=sum(PDF,1)'*dvl;                              % PDF of Xl (/m)
    PDF_vl(:,j+1)=sum(PDF,2)*dxl;                               % PDF of Vl (s/m)
    fprintf('Solution of DR-PDEE at the %i-th step is in progress\n',j)
end
varX=sum(xl'.^2.*PDF_xl).*dxl;                  % Variance of Xl (m2)
varV=sum(vl'.^2.*PDF_vl).*dvl;                  % Variance of Vl (m2/s2)
kurtX=sum(xl'.^4.*PDF_xl)*dxl./varX.^2;         % Kurtosis of Xl
kurtV=sum(vl'.^4.*PDF_vl)*dvl./varV.^2;         % Kurtosis of Vl
save('Resp_BW_Gwn.mat','PDF_xl','PDF_vl','xl','vl','varX','varV','kurtX','kurtV','Aeff','XXfit','VVfit','Xl','Vl','fl','t','dt','dxl','dvl');
%% Comprison with MCS
load('Resp_BW_Gwn.mat');
load('MCS_BW_Gwn.mat');
%% Figures
figure(1),hold on           % EDC at a typical instant
    T0=10;
    scatter3(Xl(:,T0/dt),Vl(:,T0/dt),fl(:,T0/dt))
    mesh(XXfit,VVfit,Aeff(:,:,T0/dt))
    xlabel('\itx\rm_{10} (m)'),ylabel('\itv\rm_{10} (m/s)'),zlabel('\ita\rm^{(eff)} (m/s^2)')
    legend('Original data','Effective drift coefficient')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);
figure(2),hold on           % PDF of Xl at a typical instant
    i_typ=3;
    T0=Ttyp(i_typ);
    plot(xl,PDF_xl(:,T0/dt+1)','r','LineWidth',1.2)
    histogX=reshape([x_bar(1:end-1,i_typ)';x_bar(1:end-1,i_typ)';x_bar(2:end,i_typ)';x_bar(2:end,i_typ)'],4*N_bar,1);
    histogPx=reshape([ones(1,N_bar)*1e-10;PDF_xl_mcs(:,i_typ)';PDF_xl_mcs(:,i_typ)';ones(1,N_bar)*1e-10],4*N_bar,1);
    plot(histogX,histogPx,'b','LineWidth',0.5)
    xlabel('\itx\rm_{10} (m)'),ylabel('PDF (m^{-1})')
    legend('Physically-driven DR-PDEE','MCS with 10^6 samples')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);
figure(3),hold on           % PDF of Vl at a typical instant
    i_typ=3;
    T0=Ttyp(i_typ);
    plot(vl,PDF_vl(:,T0/dt+1)','r','LineWidth',1.2)
    histogV=reshape([v_bar(1:end-1,i_typ)';v_bar(1:end-1,i_typ)';v_bar(2:end,i_typ)';v_bar(2:end,i_typ)'],4*N_bar,1);
    histogPv=reshape([ones(1,N_bar)*1e-10;PDF_vl_mcs(:,i_typ)';PDF_vl_mcs(:,i_typ)';ones(1,N_bar)*1e-10],4*N_bar,1);
    plot(histogV,histogPv,'b','LineWidth',0.5)
    xlabel('\itv\rm_{10} (m/s)'),ylabel('PDF (s/m)')
    legend('Physically-driven DR-PDEE','MCS with 10^6 samples')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);
figure(4),hold on           % CDF of Xl at a typical instant
    i_typ=3;
    T0=Ttyp(i_typ);
    plot(xl,cumsum(PDF_xl(:,T0/dt+1))'*dxl,'r','LineWidth',0.5)
    stairs(sort(Xsample(:,i_typ)),1/N:1/N:1,'b--','LineWidth',1.2)
    xlabel('\itx\rm_{10} (m)'),ylabel('CDF')
    legend('Physically-driven DR-PDEE','MCS with 10^6 samples')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);
figure(5),hold on           % CDF of Vl at a typical instant
    i_typ=3;
    T0=Ttyp(i_typ);
    plot(vl,cumsum(PDF_vl(:,T0/dt+1))'*dvl,'r','LineWidth',0.5)
    stairs(sort(Vsample(:,i_typ)),1/N:1/N:1,'b--','LineWidth',1.2)
    xlabel('\itv\rm_{10} (m/s)'),ylabel('CDF')
    legend('Physically-driven DR-PDEE','MCS with 10^6 samples')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);
figure(6),hold on           % Standard deviation history of Xl
    plot(t,sqrt(varX),'r','LineWidth',0.5)
    plot(t,sqrt(varX_mcs),'b--','LineWidth',1.2)
    xlabel('\itt\rm (s)'),ylabel('StD of \itx\rm_{10} (m)')
    legend('Physically-driven DR-PDEE','MCS with 10^6 samples')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);
figure(7),hold on           % Standard deviation history of Vl
    plot(t,sqrt(varV),'r','LineWidth',0.5)
    plot(t,sqrt(varV_mcs),'b--','LineWidth',1.2)
    xlabel('\itt\rm (s)'),ylabel('StD of \itv\rm_{10} (m/s)')
    legend('Physically-driven DR-PDEE','MCS with 10^6 samples')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);
figure(8),hold on           % Kurtosis history of Xl
    plot(t,kurtX,'r','LineWidth',0.5)
    plot(t,kurtX_mcs,'b--','LineWidth',1.2)
    xlabel('\itt\rm (s)'),ylabel('Kurt of \itx\rm_{10}')
    legend('Physically-driven DR-PDEE','MCS with 10^6 samples')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);
figure(9),hold on           % Kurtosis history of Vl
    plot(t,kurtV,'r','LineWidth',0.5)
    plot(t,kurtV_mcs,'b--','LineWidth',1.2)
    xlabel('\itt\rm (s)'),ylabel('Kurt of \itv\rm_{10}')
    legend('Physically-driven DR-PDEE','MCS with 10^6 samples')
    box(gca,'on');
    grid(gca,'on');
    set(gca,'FontName','Times New Roman','FontSize',16);