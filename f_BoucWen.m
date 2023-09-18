function y=f_BoucWen(x,M,C,k0,bwm)
% Function for solution of MDOF equation of motion with Bouc-Wen restoring force model
% Input 1: x - Responses (displacements, velocities, hysteretic displacements, and hysteretic-dissipated energies) at the current instant [4*dof by N vector]
% Input 2: M - Lumped mass (kg) [dof by dof matrix]
% Input 3: C - Rayleigh damping (kg/s) [dof by dof matrix]
% Input 4: k0 - Initial lateral inter-story stiffness of each story (N/m) [dof by 1 vector]
% Input 5: bwm - Parameters of Bouc-Wen model [13 by 1 vector]
% Output: y - Derivative of responses at the current instant [4*dof by N vector]
    [~,N]=size(x);              % Number of samples
    dof=length(k0);             % Degrees of freedom
    y=zeros(4*dof,N);           % Pre-allocated storage for derivative of responses at the current instant
    alpha=bwm(1);           % Ratio of linear to nonlinear response
    A=bwm(2);               % Basic hysteresis shape controls of Bouc-Wen model
    n=bwm(3);               % Sharpness of yield of Bouc-Wen model
    beta=bwm(4);            % Basic hysteresis shape controls of Bouc-Wen model (/m)
    gamma=bwm(5);           % Basic hysteresis shape controls of Bouc-Wen model (/m)
    d_nu=bwm(6);            % Strength degradation of Bouc-Wen model (/m2)
    d_eta=bwm(7);           % Stiffness degradation of Bouc-Wen model (/m2)
    p=bwm(8);               % Pinching slope of Bouc-Wen model (/m2)
    q=bwm(9);               % Pinching initiation of Bouc-Wen model
    d_psi=bwm(10);          % Pinching rate of Bouc-Wen model (/m2)
    lambda=bwm(11);         % Pinching severity/rate interaction of Bouc-Wen model
    zeta_s=bwm(12);         % Measure of total slip of Bouc-Wen model
    psi=bwm(13);            % Pinching magnitude of Bouc-Wen model (m)
    X=x(1:dof,:);                        	% Displacement of each story (m) [dof by N vector]
    V=x(dof+1:2*dof,:);                  	% Velocity of each story (m/s) [dof by N vector]
    Z=x(2*dof+1:3*dof,:);                	% Imaginary hysteretic displacement of each story (m) [dof by N vector]
    epsilon=x(3*dof+1:4*dof,:);            	% Hysteretic-dissipated energiy of each story (m2) [dof by N vector]
    DX=X-[zeros(1,N);X(1:dof-1,:)];        	% Inter-story drift of each story (m) [dof by N vector]
    DV=V-[zeros(1,N);V(1:dof-1,:)];         % Inter-story velocity of each story (m/s) [dof by N vector]
    DG=k0.*(alpha*DX+(1-alpha)*Z);          % Inter-story restoring force of each story (N) [dof by N vector]
    G=DG-[DG(2:dof,:);zeros(1,N)];        	% Restoring force of each story (N) [dof by N vector]
    nu=1+d_nu*epsilon;                                       	% Strength degradation shape of each story [dof by N vector]
    eta=1+d_eta*epsilon;                                       	% Stiffness degradation shape of each story [dof by N vector]
    zeta1=zeta_s*(1-exp(-p*epsilon));                           % Controling progress od pinching of each story [dof by N vector]
    zeta2=(psi+d_psi*epsilon).*(lambda+zeta1);                  % Controling progress od pinching of each story (m) [dof by N vector]
    Zu=(A./nu/(beta+gamma)).^(1/n);                             % Ultimate value of hysteretic displacement of each story (m) [dof by N vector]
    h=1-zeta1.*exp(-(Z.*sign(DV)-q*Zu).^2./zeta2.^2);           % Pinching shape of each story [dof by N vector]
    y(1:dof,:)=V;                                                                                       % Derivative of displacement of each story (m/s) [dof by N vector]
    y(dof+1:2*dof,:)=-M\(C*V+G);                                                                        % Derivative of velocity of each story (m/s2) [dof by N vector]
    y(2*dof+1:3*dof,:)=h.*(A*DV-nu.*(beta*abs(DV).*abs(Z).^(n-1).*Z+gamma*DV.*abs(Z).^n))./eta;         % Derivative of imaginary hysteretic displacement of each story (m/s) [dof by N vector]
    y(3*dof+1:4*dof,:)=Z.*DV;                                                                           % Derivative of hysteretic-dissipated energiy of each story (m2/s) [dof by N vector]
end