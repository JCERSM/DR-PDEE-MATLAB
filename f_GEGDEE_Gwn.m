function PDFtemp=f_GEGDEE_Gwn(PDF,xl,vl,Aeff,D,dt)
% Function for solution of DR-PDEE via path integral solution (Gaussian white excitaion)
% Input 1: PDF - Transient joint PDF of Xl and Vl at the last step [(Nvl+1) by (Nxl+1) matrix] PDF(i,j)=p[x(j),v(i)]
% Input 2: xl - Displacement grid sequence [(Nxl+1) by 1 vector]
% Input 3: vl - Velocity grid sequence [(Nvl+1) by 1 vector]
% Input 4: Aeff - Effective drift coefficient [(Nvl+1) by (Nxl+1) matrix]
% Input 5: D - Intensity of Gaussian white noise
% Input 6: dt - Time step
% Output: PDFtemp - Transient joint PDF of Xl and Vl at the next step [(Nvl+1) by (Nxl+1) matrix] PDFtemp(i,j)=p[x(j),v(i)]
    Nxl=length(xl)-1;                   % Number of displacement grids
    Nvl=length(vl)-1;                   % Number of velocity grids
    dxl=(xl(end)-xl(1))/Nxl;            % Displacement grid size
    dvl=(vl(end)-vl(1))/Nvl;            % Velocity grid size
    PDF_s=zeros(Nvl+1,Nxl+1);           % Pre-allocated storage for transient joint PDF of Xl and Vl on irregular grid [(Nvl+1) by (Nxl+1) matrix] PDF(i,j)=p[x(j)-dt*v(i),v(i)]
    PDFtemp=zeros(Nvl+1,Nxl+1);         % Pre-allocated storage for transient joint PDF of Xl and Vl at the next step [(Nvl+1) by (Nxl+1) matrix] PDF(i,j)=p[x(j),v(i)]
    Aeff_s=zeros(Nvl+1,Nxl+1);          % Pre-allocated storage for EDC on irregular grid [(Nvl+1) by (Nxl+1) matrix]
    xk=ones(Nvl+1,1)*xl'-vl*ones(1,Nxl+1)*dt;           % Irregular grid [(Nvl+1) by (Nxl+1) matrix] 
    xk_m=xl(end)-vl(1)*dt;                              % Domain of displacement for irregular grid 
    for k=1:Nvl+1
        PDF_s(k,:)=interp1([-xk_m,xl',xk_m],[0,PDF(k,:),0],xk(k,:),'pchip');            % Transient joint PDF of Xl and Vl on irregular grid
        Aeff_s(k,:)=interp1([-xk_m,xl',xk_m],[0,Aeff(k,:),0],xk(k,:));                  % EDC on irregular grid
    end
    p1=max(sum(PDF_s,2)*dxl,1e-10);                                     % Normalization for PDF_s [(Nvl+1) by 1 vector]
    p2=sum(PDF,2)*dxl;                                                  % Normalization for PDF_s [(Nvl+1) by 1 vector]
    PDF_s=PDF_s.*(1+heaviside(abs(p1-p2)-1e-6).*(p2./p1-1));            % Normalization for PDF_s [(Nvl+1) by (Nxl+1) matrix]
    for l=1:Nxl+1
        mu_vl=ones(Nvl+1,1)*(vl'+Aeff_s(:,l)'*dt);              % Expectation of velocity at the next step [(Nvl+1) by (Nvl+1) matrix]
        TPD=normpdf(vl*ones(1,Nvl+1),mu_vl,sqrt(D*dt));         % TPD at x(l) [(Nvl+1) by (Nvl+1) matrix]
        p=ones(Nvl+1,1)*sum(TPD)*dvl;           % Normalization for TPD [(Nvl+1) by (Nvl+1) matrix]
        p(p==0)=1;                              % Normalization for TPD
        TPD=TPD./p;                             % Normalization for TPD [(Nvl+1) by (Nvl+1) matrix]
        PDFtemp(:,l)=TPD*PDF_s(:,l)*dvl;            % Transient joint PDF of Xl and Vl at the next step
    end
    q=sum(sum(PDFtemp))*dxl*dvl;            % Normalization for PDFtemp
    if abs(q-1)>1e-5
        PDFtemp=PDFtemp/q;                  % Normalization for PDFtemp [(Nvl+1) by (Nxl+1) matrix]
    end
end