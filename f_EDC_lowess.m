function Aeff=f_EDC_lowess(xl,vl,X,V,F,alpha)
% Function for identification of effective drift coefficients via locally weighted smooth scatterplots
% Input 1: xl - Displacement grid sequence [(Nxl+1) by 1 vector]
% Input 2: vl - Velocity grid sequence [(Nvl+1) by 1 vector]
% Input 3: X - Data of displacement [N by 1 vector]
% Input 4: V - Data of velocity [N by 1 vector]
% Input 5: F - Data of drift force [N by 1 vector]
% Input 6: alpha - Local smooth scale in (0,1]
% Output: Aeff - Effective drift coefficient [(Nvl+1) by (Nxl+1) matrix]
    N=length(F);                            % Number of data samples
    Nxl=length(xl)-1;                       % Number of displacement grids
    Nvl=length(vl)-1;                       % Number of velocity grids
    n=(Nxl+1)*(Nvl+1);                      % Number of grid points
    xp=reshape(ones(Nvl+1,1)*xl',1,n);      % Displacement grid points [1 by n vector]
    vp=reshape(vl*ones(1,Nxl+1),1,n);       % Velocity grid points [1 by n vector]
    r=round(alpha*N);                       % Number of local used data samples
    sigmaX=max(std(X),1e-10);               % Standard deviation of displacement samples
    sigmaV=max(std(V),1e-10);               % Standard deviation of velocity samples
    Mdlkdt=KDTreeSearcher([X./sigmaX,V./sigmaV]);                   % Build a KD-Tree structure
    [ind,d]=knnsearch(Mdlkdt,[xp',vp']./[sigmaX,sigmaV],'K',r);     % KNN query: ind - Index of locally selected points [n by r matrix]; d - distance to locally selected points [n by r matrix]
    Xknn=X(ind');                                   % Data of displacement at locally selected points [r by n matrix]
    Vknn=V(ind');                                   % Data of velocity at locally selected points [r by n matrix]
    Fknn=F(ind');                                   % Data of drift force at locally selected points [r by n matrix]
    h=d(:,end)';                                    % Bandwidth round each grid point [1 by n vector]
    w=heaviside(1-d'./h).*(1-(d'./h).^3).^3;        % Weight [r by n matrix]
    w0=sum(w);                                      % Coefficient used for locally weighted least square [1 by n vector]
    wX=sum(w.*Xknn);                                % Coefficient used for locally weighted least square [1 by n vector]
    wV=sum(w.*Vknn);                                % Coefficient used for locally weighted least square [1 by n vector]
    wX2=sum(w.*Xknn.^2);                            % Coefficient used for locally weighted least square [1 by n vector]
    wXV=sum(w.*Xknn.*Vknn);                         % Coefficient used for locally weighted least square [1 by n vector]
    wV2=sum(w.*Vknn.^2);                            % Coefficient used for locally weighted least square [1 by n vector]
    wF=sum(w.*Fknn);                                % Coefficient used for locally weighted least square [1 by n vector]
    wXF=sum(w.*Xknn.*Fknn);                         % Coefficient used for locally weighted least square [1 by n vector]
    wVF=sum(w.*Vknn.*Fknn);                         % Coefficient used for locally weighted least square [1 by n vector]
    D0=w0.*wX2.*wV2+2*wX.*wV.*wXV-w0.*wXV.^2-wX.^2.*wV2-wV.^2.*wX2;                     % Determinant operation for locally weighted least square [1 by n vector]
    D1=wX2.*wV2.*wF+wV.*wXV.*wXF+wX.*wXV.*wVF-wXV.^2.*wF-wX.*wV2.*wXF-wV.*wX2.*wVF;     % Determinant operation for locally weighted least square [1 by n vector]
    D2=wV.*wXV.*wF+w0.*wV2.*wXF+wX.*wV.*wVF-wX.*wV2.*wF-wV.^2.*wXF-w0.*wXV.*wVF;        % Determinant operation for locally weighted least square [1 by n vector]
    D3=wX.*wXV.*wF+wX.*wV.*wXF+w0.*wX2.*wVF-wV.*wX2.*wF-w0.*wXV.*wXF-wX.^2.*wVF;        % Determinant operation for locally weighted least square [1 by n vector]
    beta=[D1;D2;D3]./max(D0,1e-10);                     % Local regression coefficient [3 by n matrix]
    Aeff_p=beta(1,:)+beta(2,:).*xp+beta(3,:).*vp;       % Effective drift coefficient at each grid point [1 by n vector]
    Aeff=reshape(Aeff_p,Nvl+1,Nxl+1);                   % Effective drift coefficient [(Nvl+1) by (Nxl+1) matrix]
end