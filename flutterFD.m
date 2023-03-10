function [Ucr,wCr,U] = flutterFD(Bridge,varargin)
% [Ucr,fCr] = flutterFD(Bridge,varargin) computes the flutter velocity of a
% single-span suspension bridge using a frequency-domain and a multi-mode
% approach. The method used is from [1].
% 
% Inputs:
% Bridge: structure with 13 fields:
% Bridge.wn: [2 x Nmodes] matrix of vertical and torsional eigen frequencies  (rad/s)
% Bridge.phi : [2 x Nyy x Nmodes] matrix of vertical and torsional mode shapes
% Bridge.zetaStruct: [2 x Nmodes] matrix of structural modal damping ratios
% Bridge.L: [1 x 1] scalar: length of main span (m) 
% Bridge.B: [1 x 1] scalar: deck width (m)
% Bridge.D: [1 x 1] scalar: Deck height (m)
% Bridge.m: [1 x 1] scalar: lineic mass of the girder + the two main cables (kg/m) 
% Bridge.m_theta: [1 x 1] scalar:  torsional mass (kg.m^2/m )
% Bridge.Cd: [1 x 1] scalar:  drag coefficient
% Bridge.dCd: [1 x 1] scalar: first derivative of the drag coefficient
% Bridge.Cl: [1 x 1] scalar:  lift coefficient
% Bridge.dCl: [1 x 1] scalar: first derivative of the lift coefficient
% Bridge.Cm: [1 x 1] scalar: Overturning moment coefficient
% Bridge.dCm: [1 x 1] scalar: first derivative of the overturning moment coefficient
% 
% Outputs:
% Ucr: [1 x 1] scalar: Estimate of critical wind speed (m/s)
% wcr: [M x 2 x Nmodes] matrix of vertical and torsional eigen frequencies
% for increasing wind velocities (rad/s), where M is the number of
% iteration until the flutter velocity is reached.
% U : [1 x M] vector of mean wind velocities tested (m/s).
%% Author Info
% E. Cheynet - University of Stavanger - Norway - last modified 22.04.2018
% 
%% References 
% [1] Jain, A., Jones, N. P., & Scanlan, R. H. (1996). 
% Coupled aeroelastic and aerodynamic response analysis of long-span bridges.
% Journal of Wind Engineering and Industrial Aerodynamics, 60, 69-80.

%% Inputparseer (optional parameters)
p = inputParser();
p.CaseSensitive = false;
p.addOptional('Nfreq',2000); % Number of frequency step
p.addOptional('Niter',100); % Number of iterations for the mean wind velocity
p.addOptional('m12',0); % structural coupling between lateral and vertical motion (percentage of Bridge.m)
p.addOptional('m13',0); % structural coupling between lateral and torsional motion  (percentage of Bridge.m)
p.addOptional('m23',0); % structural coupling between vertical and torsional motion  (percentage of Bridge.m)
p.addOptional('kTheta',1/4); % bridge constant introducing an additional coupling between the vertical and torsional motion. cf aerodynamic of streamlined bridge
p.addOptional('rho',1.25); % air density
p.addOptional('fmin',[]); % Minimal frequency for the computation of the critical velocity in the frequency domain
p.addOptional('fmax',[]); % Maximal frequency for the computation of the critical velocity in the frequency domain
p.addOptional('Umin',5); % Lowest boundary of the velocity range investigated
p.addOptional('Umax',200); % upper boundary of the velocity range investigated
p.parse(varargin{:});
% shorthen the variables name
Niter = p.Results.Niter ;
Nfreq = p.Results.Nfreq ;
m12 = p.Results.m12;
m13 = p.Results.m13;
m23= p.Results.m23;
kTheta = p.Results.kTheta;
rho = p.Results.rho;
fmin = p.Results.fmin;
fmax = p.Results.fmax;
Umin = p.Results.Umin;
Umax = p.Results.Umax;
if abs(m12)>=100 ||  abs(m13)>=100 || abs(m13) >=100,
    warning(' m12, m13 or m23 are larger than 1. Error may arise');
end
%% The variable names are shortened
Cd = Bridge.Cd;
dCd = Bridge.dCd;

Cl = Bridge.Cl;
dCl = Bridge.dCl;

Cm =Bridge.Cm;
dCm =Bridge.dCm;

if ~isscalar(Cd)||~isscalar(Cl)||~isscalar(Cm)||~isscalar(dCd)||~isscalar(dCl)||~isscalar(dCm),
    error('Each of the aerodynamic coefficient should be scalar.')
end

D = Bridge.D;
B = Bridge.B;
[Ndof,Nmodes,Nyy] = size(Bridge.phi);
x = linspace(0,Bridge.L,Nyy);
wn = Bridge.wn(:);
if size(wn,1)~= Ndof*Nmodes,    error('The size of wn and phi are not consistent');end
zeta = Bridge.zetaStruct(:);
Nyy = numel(x);
phi0 = Bridge.phi;
phi = reshape(phi0,[],Nyy);
newN = size(phi,1);

if isempty(fmin),fmin = 0.1*min(wn)/2/pi;end
if isempty(fmax),fmax = 10*min(wn)/2/pi;end
f = logspace(log10(fmin),log10(fmax),Nfreq);
%% Modal parameter calculation
M = diag([repmat(Bridge.m,1,Ndof-1),Bridge.m_theta]');
M(Ndof-1,Ndof)=m12*M(Ndof-1,Ndof-1);
M(Ndof,Ndof-1)=m12*M(Ndof-1,Ndof-1);

if Ndof==3
M(1,3)=m13*M(1,1);
M(3,1)=m13*M(1,1);
M(2,3)=m23*M(2,2);
M(3,2)=m23*M(2,2);
end

U = linspace(Umin,Umax,Niter);
wCr = nan(Niter,Ndof,newN/Ndof);
for oo=1:Niter,
    
    [Kae,Cae] = AeroPara(U(oo));
    [Cae_modal,Kae_modal,Mtot,K_modal,C_modal] = getModalPara(x,phi,Cae,Kae,M,wn,zeta);
    Ktot = K_modal-Kae_modal;
    Ctot = C_modal-Cae_modal;
    
    dummyF = zeros(1,newN);
    for jj=1:newN,
        [dummyF(1,jj),~] = eigenFreqSyst(Mtot(jj,jj),Ktot(jj,jj),Ctot(jj,jj));
    end
    wCr(oo,:,:) = reshape(2*pi*dummyF',Ndof,[]);
    
    
    % Test Divergence
    if any(diag(Ktot)<0)
        warning('static torsional divergence reached ! (or maybe Nfreq and Niter are not properly defined)')
        Ucr = U(oo);
        wCr(U>Ucr,:,:)= [];
        U(U>Ucr)= [];
        fprintf(['Vcr = ',num2str(Ucr,5),' m/s \n']);
        return
    end
    
    r = zeros(1,Nfreq);
    for ii=1:Nfreq,
        E = mecaAdmittance(1i*2*pi*f(ii),Mtot,Ctot,Ktot);
        r(ii) = det(E);
    end
    

       
    Re = zerocross(real(r));
    Im = zerocross(imag(r));
    C = intersect(Re,Im);
    
    if oo==Niter && isempty(C),
        warning('No critical wind velocity found')
        Ucr = nan;
        wCr = nan(newN,1);
        return
    elseif ~isempty(C)
        Ucr = U(oo);
        fprintf('Critical velocity for coupled flutter found! \n');
        fprintf(['Vcr = ',num2str(Ucr,5),' m/s \n']);
        wCr(U>Ucr,:,:)= [];
        U(U>Ucr)= [];
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%% Nested functions
    function z=zerocross(v)
        z=find(diff(v>0)~=0)+1;
    end
    function [Kae,Cae] = AeroPara(U)
        CST = 0.5*B*U*rho;
        
        if Ndof ==2,
            Cae  =-CST.*...
                [dCl+D/B*Cd,kTheta*B*(dCl+D/B*Cd) ;...
                B*dCm, kTheta*dCm*B*B];% 3 x 3
            Kae  = CST.*U.*...
                [0 ,dCl ;...
                0, dCm*B];
        else
            Cae  =-CST.*...
                [2*D/B*Cd, D/B*dCd-Cl,kTheta*B*(D/B*dCd-Cl);...
                2*Cl,dCl+D/B*Cd,kTheta*B*(dCl+D/B*Cd) ;...
                2*B*Cm, B*dCm, kTheta*dCm*B*B];% 3 x 3
            Kae  = CST.*U.*...
                [0, 0, D/B*dCd;...
                0,0 ,dCl ;...
                0, 0, dCm*B];
        end
    end
    function [Cae_modal,Kae_modal,M_modal,K_modal,C_modal] = getModalPara(x,phi,Cae,Kae,M,wn,zeta)
        Cae_modal = zeros(newN); %modal aerodynamic mass
        Kae_modal = zeros(newN); %modal aerodynamic mass
        M_modal = zeros(newN);
        M = repmat(M,round(newN/Ndof),round(newN/Ndof));    
        for pp=1:newN,
            for qq=1:newN/Ndof,
                if (pp+Ndof*qq)<=newN,
                    M(pp,pp+Ndof*qq) = 0;
                    M(pp+Ndof*qq,pp) = 0;
                end
            end
        end
        
        Cae = repmat(Cae,round(newN/Ndof),round(newN/Ndof));
        Kae = repmat(Kae,round(newN/Ndof),round(newN/Ndof));
        for pp=1:newN,
            for qq=1:newN,
                M_modal(pp,qq)  = trapz(x,phi(pp,:).*phi(qq,:).*M(pp,qq));
                Cae_modal(pp,qq) = trapz(x,phi(pp,:).*phi(qq,:).*Cae(pp,qq));
                Kae_modal(pp,qq) = trapz(x,phi(pp,:).*phi(qq,:).*Kae(pp,qq));
            end
        end
        K_modal = diag(wn).^2*M_modal;
        C_modal = 2.*diag(wn)*M_modal*diag(zeta);
        
    end
    function [E] = mecaAdmittance(w,M,C,K)
        E=K+ C.*w+M.*w.^2;
    end
    function [fC,zeta] = eigenFreqSyst(M,K,C)
        
        % INPUT
        % M is a [newN x newN] matrix
        % K is a [newN x newN] matrix
        % C is a [newN x newN] matrix
        
        % OUTPUT
        % fn is [1 x newN]
        % zeta is [ newN x newN]
        
        A = [zeros(unique(size(M))),eye(unique(size(M)));...
            -inv(M)*K, -inv(M)*C];
        
        [~,Di] = eig(A);
        ReD=real(Di);
        ImD = imag(Di);
        fC= diag(abs(Di(1:2:end,1:2:end))/(2*pi));
        zeta = -ReD(1:2:end,1:2:end)/abs(Di(1:2:end,1:2:end));
        
    end
end