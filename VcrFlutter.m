function [Vcr] = VcrFlutter(B,m,m_theta,fz,ftheta, varargin)
% 
% [Vcr] = VcrFlutter(B,Mg,Jo,fz,ftheta,method) computes the critical
% flutter velocity for a single-span suspension bridge. The case of coupled
% flutter vertical-torsion is considered here for a single-mode approach.
%
% Input
% B: [1 x 1] scalar:  Bridge deck width (m)
% m: [1 x 1] scalar: lineic mass of the girder + the two main cables (kg/m) 
% m_theta: [1 x 1] scalar: torsional mass of the girder (kg.m^2/m)
% fz: [1 x 1] scalar:  Eigen frequency of the vertical eigenmode selected
% ftheta: [1 x 1] scalar:  Eigen frequency of the torsional eigenmode selected
% Optional: 
%           'method': string: 'Selberg1', 'Selberg2' or 'Rocard'
%           'rho':[1 x 1] scalar: air density (default value of 1.25
%           kg/m^3)
% Output
% Vcr: [1 x 1] scalar: Estimate of critical wind speed (m/s)
%
%% Author Info
% Etienne Cheynet - last modified 22.04.2018

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('method','Selberg1');
p.addOptional('rho',1.25);
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
method = p.Results.method ; %
rho = p.Results.rho;
gamma = fz./ftheta;
DD = (1-gamma.^2);

switch method
    case 'Selberg1'
        Rg = sqrt (m_theta/m);
        A = m.*Rg./(rho*B.^3).*DD;
        Vcr = 3.72*B*ftheta.*sqrt(A);
    case 'Selberg2'
        % E. Selberg, E. Hjorth-Hansen, Aerodynamic stability and related aspects of suspension bridges, Paper no. 20,1961
        nu = 8*m_theta./(B.^2*m);
        mu = pi*rho*B.^2./2/m;
        Vcr = 0.44*2*pi.*ftheta*B.*sqrt(DD.*sqrt(nu)./mu);
    case 'Rocard'
        % Y. Rocard, Instabilite des ponts suspendus dans le vent勇xperiences sur mode le reduit, Nat. Phys. Lab.Paper 10 (1963) Teddington, England.
        nu = 8*m_theta./(B.^2*m);
        mu = pi*rho*B.^2./2/m;
        Vcr =  0.443*2*pi.*ftheta*B.*sqrt(DD./mu.*(2*nu)./(1+nu));
end


end

