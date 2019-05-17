%========================================================================
% DESCRIPTION: 
% Investigation of the dynamics of a beam with elastic dry friction
% nonlinearity using NLvib and simulated measurements of the backbone
%========================================================================

clearvars;
close all;
clc;

srcpath = '~/src/matlab/nlvib/SRC';
addpath(genpath(srcpath));
savedata = true;

%% Define system (undamped)

% Properties of the beam
len = 0.7;                % length in m
height = .03;       % height in the bending direction in m
thickness = height;   % thickness in the third dimension in m
E = 185e9;              % Young's modulus in Pa
rho = 7830;             % density in kg/m^3
BCs = 'clamped-free';   % constraints

% Setup one-dimensional finite element model of an Euler-Bernoulli beam
n_nodes = 8;            % number of equidistant nodes along length
beam = FE_EulerBernoulliBeam(len,height,thickness,E,rho,...
    BCs,n_nodes);
n = beam.n;

% Apply elastic Coulomb element at node 4 in translational direction
nl_node = 4; % index includes clamped node
dir = 'trans';
kt = 1.3e6; % tangential stiffness in N/m
muN = 1;    % friction limit force in N
add_nonlinear_attachment(beam,nl_node,dir,'elasticdryfriction',...
    'stiffness',kt,'friction_limit_force',muN);

%% Modal analysis the linearized system

% Modes for free sliding contact
[PHI_free,OM2] = eig(beam.K,beam.M);
om_free = sqrt(diag(OM2));
% Sorting
[om_free,ind] = sort(om_free);
PHI_free = PHI_free(:,ind);

% Modes for fixed contact
K_ex = eye(length(beam.M));
inl = find(beam.nonlinear_elements{1}.force_direction);
K_ex(inl,inl) = kt;
K_ex = beam.K + K_ex;
[PHI_fixed,OM2] = eig(K_ex,beam.M);
om_fixed = sqrt(diag(OM2));
% Sorting
[om_fixed,ind] = sort(om_fixed);
PHI_fixed = PHI_fixed(:,ind);

%% Define linear modal damping

% desired modal damping ratios for first two modes
D1 = 0.008;
D2 = 0.002;

% define Rayleigh damping based on modal damping ratios
beta = 2*(D2*om_fixed(2)-om_fixed(1)*D1)/(om_fixed(2)^2-om_fixed(1)^2);
alpha = 2*om_fixed(1)*D1 - beta*om_fixed(1)^2;

beam.D = alpha*beam.M + beta*K_ex;

% mass-normalized mode shapes
qq =diag(PHI_fixed'*beam.M*PHI_fixed);
PHI_fixed = PHI_fixed./repmat(sqrt(qq)',n,1);

%% Nonlinear modal analysis using harmonic balance
analysis = 'NMA';

% Analysis parameters
H = 9;             % harmonic order
Ntd = 2^10;         % number of time samples per period
imod = 1;           % mode to be analyzed
log10a_s = -5.7;    % start vibration level (log10 of modal mass)
log10a_e = -3;       % end vibration level (log10 of modal mass)
inorm = 2;          % coordinate for phase normalization

% Initial guess vector x0 = [Psi;om;del], where del is the modal
% damping ratio, estimate from underlying linear system
om = om_fixed(imod); phi = PHI_fixed(:,imod);
Psi = zeros((2*H+1)*n,1);
Psi(n+(1:n)) = phi;
x0 = [Psi;om;D(1)];

% Solve and continue w.r.t. Om
ds = .05;
fscl = mean(abs(beam.K*phi));
X_HB = solve_and_continue(x0,...
    @(X) HB_residual(X,beam,H,Ntd,analysis,inorm,fscl),...
    log10a_s,log10a_e,ds);

% Interpret solver output
Psi_HB = X_HB(1:end-3,:);
om_HB = X_HB(end-2,:);
del_HB = X_HB(end-1,:);
log10a_HB = X_HB(end,:);
a_HB = 10.^log10a_HB;
Q_HB = Psi_HB.*repmat(a_HB,size(Psi_HB,1),1);
% fundamental harmonic motion
Y_HB_1 = Q_HB(n+(1:n),:)-1i*Q_HB(2*n+(1:n),:);

if savedata
    save('nma.mat','X_HB','H')
end
