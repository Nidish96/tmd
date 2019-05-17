%========================================================================
% DESCRIPTION: 
% Investigation of the dynamics of an Euler-Bernoulli beam
% with regularized Coulomb dry friction nonlinearity.
%========================================================================

clearvars;
close all;
clc;

srcpath = '~/src/matlab/nlvib/SRC';
addpath(genpath(srcpath));
%% Define system

% Properties of the beam
len = 0.7;              % length in m
height = .03;           % height in the bending direction in m
thickness = height;     % thickness in the third dimension in m
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

% Vectors recovering deflection at tip and nonlinearity location
T_nl = beam.nonlinear_elements{1}.force_direction';
T_tip = beam.L(end-1,:);

dKfric = kt*(T_nl*T_nl');

%% Modal analysis the linearized system

% Modes for free sliding contact
[PHI_free,OM2] = eig(beam.K,beam.M);
om_free = sqrt(diag(OM2));
% Sorting
[om_free,ind] = sort(om_free);
PHI_free = PHI_free(:,ind);

% Modes for fixed contact
K_ex = eye(length(beam.M));
inl = find(T_nl);
K_ex(inl,inl) = kt;
K_ex = beam.K + K_ex;
[PHI_fixed,OM2] = eig(K_ex,beam.M);
om_fixed = sqrt(diag(OM2));
% Sorting
[om_fixed,ind] = sort(om_fixed);
PHI_fixed = PHI_fixed(:,ind);

% constraining matrix - alternative way to fix dofs
B = eye(length(beam.M)); B(:,inl) = [];
% [PHI_fixed,OM2] = eig(B'*beam.K*B,B'*beam.M*B);
% om_fixed = sqrt(diag(OM2));
% % Sorting
% [om_fixed,ind] = sort(om_fixed); PHI_fixed = B*PHI_fixed(:,ind);

%% Define linear modal damping
% desired modal damping ratios for first two modes
D1 = 0.008;
D2 = 0.002;

% define Rayleigh damping based on modal damping ratios
beta = 2*(D2*om_fixed(2)-om_fixed(1)*D1)/(om_fixed(2)^2-om_fixed(1)^2);
alpha = 2*om_fixed(1)*D1 - beta*om_fixed(1)^2;

beam.D = alpha*beam.M + beta*K_ex;
%% Nonlinear frequency response analysis using Harmonic Balance
analysis = 'FRF';

% Analysis parameters
H = 13;                      % harmonic order
N = 2^10;                    % number of time samples per period
% freq interval. We go from high to low, as the peak is bent towards left
% Om_s = 1.5*om_fixed(1);     % start frequency
% Om_e = .7*om_free(1);       % end frequency
Om_e = 500;
Om_s = 200;

% Apply forcing to free end of beam in translational direction, with
% magnitude 'fex'
inode = 8; % node for external excitation. equal to n_nodes.
dir = 'trans';

% create force vector and find index. The actual value is overwritten later
fex = 1;
add_forcing(beam,inode,dir,fex);
idx_f = find(beam.Fex1);

% Excitation levels
exc_lev = fliplr([0.01, 0.05, 0.1, 0.2, 0.3, 0.6]);
Om = cell(size(exc_lev));
for iex=1:length(exc_lev)
    % update force. Note that calling add_force is additive -> here we
    % overwrite
    fex = exc_lev(iex);
    beam.Fex1(idx_f) = fex;

    % Initial guess (from underlying linear system)
    Fex1 = beam.Fex1;
    % Q1 = B*((B'*(-Om_s^2*beam.M + 1i*Om_s*beam.D + K_ex)*B)\(B'*Fex1));
    Q1 = (-Om_s^2*beam.M + 1i*Om_s*beam.D + beam.K + dKfric)\Fex1;
    
    x0 = zeros((2*H+1)*size(Q1,1),1);
    x0(size(Q1,1)+(1:2*size(Q1,1))) = [real(Q1);-imag(Q1)];

    % Solve and continue w.r.t. Om
    ds = 1;
    % Set options of path continuation
    % flag = 0: no actual continuation, just sequential continuation
    % stepadatp = 0: no automatic step size adjustment
    
    % introduce a 'scaling' depending on the excitation level. For this,
    % the linear resonant response level, assuming sticking conditions, and
    % applying a truncation to a single mode, is determined analytically.
    % The scaling vector for the displacement unknowns is then defined
    % accordingly. In my experience, this makes a lot of sense if you want
    % to cover a large range of excitation levels.
    aref = abs((PHI_fixed(:,1)'*Fex1)/(2*D1*om_fixed(1)^2));
    Dscale = [1e1*aref*ones(length(x0),1); (Om_s+Om_e)/2];
    Sopt = struct('Dscale',Dscale,'dynamicDscale',1);
    X_HB = solve_and_continue(x0,...
        @(X) HB_residual(X,beam,H,N,analysis),...
        Om_s,Om_e,ds,Sopt);

    % Interpret solver output
    OM_HB{iex} = X_HB(end,:);
    Q_HB{iex} = X_HB(1:end-1,:);

    % Determine harmonics and root-mean-square value of tip displacement
    Qtip_HB = kron(eye(2*H+1),T_tip)*Q_HB{iex};
    Qtip_rms_HB{iex} = sqrt(sum(Qtip_HB.^2))/sqrt(2);
end

save('hb.mat','OM_HB','Q_HB','Qtip_HB','Qtip_rms_HB','beam','exc_lev',...
    'H','N','ds','n_nodes')
