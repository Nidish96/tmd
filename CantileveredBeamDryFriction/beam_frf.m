%========================================================================
% DESCRIPTION: 
% Investigation of the dynamics of an Euler-Bernoulli beam
% with regularized Coulomb dry friction nonlinearity.
%
%========================================================================

clearvars;
close all;
clc;

srcpath = '~/src/matlab/nlvib/SRC';
addpath(genpath(srcpath));
%% Define system

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

% mass-normalized mode shapes
qq =diag(PHI_fixed'*beam.M*PHI_fixed);
PHI_fixed = PHI_fixed./repmat(sqrt(qq)',n,1);

cc = diag(PHI_fixed'*beam.D*PHI_fixed);
D = cc./(2*om_fixed); % modal damping ratios

%% Nonlinear frequency response analysis using Harmonic Balance
analysis = 'FRF';

% Analysis parameters
H = 9;                      % harmonic order
N = 2^9;                    % number of time samples per period
% freq interval. We go from high to low, as the peak is bent towards left
Om_s = 1.5*om_fixed(1);     % start frequency
Om_e = .9*om_free(1);       % end frequency

% Apply forcing to free end of beam in translational direction, with
% magnitude 'fex'
inode = 8; % node for external excitation. equal to n_nodes.
dir = 'trans';

% create force vector and find index. The actual value is overwritten later
fex = 1;
add_forcing(beam,inode,dir,fex);
idx_f = find(beam.Fex1);

% Excitation levels
exc_lev = fliplr([0.5]);
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
    ds = 50;
    % Set options of path continuation
    % flag = 0: no actual continuation, just sequential continuation
    % stepadatp = 0: no automatic step size adjustment
    Dscale = [1e-6*ones(length(x0),1);(Om_s+Om_e)/2];
    Sopt = struct('stepadatp', 1,'errmax',4,...%'eps',1e-8,...
        'stoponerror',1,'noconv_stepcor','red',...
        'Dscale',Dscale,'dynamicDscale',1,'jac','full','stepmax',1e4);
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

% % Determine linear references
% OM = linspace(min(OM_HB),max(OM_HB),1e3);
% Qtip_fixed = zeros(size(OM));
% Qtip_free = zeros(size(OM));
% for i=1:length(OM)
%     Q1tmp = B*((B'*(-OM(i)^2*beam.M + 1i*OM(i)*beam.D + beam.K)*B)\...
%         (B'*Fex1));
%     Qtip_fixed(i) = abs(T_tip*Q1tmp)/sqrt(2);
%     Q1tmp = (-OM(i)^2*beam.M + 1i*OM(i)*beam.D + beam.K)\...
%         (Fex1);
%     Qtip_free(i) = abs(T_tip*Q1tmp)/sqrt(2);
% end

%% NMA results (for backbone calculation)

nma = load('nma.mat');
H_nma = nma.H;
X_HB_nma = nma.X_HB;
Psi_HB = X_HB_nma(1:end-3,:);
om_HB_nma = X_HB_nma(end-2,:);
del_HB_nma = X_HB_nma(end-1,:);
log10a_NMA = X_HB_nma(end,:);
a_NMA = 10.^log10a_NMA;
Q_HB_nma = Psi_HB.*repmat(a_NMA,size(Psi_HB,1),1);

% Determine harmonics and root-mean-square value of tip displacement
Qtip_HB_nma = kron(eye(2*H_nma+1),T_tip)*Q_HB_nma;
Qtip_rms_HB_nma = sqrt(sum(Qtip_HB_nma.^2))/sqrt(2);

%% Compare frequency response of both Harmonic Balance and Shooting methods
figure; hold on;
% plot(OM,Qtip_fixed,'-','color',.75*[1 1 1]);
plot(om_HB_nma,Qtip_rms_HB_nma,'k-');
for iex=1:length(exc_lev)
    plot(OM_HB{iex},Qtip_rms_HB{iex},'g-');
end
set(gca,'yscale','log');
xlabel('excitation frequency'); ylabel('tip displacement amplitude');
legend('fixed contact','HB');


%% Comparison against direct time step integration

% Select resonance
[~,ind] = max(Qtip_rms_HB{1});
Om = OM_HB{1}(ind);

% Synthesize time history of harmonic balance
tau = linspace(0,2*pi,N+1)';
DFT = [ones(N+1,1) zeros(N+1,2*H)];
DFT(:,2:2:end) = cos(tau*(1:H));
DFT(:,3:2:end) = sin(tau*(1:H));
q_HB = DFT*reshape(Q_HB{1}(:,ind),beam.n,2*H+1)';
DER = zeros(2*H,1);
DER(2:2:end) = 1:H;
DER = ( diag(DER,1) - diag(DER,-1) )*Om;
u_HB = DFT*DER*reshape(Q_HB{1}(:,ind),beam.n,2*H+1)';

% Time interval and discretization
ts = 0; te = 2*pi/Om*20;
maxstep = 2*pi/Om/20;%2*pi/om_sticking(end)/5;

% Initial values
q0 = zeros(beam.n,1);
u0 = zeros(beam.n,1);

% Translate properties
M = beam.M;
D = beam.D;
K = beam.K;

% Run simulation using ODE45
tic
[t,Y] = ode45(@(t,y) [y(n+1:2*n);M\(Fex1*cos(Om*t) - K*y(1:n) - ...
    D*y(n+1:2*n) - T_nl'*elastic_dry_friction(kt,muN,T_nl*y(1:n)) )],...
    [ts te],[q0;u0]);%,odeset('maxstep',maxstep));
disp(['Conventional direct time step integration for a single (!) ' ...
    'frequency from homogeneous initial conditions required ' ...
    num2str(toc) ' s.']);

% Select results of last few periods
ILP = t>(t(end)-5*2*pi/Om);

% Phase portrait of tip deflection
figure; hold on;
plot(T_tip*q_HB',T_tip*u_HB','g-','linewidth',2);
plot(T_tip*Y(ILP,1:n)',T_tip*Y(ILP,n+1:2*n)','k--');
xlabel('q_{tip}'); ylabel('u_{tip}');
legend('HB','time integration');

% Phase portrait of deflection at nonlinearity location
figure; hold on;
plot(T_nl*q_HB',T_nl*u_HB','g-','linewidth',2);
plot(T_nl*Y(ILP,1:n)',T_nl*Y(ILP,n+1:2*n)','k--');
xlabel('q_{fric}'); ylabel('u_{fric}');
legend('HB','time integration');

function fnl=elastic_dry_friction(k,fc,qnl)
% Set stiffness of elastic dry friction element (aka
% Jenkins element)

% Initialize force vector and Coulomb slider position
fnl = zeros(size( qnl ));
qsl = zeros(size( qnl ));

% Iterative computation of force
for ij = 2:length(fnl)
    % Predictor step
    fnl(ij) = k*( qnl(ij) - qsl(ij-1) );
    
    % Corrector step
    if(abs(fnl(ij)) >= fc)
        fnl(ij) = fc*sign(fnl(ij));
        qsl(ij) = qnl(ij) - fnl(ij)/k;
    else
        qsl(ij) = qsl(ij-1);
    end
end
end
