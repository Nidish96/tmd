% Compare periodic solution with that found by HB
%
% This requires that the odesys is written in such a way that it can handle
% adaptive time steps. This is not the case if the multisine is calculated
% a priori, fx. using PNLSS
%
% We compare the periodic solutions found variable step ode45 and fixed
% step ode2(modified Euler's method, second order)


close all
clearvars

srcpath = '~/src/matlab/nlvib';
addpath(genpath(srcpath));
srcpath = '../';
addpath(srcpath);

figpath = '../fig/';

%% Define system
% Fundamental parameters
Dmod = [.38 .12 .09 .08 .08]*.01;
Nmod = 5;
setup = 'New_Design_Steel';
thickness = .001;
[L,rho,E,om,PHI,~,gam] = beams_for_everyone(setup,Nmod,thickness);
PHI_L2 = PHI(L/2);

% Properties of the underlying linear system
M = eye(Nmod);
D = diag(2*Dmod(:).*om(:));
K = diag(om.^2);

% load nonlinear coefficients (can be found e.g. analytically)
fname = ['beam_New_Design_Steel_analytical_5t_' ...
    num2str(thickness*1000) 'mm.mat'];
[p, E] = nlcoeff(fname, Nmod);

% Fundamental harmonic of external forcing
Fex1 = gam;

nz = size(p,1);
n = Nmod;

%% select periodic solution at fold bifurcation (resonance/upper peak).
N = 1e3;
frf = load_hb('frf.mat',N,PHI_L2);
H = frf.H;

% Select resonance
[~,ind] = max(frf.arms{1});
Om = frf.Om{1}(ind);

% Synthesize time history of harmonic balance
tau = linspace(0,2*pi,N+1)';
DFT = [ones(N+1,1) zeros(N+1,2*H)];
DFT(:,2:2:end) = cos(tau*(1:H));
DFT(:,3:2:end) = sin(tau*(1:H));
q_HB = DFT*reshape(frf.Qsc(:,ind),n,2*H+1)';
DER = zeros(2*H,1);
DER(2:2:end) = 1:H;
DER = ( diag(DER,1) - diag(DER,-1) )*Om;
u_HB = DFT*DER*reshape(frf.Qsc(:,ind),n,2*H+1)';


%% periodic ode solution

% Time interval and discretization - 50 periods. Change if you change init
% cond. At least remember to check that we reached steady state!
ts = 0; te = 2*pi/Om*50;
maxstep = 2*pi/Om/20;

% Initial values - Important! Ensure we reach the right steady state
% solution - and not the lower one, which will happen if we use zero as
% initial guess. As a bonus, we reach steady state fast!
q0 = q_HB(1,:)';
u0 = u_HB(1,:)';
% q0 =  zeros(n,1);
% u0 =  zeros(n,1);

famp = Fex1*frf.exc_lev(1);
fex = @(t) cos(Om*t);
% Run simulation using variable time-step RK45(dopri5)
par = struct('M',M,'C',D,'K',K,'p',p,'E',E,'fex',fex, 'amp', famp);
tic
[t,Y] = ode45(@(t,y) sys(t,y, par), [ts te],[q0;u0]);
                        %,odeset('maxstep',maxstep));
disp(['ode45 required ' num2str(toc) ' s.']);

% Select results of last few periods
ILP = t>(t(end)-5*2*pi/Om);

%% fixed step ode3
tic
Nt = 1e5;
h = 1/1e6;
t2 = (0:Nt-1)*h;

% we need to give the time step h as input to cos; because ode2_modif calls
% using time index instead of actual time.
fex = @(t) cos(h*Om*t);
par = struct('M',M,'C',D,'K',K,'p',p,'E',E,'fex',fex, 'amp', famp);

Y2 = ode2_modif(@(t,y) sys(t,y, par), h,[q0;u0], Nt);
disp(['ode2 required ' num2str(toc) ' s.']);

% Select results of last few periods
ILP2 = t2>(t2(end)-3*2*pi/Om);

%% visually check we reached steady state
figure;
plot(t, PHI_L2 *Y(:,1:n)','k-');
xlabel('t'); ylabel('q');
legend('time integration');
% export_fig([figpath, 'ode45_resonance.pdf'])

figure;
plot(t2, PHI_L2 *Y2(:,1:n)','k-');
xlabel('t'); ylabel('q');
legend('time integration');
% export_fig([figpath, 'ode2_resonance.pdf'])

%% Phase portrait of middle deflection
figure; hold on;
plot(PHI_L2 *q_HB',PHI_L2 *u_HB','g-','linewidth',2);
plot(PHI_L2 *Y(ILP,1:n)',PHI_L2 *Y(ILP,n+1:2*n)','k--','linewidth',1);
plot(PHI_L2 *Y2(ILP2,1:n)',PHI_L2 *Y2(ILP2,n+1:2*n)','r-.','linewidth',0.5);
xlabel('q'); ylabel('u');
legend('HB','ode45','ode2');
% export_fig([figpath, 'phase_portrait.pdf'])

