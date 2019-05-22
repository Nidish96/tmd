
close all
clearvars

% multisine formulated in freq domain. Requres fixed step solver.

%% Define system
% Fundamental parameters
Dmod = [.38 .12 .09 .08 .08]*.01;
Nmod = 5;
setup = 'New_Design_Steel';
thickness = .001;
[L,rho,E,om,PHI,~,gam] = beams_for_everyone(setup,Nmod,thickness);
PHI_L2 = PHI(L/2);

% load nonlinear coefficients (can be found e.g. analytically)
fname = ['beam_New_Design_Steel_analytical_5t_' ...
    num2str(thickness*1000) 'mm.mat'];
[p, E] = nlcoeff(fname, Nmod);

% Properties of the underlying linear system
M = eye(Nmod);
D = diag(2*Dmod(:).*om(:));
K = diag(om.^2);

% Fundamental harmonic of external forcing
Fex1 = gam;

nz = size(p,1);
n = Nmod;

%% Time integration parameters.
fs = 2000;               % working sampling frequency.
fsint = 1e7;
upsamp = fsint/fs;

% upsamp = 20;            % upsampling factor to ensure integration accuracy.
% fsint = fs*upsamp;      % integration sampling frequency.
h = 1/fsint;            % integration time step.

%% Excitation signal design. 
P = 1;                  % number of excitation periods.
N = 8192;               % number of points per period.
Nint = N*upsamp;        % number of points per period during integration.

fmin = 5;               % excitation bandwidth.
fmax = 400;

Amp = 5;                 % excitation amplitude.
Pfilter = 1;            % extra period to avoid edge effects during low-pass filtering (see line 59).
P = P + Pfilter;

U = zeros(Nint,1);      % definition of the multisine excitation.
fres = fsint/Nint;
exclines = 1:ceil(fmax/fres);
exclines(exclines < floor(fmin/fres)) = [];

U(exclines+1) = exp(complex(0,2*pi*rand(size(exclines))));
u = 2*real(ifft(U));
u = Amp*u/std(u);
u = repmat(u,[P 1]);

famp = Fex1*frf.exc_lev(1);
par = struct('M',M,'C',D,'K',K,'p',p,'E',E,'fex',u, 'amp', famp);
tic
Nt = length(u);
t = (0:Nt-1)*h;

% y = ode2_modif(@(t,y) sys(t,y, par), h,[q0;u0], Nt);
                        %,odeset('maxstep',maxstep));
disp(['Conventional direct time step integration for a single (!) ' ...
    'frequency from homogeneous initial conditions required ' ...
    num2str(toc) ' s.']);


%% Low-pass filtering and downsampling.
drate = factor(upsamp);        % prime factor decomposition.
for k=1:length(drate),
    y = decimate(y,drate(k),'fir');    
end %k
u = downsample(u,upsamp);

%% Removal of the last simulated period to eliminate the edge effects due to the low-pass filter.
y = y(1:(P-1)*N);
u = u(1:(P-1)*N,:);
P = P-1;


