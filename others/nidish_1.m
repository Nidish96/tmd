clc
clear all

addpath('./NLvib_v1.1/SRC/')
addpath('./NLvib_v1.1/SRC/MechanicalSystems/')

%% System definition
len = 0.70;
hgt = 0.03;
thk = hgt;
E   = 185e9;
rho = 7830.0;
BCs = 'clamped-free';

Nn = 8;
beam = FE_EulerBernoulliBeam(len, hgt, thk, E, rho, BCs, Nn);

% Nonlinearity
Nnl = 4;
dir = 'trans';
kt  = 1.3e6;
muN = 1;
add_nonlinear_attachment(beam, Nnl, dir, 'elasticdryfriction', ...
    'stiffness', kt, 'friction_limit_force', muN, ...
    'ishysteretic', true);

%% Linearized limit cases
% Slipped
[Vsl, Dsl] = eig(beam.K, beam.M);
[Dsl, si] = sort(sqrt(diag(Dsl)));
Vsl = Vsl(:, si);  Vsl = Vsl./sqrt(diag(Vsl'*beam.M*Vsl)');

% Stuck
Knl = zeros(size(beam.M));
nli = find(beam.nonlinear_elements{1}.force_direction);
Knl(nli, nli) = kt;
Kst = beam.K + Knl;
[Vst, Dst] = eig(Kst, beam.M);
[Dst, si] = sort(sqrt(diag(Dst)));
Vst = Vst(:, si); Vst = Vst./sqrt(diag(Vst'*beam.M*Vst));

%% Rayleigh damping
% Desired
zs = [8e-3; 2e-3];
ab = [ones(length(zs),1) Dst(1:length(zs)).^2]\(2*zs.*Dst(1:length(zs)));

beam.D = ab(1)*beam.M + ab(2)*Kst;

Zetas = diag(Vst'*beam.D*Vst)./(2*Dst);

%% Frequency Response
analysis = 'FRF';
Nh = 9;
Nhc = 2*Nh+1;
Nd = size(beam.M, 1);
Nt = 2^10;

Ws = 200;
We = 500;
ds = (We-Ws)/100;
dsmax = (We-Ws)/50;

Dscale = [1e-4*ones(Nd*Nhc,1); (Dst(1)+Dsl(1))/2];

Fas = [0.01 0.05 0.1 0.2 0.3 0.6];
Xcont = cell(size(Fas));
for k=1:length(Fas)
    beam.Fex1(end-1) = Fas(k);  % Forcing from last node
    H1tmp = ((Kst-Ws^2*beam.M)+1j*(Ws*beam.D))\beam.Fex1;
    X0 = zeros(Nd*Nhc, 1); X0(Nd+(1:2*Nd)) = [real(H1tmp); -imag(H1tmp)];
    
    Sopt = struct('jac','full','stepmax',1000,'MaxFfunEvals',500, ...
        'dsmax', dsmax, 'Dscale', Dscale, 'dynamicDscale', 1); %, ...
%         'parametrization', 'normal');
    Xcont{k} = solve_and_continue(X0, ...
        @(X) HB_residual(X, beam, Nh, Nt, analysis), Ws, We, ds, Sopt);
end

%% EPMC NMA
analysis = 'nma';
imod = 1;
log10a_s = -6;
log10a_e = -3.8;
inorm = Nd-1;

X0 = zeros(Nd*Nhc+2, 1);
X0(Nd+(1:Nd)) = Vst(:,imod);
X0(end-1) = Dst(imod);
X0(end) = 2*Zetas(imod)*Dst(imod);

Dscale = [1e-4*ones(Nd*Nhc,1); (Dst(imod)+Dsl(imod))/2; 2*Zetas(imod)*Dst(imod); 1.0];

beam.Fex1 = beam.Fex1*0;

ds = 0.05;
dsmax = abs(log10a_e-log10a_s)/5;
Sopt = struct('jac','full','stepmax',1000,'MaxFfunEvals',500, ...
	'dsmax', dsmax); %, ...
%         'parametrization', 'normal');
fscl = mean(abs(beam.K*Vst(:,imod)));
Xbb = solve_and_continue(X0, ...
    @(X) HB_residual(X, beam, Nh, Nt, analysis, inorm, fscl), ...
    log10a_s, log10a_e, ds, Sopt);

%% Plotting Responses
figure(1)
clf()
aa = gobjects(length(Fas)+1,1);
for k=1:length(Fas)
    aa(k) = semilogy(Xcont{k}(end,:)/(2*pi), sqrt([1 0.5*ones(1,Nhc-1)]*Xcont{k}((Nd-1):Nd:end-1,:).^2), ...
        '-', 'LineWidth', 2); hold on
    legend(aa(k), sprintf('F = %.2f N', Fas(k)));
end
aa(end) = semilogy(Xbb(end-2,:)/(2*pi), 10.^Xbb(end,:).*sqrt([1 0.5*ones(1,Nhc-1)]*Xbb((Nd-1):Nd:end-2,:).^2), ...
        'k--', 'LineWidth', 2)
legend(aa(end), 'EPMC')
legend(aa(1:end), 'Location', 'northeast')
xlabel('Frequency (Hz)')
ylabel('RMS Response Amplitude (m)')
axis tight
savefig('./FIGS/FRESP_AMP.fig')
% print('./FIGS/FRESP_AMP.eps', '-depsc')

figure(2)
clf()
for k=1:length(Fas)
    plot(Xcont{k}(end,:)/(2*pi), atan2d(-Xcont{k}(2*Nd+(Nd-1),:),Xcont{k}(Nd+(Nd-1),:)), ...
        '-', 'LineWidth', 2); hold on
end
% plot(Xbb(end-2,:)/(2*pi), atan2d(-Xbb(2*Nd+(Nd-1),:),Xbb(Nd+(Nd-1),:)), ...
%         'k--', 'LineWidth', 2); hold on
xlabel('Frequency (Hz)')
ylabel('Response Phase (degs)')
axis tight
savefig('./FIGS/FRESP_PHASE.fig')
% print('./FIGS/FRESP_PHASE.eps', '-depsc')

%% Deflection Shapes
k = 1;

for k=1:size(Xbb,2)
tmp = [zeros(2,Nhc); reshape(Xbb(1:end-3,k), Nd, Nhc)];
phi = atan2(tmp(end-1,3), tmp(end-1,2));
Defs = tmp(:, 1) + sum(tmp(:, 2:2:end)*cos(phi)+tmp(:, 3:2:end)*sin(phi),2);

Xnodal = linspace(0, len, Nn);

figure(3)
clf()
plot(Xnodal, Defs(1:2:end), 'o-')
xlabel('X Coordinate (m)')
ylabel('Y Deflection (m)')
title(sprintf('Frame %d/%d', k, size(Xbb,2)))
end