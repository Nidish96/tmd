% clearvars
close all

srcpath = '~/src/matlab/nlvib/SRC';
addpath(genpath(srcpath));
srcpath = '~/src/matlab/export_fig';
addpath(srcpath);

% index of Q_HB
% columns: HB coeffs for each omega step
% rows: n*(2*HB+1)'s HB coeffs. n:DOFs. Ie fundamental harmonic motion:
% Y_HB_1 = Q_HB(n+(1:n),:)-1i*Q_HB(2*n+(1:n),:);

%% HB results for NFRC
hb = load('hb.mat');
OM_HB = hb.OM_HB;
Q_HB = hb.Q_HB;
Qtip_HB = hb.Qtip_HB;
Qtip_rms_HB = hb.Qtip_rms_HB;
beam = hb.beam;
exc_lev = hb.exc_lev;
H = hb.H;
N = hb.N;
ds = hb.ds;

% Vectors recovering deflection at tip and nonlinearity location
T_nl = beam.nonlinear_elements{1}.force_direction';
T_tip = beam.L(end-1,:);
Nhc = 2*H+1;
n_dof = size(beam.M, 1);
idx_tip = find(T_tip);

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

%% set default line style
set(groot,'defaultAxesColorOrder',[1 0 0;0 1 0;0 0 1],...
    'defaultAxesLineStyleOrder','-|--|:')
% set(gca,'ColorOrder',[1 0 0;0 1 0;0 0 1],...
%     'LineStyleOrder','-|--|:')

%% Frequency response of Harmonic Balance and backbone from NMA
figure; hold on;
% plot(OM,Qtip_fixed,'-','color',.75*[1 1 1]);
plot(om_HB_nma,Qtip_rms_HB_nma,'k--', 'LineWidth',2);

slabel = cell(1,length(exc_lev));
for k=1:length(exc_lev)
    plot(OM_HB{k},Qtip_rms_HB{k}, 'LineWidth',2);
    slabel{k} = sprintf('F = %0.2f N',exc_lev(k));
end
set(gca,'yscale','log');
xlabel('excitation frequency (rad/s)');
ylabel('tip displacement amplitude (m)');
set(gca,'xlim',[150 650],'ylim',[10^(-8) 10^(-4)]);
legend(['nm - backbone',slabel],'Location','ne');
title(sprintf('HB par: H:%d, N:%d, ds:%0.2f',H,N,ds))
axis tight
export_fig('fig/cantilever_frf.pdf')

%% Phase plot
% phase between fundamental harmonic motion, ie c1*cos+c2*isin
% phase = atan2(c2,c1)

figure;
hold on
for k=1:length(exc_lev)
    plot(OM_HB{k}, atan2d(-Q_HB{k}(2*n_dof+idx_tip,:),Q_HB{k}(1*n_dof+idx_tip,:)), ...
        'LineWidth', 2)
end
xlabel('Frequency (rad/s)')
ylabel('Response Phase (degs)')
legend([slabel],'Location','ne');
axis tight
export_fig('fig/cantilever_phase.pdf')

%% Deflection Shapes
Nhc_nma = H_nma*2+1;
% Properties of the beam
len = 0.7;              % length in m
n_nodes = 8;

Xnodal = linspace(0, len, n_nodes);
% for each omega step, plot latteral deflection.
k = 1;
for k=1:size(Psi_HB,2)
    tmp = [zeros(2,Nhc_nma); reshape(Psi_HB(:,k), n_dof, Nhc_nma)];
    phi = atan2(tmp(idx_tip,3), tmp(idx_tip,2));
    Defs = tmp(:, 1) + sum(tmp(:, 2:2:end)*cos(phi) + ...
                           tmp(:,3:2:end)*sin(phi),2);

    figure(3)
    clf()
    plot(Xnodal, Defs(1:2:end), 'o-')
    xlabel('X Coordinate (m)')
    ylabel('Y Deflection (m)')
    title(sprintf('Freq %0.2f - Frame %d/%d', om_HB_nma(k),k, size(Psi_HB,2)))
end

%% remove default linestyle
set(groot,'defaultAxesLineStyleOrder','remove')
set(groot,'defaultAxesColorOrder','remove')
