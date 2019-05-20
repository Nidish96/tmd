
close all
clearvars

srcpath = '~/src/matlab/nlvib/SRC';
addpath(genpath(srcpath));


% Fundamental parameters
Dmod = [.38 .12 .09 .08 .08]*.01;
Nmod = 5;
setup = 'New_Design_Steel';
thickness = .001;
[L,rho,E,om,PHI,~,gam] = beams_for_everyone(setup,Nmod,thickness);
PHI_L2 = PHI(L/2);

% time points in period
Ntd = 1e3;
% dimensionless time - exclude last point, ie. tau = [0, 2pi[
tau = linspace(0,2*pi,Ntd+1); tau = tau(1:end-1);


set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 

%% FRF
frf = load('frf.mat');

n = frf.n;
% different ext. levels
nex = length(frf.exc_lev);

a_rms = cell(1,nex);
a_diff = cell(1,nex);
Om_frf = cell(1,nex);
for j = 1:nex
    % Interpret solver output
    Om_frf{j} = frf.X{j}(end,:);
    Q_HB = frf.X{j}(1:end-1,:);
    
    % two ways to define amplitude (at x = L/2)
    % 1) difference between min/max amplitude in time domain or
    % 2) RMS, as found by persevals identity, eq. 2.18 p 21 in (Malte book)
    % https://sci-hub.tw/10.1007/978-3-030-14023-6
    
    % w_L2: [Time points, cont. step]
    w_L2 = zeros(Ntd,length(Om_frf{j}));
    a_rms_loc = 0;
    % loop over all modes(or DOF in the discrete case)
    for k = 1:n
        % 1) convert Q_sc -> Q_ce for finding the periodic sol in time
        % include all harmonics
        Qc = [Q_HB(k,:);Q_HB(n+k:2*n:end,:)-1i*Q_HB(2*n+k:2*n:end,:)];
        % get displacement at center caused by each mode in time domain
        w_mode = PHI_L2(k)*real(exp(1i*tau(:)*(0:frf.H))*Qc);
        % sum up to actual displacement at center in time domain
        w_L2 = w_L2 + w_mode;

%         q = real(exp(1i*tau(:)*(0:frf.H))*Qc);
%         a_q(k,:) = (max((q))-min((q)))/2;

        % 2) RMS. include all harmonics
        Q_rms = sqrt(sum(Q_HB(k:n:end,:).^2))/sqrt(2);
        a_rms_loc = a_rms_loc + PHI_L2(k)*Q_rms;
    end
    % compute peak to peak amplitude
    a_diff{j} = (max((w_L2))-min((w_L2)))/2;     
    a_rms{j} = a_rms_loc;
end


% % Illustrate frequency response
% figure; hold on; box on
% aa = gobjects(nex,1);
% for j = 1:nex
%     aa(j) = plot(Om_frf{j}/2/pi,a_diff{j}*1000,'k-','linewidth',1.5);
%     legend(aa(j), sprintf('F = %.2f N', frf.exc_lev(j)));
% end
% legend(aa, 'Location', 'northeast')
% xlabel('$f_{\mathrm{ex}}$ in Hz');
% ylabel('$\hat{w}_{L/2}$ in mm');
% title(['thickness = ' num2str(thickness*1000) 'mm']);

%% NMA
nma = load('nma.mat');

Psi = nma.X(1:end-3,:);
Om_nma = nma.X(end-2,:);
del_HB = nma.X(end-1,:);
log10a_NMA = nma.X(end,:);
a_NMA = 10.^log10a_NMA;
Q_HB = Psi.*repmat(a_NMA,size(Psi,1),1);

w_L2_nma= zeros(Ntd,length(Om_nma)); 
a_q_NMA = zeros(n,length(Om_nma));
a_rms_nma = 0;
% loop over all modes
for k = 1:n
    % get complex coefficients
    Qc = [Q_HB(k,:);Q_HB(n+k:2*n:end,:)-1i*Q_HB(2*n+k:2*n:end,:)];
    w_mode = PHI_L2(k)*real(exp(1i*tau(:)*(0:nma.H))*Qc);
    w_L2_nma= w_L2_nma + w_mode;
    
    q_NMA = real(exp(1i*tau(:)*(0:nma.H))*Qc);
    a_q_NMA(k,:) = (max((q_NMA))-min((q_NMA)))/2;
    
    Q_rms = sqrt(sum(Q_HB(k:n:end,:).^2))/sqrt(2);
    a_rms_nma = a_rms_nma + PHI_L2(k)*Q_rms;
end

a_diff_nma= (max((w_L2_nma))-min((w_L2_nma)))/2;

% Determine total energy in the system from the displacement and velocity
% at t=0
energy = zeros(size(a_NMA));
% loop over countinuation parameter(omega)
for i=1:length(Om_nma)
    Qi = reshape(Q_HB(:,i),n,2*H+1);
    q0 = Qi(:,1) + sum(Qi(:,2:2:end),2);
    u0 = sum(Qi(:,3:2:end),2)*om_HB(i);
    energy(i) = 1/2*u0'*beam.M*u0 + 1/2*q0'*beam.K*q0 + knl*q0(n-1)^4/4;
end

%% set default line style
% https://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
set(groot,'defaultAxesColorOrder',[1 0 0;0 1 0;0 0 1;0 0.5 0;0.75 0 0.75],...
    'defaultAxesLineStyleOrder','--|:|-')

%% FEP

% Modal frequency vs. amplitude
figure; hold on;
plot(log10(energy),Om_nma/(2*pi),'k-o');
xlabel('log10(energy)'); ylabel('modal frequency in Hz');
% set(gca,'ylim',[20 50]);

%% NFRC
% min/max amplitude
figure;
hold on
plot(Om_nma/2/pi,1000*abs(a_diff_nma),'color',[1 .2 .3],'linewidth',2)

slabel = cell(1,nex);
for j = 1:nex
    aa(j) = plot(Om_frf{j}/2/pi,a_diff{j}*1000,'k-','linewidth',1.5);
    slabel{j} = sprintf('F = %0.2f N',frf.exc_lev(j));
end
legend(['nm - backbone', slabel], 'Location', 'northeast')

xlim([200 400])
xlabel('$f_{\mathrm{ex}}$ in Hz');
ylabel('$\hat{w}_{L/2}$ in mm');
title(['thickness = ' num2str(thickness*1000) 'mm']);


% RMS amplitude
figure;
hold on
plot(Om_nma/2/pi,1000*abs(a_rms_nma),'color',[1 .2 .3],'linewidth',2)

slabel = cell(1,nex);
for j = 1:nex
    aa(j) = plot(Om_frf{j}/2/pi,a_rms{j}*1000,'k-','linewidth',1.5);
    slabel{j} = sprintf('F = %0.2f N',frf.exc_lev(j));
end
legend(['nm - backbone', slabel], 'Location', 'northeast')

xlim([200 400])
xlabel('$f_{\mathrm{ex}}$ in Hz');
ylabel('$RMS \hat{w}_{L/2}$ in mm');
title(['thickness = ' num2str(thickness*1000) 'mm']);

%% bending modes

figure;
plot(Om_nma/2/pi,a_q_NMA(1:5,:),'linewidth',1.5)
legend('$q_1$ (first bending)','$q_2$ (second bending)',...
    '$q_3$ (third bending)','$q_4$','$q_5$')
xlabel('$\omega$ in Hz');
ylabel('modal coordinates $q$ of NMA')
xlim([min(Om_nma/2/pi) max(Om_nma/2/pi)]);
title(['thickness = ' num2str(thickness*1000) 'mm']);

%% remove default linestyle
set(groot,'defaultAxesLineStyleOrder','remove')
set(groot,'defaultAxesColorOrder','remove')