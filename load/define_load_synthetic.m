loads = struct;

%% OPTION 1: simulate synthetic load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % harmonic load at top of structure 
% params = [100, 0.2];
% dofs = 10;
% loads.p1 = load_Harmonic(params(1), params(2), dt_f, T, 10);
% 
% 
% % GP distributed along height
% hp = [50, 1];
% dofs = 1:ndof;
% loads.p2 = load_Matern52(hp, dt_f, T, dofs);
% 
% 
% % combine 
% np = length(fieldnames(loads));
% nsteps = length(loads.p1.time);
% 
% % define force time history
% loads.f = zeros(np, nsteps);
% loads.f(1,:) = loads.p1.force;
% loads.f(2,:) = loads.p2.force;
% 
% % define load indicator matrices
% Sp = zeros(ndof, np);
% Sp(loads.p1.loc, 1) = 1;
% Sp(loads.p2.loc, 2) = 1;


%% OPTION 2: load from file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("sample_input.mat"); % loads, Sp, nsteps


%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% combined load at top of tower
time = dt_f:dt_f:nsteps*dt_f;
comb_f = loads.f(1,:);
for i = 2:size(loads.f,1)
    comb_f = comb_f + loads.f(i,:);
end


fig = figure;
set(gcf, 'Position', [100, 100, 650, 400]);

% time series
sp1 = subplot(2,2,1);
sp1.Position = sp1.Position + [0 0 0.1 0];
plot(time, comb_f, 'k'); hold on 
xlim([0 T]);
xlabel('time (t)');
ylabel('f(t)');

% histogram
sp2 = subplot(2,2,2);
sp2.Position = sp2.Position + [0.1 0 -0.1 0];
h = histogram(comb_f, 30, 'FaceColor', 'k');
% h.Normalization = 'pdf';
h.Orientation= 'horizontal';
xlabel('count');
ylabel('f(t)');

% frequency
subplot(2,2,[3,4]);
win = hann(fix((length(time))/4));
[psd,f] = pwelch(comb_f,win,fix(length(win)/2),[],fs_f,'power');
semilogy(f, psd, 'k', 'Linewidth', 1.2); hold on
xlabel('freq. (Hz)'); xlim([0 2]);
ylabel('PSD');


