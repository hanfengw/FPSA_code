
% Hanfeng Wang 12.23.2022
% For "Field Programmable Spin Arrays for Scalable Quantum Repeaters" paper
% This code is for dealing with the data for before and after cross-talk
% elimination. For the after cross-talk elimination part, due to the
% symmetry consideration, we only consider half of the data and flip them
% to get another half. The simulation itself looks not symmetric due to the
% simulation error for Ex and Ez and the linecut through the mesh. 


before_norm(:,1) = (before(:,1)-0.001)/182e-6;
before_norm(:,2) = before(:,2);

aft_norm2(:,1) = (after2(:,1)-0.001)/182e-6;
aft_norm2(:,2) = after2(:,2);

before_norm_2(:,1) = before_norm(:,1);
before_norm_2(:,2) = before_norm(:,2)/1e6;

semilogy(before_norm_2(:,1),before_norm_2(:,2),'linewidth',2,'color',[0.8500 0.3250 0.0980])

hold on

clear aft_v2

aft_v2(1:4039,1) = aft_norm2(1:4039,1);
aft_v2(1:4039,2) = aft_norm2(1:4039,2);

aft_v2(4039:8077,1) = -flip(aft_norm2(1:4039,1));
aft_v2(4039:8077,2) = flip(aft_norm2(1:4039,2));   


aft_v2_2(:,1) = aft_v2(:,1);
aft_v2_2(:,2) = aft_v2(:,2)/1e6*25/0.5161;  %% The comsol data is normalized by voltage, change the voltage here

semilogy(aft_v2_2(:,1),aft_v2_2(:,2),'linewidth',2,'color',[0 0.4470 0.7410])

semilogy([-4.5 4.5],[1.7 1.7],'--','linewidth',2)

set(gca,'Fontsize',14)


legend('Without CE','With CE','F = 0.99')

grid on

set(gca,'xlim',[-4.5 4.5])

set(gca,'ylim',[1e-5 4e3])

yticks([1e-4 1e-2 1e0 1e2])

xticks([-4:1:4])

box off

ax1 = gca;

ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
     'XAxisLocation','top',...
     'YAxisLocation','right',...
     'Color','none');

ax2.XLim = [-4.5 4.5];

ax2.YLim = [1e-6 4e3]*0.17;

ax2.YScale = 'log';

ax2.YTick = [1e-6 1e-4 1e-2 1e0 1e2];


