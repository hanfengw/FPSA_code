% Hanfeng Wang 12.23.2022
% For "Field Programmable Spin Arrays for Scalable Quantum Repeaters" paper

% This code is for cross-talk elimination process for the strain-based
% FPSA. We perform a Green function-like way to make the cross-talk
% elimination process as discussed in the main text. This code is for both
% 1D and 2D. 

%%

beta_norm(:,1) = beta1Dv3(:,1)/244e-9;
beta_norm(:,2) = beta1Dv3(:,2)/1e9;

pp = 0;

beta_point(:,1) = 1:1:36;

for ii = 1:length(beta_norm(:,1)) 
    
    if beta_norm(ii,1)>=pp+0.5
        beta_point(pp+1,2) = beta_norm(ii,2);
        pp = pp + 1;
    end

end


gamma_norm(:,1) = gamma1Dv3(:,1)/244e-9;
gamma_norm(:,2) = gamma1Dv3(:,2)/1e9;

pp = 0;

gamma_point(:,1) = 1:1:36;

for ii = 1:length(gamma_norm(:,1)) 
    
    if gamma_norm(ii,1)>=pp+0.5
        gamma_point(pp+1,2) = gamma_norm(ii,2);
        pp = pp + 1;
    end

end






%%

for ii = 1:14

    Array(ii,1:7) = beta_point(20-ii:2:32-ii,2);
    
    Array(ii,8:14) = gamma_point(20-ii:2:32-ii,2);

end


% for ii = 1:10
% 
%     Array(ii+10,1:10) = beta_point2(13-ii:1:22-ii,2);
%     
%     Array(ii+10,11:20) = gamma_point2(13-ii:1:22-ii,2);
% 
% end

map = inv(Array');

%%

map1d = map(:,4);

nq = 1706;

xx = linspace(11,25,nq)-17.5;



sum1_1D = zeros(1,nq);

sum2_1D = zeros(1,nq);

beta_inter = zeros(14,nq);

gamma_inter = zeros(14,nq);

for ii = 1:14

    beta_inter(ii,:) = interp1(beta_norm(:,1),beta_norm(:,2),linspace(11,25,nq)+8-ii);

    gamma_inter(ii,:) = interp1(gamma_norm(:,1),gamma_norm(:,2),linspace(11,25,nq)+8-ii);



    sum1_1D = sum1_1D + map1d(ii) * beta_inter(ii,:);

    sum2_1D = sum2_1D + map1d(ii) * gamma_inter(ii,:);

end

semilogy(xx(1:2:end)/2,85*sqrt(beta_inter(8,1:2:end).^2+gamma_inter(8,1:2:end).^2),'linewidth',2,'color',[0.8500 0.3250 0.0980])

hold on

semilogy(xx(1:2:end)/2,3.27*85*sqrt(sum1_1D(1:2:end).^2+sum2_1D(1:2:end).^2),'linewidth',2,'color',[0 0.4470 0.7410])

set(gca,'xlim',[-6.2 6.2]/2)

hold on

semilogy([-6.2 6.2],[11 11],'--','linewidth',2)




set(gca,'Fontsize',15)

% hold on
% 
% plot(linspace(0,20,2395)+1,beta_inter)

legend('Without CE','With CE','F = 0.99')

grid on

set (gca,'position',[0.1,0.1,0.8,0.5]);

set(gca,'ylim',[0.01 8000])

xticks([-12:2:12])

yticks([0.01 1 100])


%%



beta2D_norm(:,1) = beta2Dv3(:,1)/244e-9;
beta2D_norm(:,2) = beta2Dv3(:,2)/244e-9;
beta2D_norm(:,3) = beta2Dv3(:,3)/1e9;

gamma2D_norm(:,1) = gamma2Dv3(:,1)/244e-9;
gamma2D_norm(:,2) = gamma2Dv3(:,2)/244e-9;
gamma2D_norm(:,3) = gamma2Dv3(:,3)/1e9;


% [xq,yq] = meshgrid(-6:1e-2:6, -3:1e-2:3);
% vq = griddata(beta2D_norm(:,1),beta2D_norm(:,2),beta2D_norm(:,3),xq,yq);

sum1 = zeros(601,1201);
sum2 = zeros(601,1201);


for ii = 1:14
    ii

    [xq,yq] = meshgrid((-6:1e-2:6)+8-ii, (-3:1e-2:3));

    vq11 = griddata(beta2D_norm(:,1),beta2D_norm(:,2),beta2D_norm(:,3),xq,yq);

    vq22 = griddata(gamma2D_norm(:,1),gamma2D_norm(:,2),gamma2D_norm(:,3),xq,yq);
 
%     beta_inter(ii,:) = interp1(beta_norm(:,1),beta_norm(:,2),linspace(5,15,1255)+6-ii);
% 
%     gamma_inter(ii,:) = interp1(gamma_norm(:,1),gamma_norm(:,2),linspace(5,15,1255)+6-ii);



    sum1 = sum1 + map(ii,4) * vq11;
    sum2 = sum2 + map(ii,4) * vq22;



end

%% Fig2e

[xq,yq] = meshgrid(-6:1e-2:6, -0.84:1e-2:0.84);

A = sqrt(sum1.^2+sum2.^2);





h = pcolor(xq+1,yq,3.27*A(217:385,:)*85);

set(h, 'EdgeColor', 'none');



shading interp

colormap hot



colorbar

caxis([1 1000])

set(gca,'xlim',[-5 5.5])
set(gca,'Fontsize',15)
set(gca, 'ColorScale', 'log')


% set (gca,'position',[0.1,0.1,0.8,0.2]);

%%

[xq,yq] = meshgrid(-6:1e-2:6, -0.84:1e-2:0.84);



vq = griddata(beta2D_norm(:,1),beta2D_norm(:,2),beta2D_norm(:,3),xq,yq);

vq2 = griddata(gamma2D_norm(:,1),gamma2D_norm(:,2),gamma2D_norm(:,3),xq,yq);

%%

h = pcolor(xq+1,yq,sqrt(vq.^2+vq2.^2)*85);

set(h, 'EdgeColor', 'none');



shading interp

colormap hot



colorbar

caxis([1 1000])

set(gca,'xlim',[-5 5.5])
set(gca,'Fontsize',15)
set(gca, 'ColorScale', 'log')





%%

plot(es12_norm(:,1),es12_norm(:,2))

%%

% plot(betaafterCE(:,1)/244e-9,log10(sqrt(betaafterCE(:,2).^2 + gammaafterCE(:,2).^2 )/1e9),'linewidth',1.5)

semilogy(betagamma(1:2:end,1)/244e-9-10,betagamma(1:2:end,2)/1e9*85*10,'linewidth',2)

hold on

semilogy(betagammav3(1:2:end,1)/244e-9-10,betagammav3(1:2:end,2)/1e9*85*3.5201,'linewidth',2)




set(gca,'xlim',[-7 7])

set(gca,'ylim',[10^(-0.5) 1e4])

set(gca,'Fontsize',15)

legend('Without CE','With CE')

%%

plot(betasquare(:,1)/244e-9,log10(sqrt(betasquare(:,2))/1e9),'linewidth',1.5)

set(gca,'xlim',[5 16])




