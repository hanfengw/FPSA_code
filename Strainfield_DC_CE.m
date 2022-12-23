% Hanfeng Wang 12.23.2022
% For "Field Programmable Spin Arrays for Scalable Quantum Repeaters" paper

% This code is for cross-talk elimination process for the strain-based
% FPSA. We perform a Green function-like way to make the cross-talk
% elimination process as discussed in the main text. This code is for both
% 1D and 2D. This code is for DC. 

%%


es1_norm(:,1) = es1(:,1)/244e-9;
es1_norm(:,2) = es1(:,2)/1e9;

pp = 0;

es1_point(:,1) = 0:1:22;

for ii = 1:length(es1_norm(:,1)) 
    
    if es1_norm(ii,1)>=pp
        es1_point(pp+1,2) = es1_norm(ii,2);
        pp = pp + 1;
    end

end


es2_norm(:,1) = es2(:,1)/244e-9;
es2_norm(:,2) = es2(:,2)/1e9;

pp = 0;

es2_point(:,1) = 0:1:22;

for ii = 1:length(es2_norm(:,1)) 
    
    if es2_norm(ii,1)>=pp
        es2_point(pp+1,2) = es2_norm(ii,2);
        pp = pp + 1;
    end

end





gs1_norm(:,1) = gs1(:,1)/244e-9;
gs1_norm(:,2) = gs1(:,2)/1e9;

pp = 0;

gs1_point(:,1) = 0:1:22;

for ii = 1:length(gs1_norm(:,1)) 
    
    if gs1_norm(ii,1)>=pp
        gs1_point(pp+1,2) = gs1_norm(ii,2);
        pp = pp + 1;
    end

end


gs2_norm(:,1) = gs2(:,1)/244e-9;
gs2_norm(:,2) = gs2(:,2)/1e9;

pp = 0;

gs2_point(:,1) = 0:1:22;

for ii = 1:length(gs2_norm(:,1)) 
    
    if gs2_norm(ii,1)>=pp
        gs2_point(pp+1,2) = gs2_norm(ii,2);
        pp = pp + 1;
    end

end







es12_norm(:,1) = es12(:,1)/244e-9;
es12_norm(:,2) = es12(:,2)/1e9;

pp = 0;

es12_point(:,1) = 0:1:22;

for ii = 1:length(es12_norm(:,1)) 
    
    if es12_norm(ii,1)>=pp
        es12_point(pp+1,2) = es12_norm(ii,2);
        pp = pp + 1;
    end

end


es22_norm(:,1) = es22(:,1)/244e-9;
es22_norm(:,2) = es22(:,2)/1e9;

pp = 0;

es22_point(:,1) = 0:1:22;

for ii = 1:length(es22_norm(:,1)) 
    
    if es22_norm(ii,1)>=pp
        es22_point(pp+1,2) = es22_norm(ii,2);
        pp = pp + 1;
    end

end





gs12_norm(:,1) = gs12(:,1)/244e-9;
gs12_norm(:,2) = gs12(:,2)/1e9;

pp = 0;

gs12_point(:,1) = 0:1:22;

for ii = 1:length(gs12_norm(:,1)) 
    
    if gs12_norm(ii,1)>=pp
        gs12_point(pp+1,2) = gs12_norm(ii,2);
        pp = pp + 1;
    end

end


gs22_norm(:,1) = gs22(:,1)/244e-9;
gs22_norm(:,2) = gs22(:,2)/1e9;

pp = 0;

gs22_point(:,1) = 0:1:22;

for ii = 1:length(gs22_norm(:,1)) 
    
    if gs22_norm(ii,1)>=pp
        gs22_point(pp+1,2) = gs22_norm(ii,2);
        pp = pp + 1;
    end

end





ezz_norm(:,1) = ezz(:,1)/244e-9;
ezz_norm(:,2) = ezz(:,2)/1e9;

pp = 0;

ezz_point(:,1) = 0:1:22;

for ii = 1:length(ezz_norm(:,1)) 
    
    if ezz_norm(ii,1)>=pp
        ezz_point(pp+1,2) = ezz_norm(ii,2);
        pp = pp + 1;
    end

end




ezz2_norm(:,1) = ezz2(:,1)/244e-9;
ezz2_norm(:,2) = ezz2(:,2)/1e9;

pp = 0;

ezz_point(:,1) = 0:1:22;

for ii = 1:length(ezz2_norm(:,1)) 
    
    if ezz2_norm(ii,1)>=pp
        ezz2_point(pp+1,2) = ezz2_norm(ii,2);
        pp = pp + 1;
    end

end







%%
Array = []

kk = 0;

for ii = 1:10

    Array(ii,1:5) = es1_point(13-ii+kk:2:21-ii+kk,2);
    
    Array(ii,6:10) = gs2_point(13-ii+kk:2:21-ii+kk,2);

    Array(ii,11:15) = es2_point(13-ii+kk:2:21-ii+kk,2);

    Array(ii,16:20) = ezz_point(13-ii+kk:2:21-ii+kk,2);

end




for ii = 1:10

    Array(ii+10,1:5) = es12_point(13-ii+kk:2:21-ii+kk,2);
    
    Array(ii+10,6:10) = gs22_point(13-ii+kk:2:21-ii+kk,2);

    Array(ii+10,11:15) = es22_point(13-ii+kk:2:21-ii+kk,2);

    Array(ii+10,16:20) = ezz2_point(13-ii+kk:2:21-ii+kk,2);

end

map = inv(Array');


%%

map1d = (map(:,3)*3 + parameter(ii_t) * map(:,8) + parameter(jj_t) * map(:,13) + ...
                parameter2(kk_t) * map(:,1)+ parameter2(iii_t) * map(:,5))*14;

nq = 1000; 

xx = linspace(6,16,nq);

sum_es1 = zeros(1,nq);

sum_es2 = zeros(1,nq);

sum_gs1 = zeros(1,nq);

sum_gs2 = zeros(1,nq);

sum_es12 = zeros(1,nq);

sum_es22 = zeros(1,nq);

sum_gs12 = zeros(1,nq);

sum_gs22 = zeros(1,nq);

sum_ezz = zeros(1,nq);

sum_ezz2 = zeros(1,nq);

es1_inter = zeros(10,nq);

es2_inter = zeros(10,nq);

gs1_inter = zeros(10,nq);

gs2_inter = zeros(10,nq);

es12_inter = zeros(10,nq);

es22_inter = zeros(10,nq);

gs12_inter = zeros(10,nq);

gs22_inter = zeros(10,nq);

ezz_inter = zeros(10,nq);

ezz2_inter = zeros(10,nq);

for ii = 1:10

    es1_inter(ii,:) = interp1(es1_norm(:,1),es1_norm(:,2),linspace(6,16,nq)+5-ii);

    es2_inter(ii,:) = interp1(es2_norm(:,1),es2_norm(:,2),linspace(6,16,nq)+5-ii);

    gs1_inter(ii,:) = interp1(gs1_norm(:,1),gs1_norm(:,2),linspace(6,16,nq)+5-ii);

    gs2_inter(ii,:) = interp1(gs2_norm(:,1),gs2_norm(:,2),linspace(6,16,nq)+5-ii);

    es12_inter(ii,:) = interp1(es12_norm(:,1),es12_norm(:,2),linspace(6,16,nq)+5-ii);

    es22_inter(ii,:) = interp1(es22_norm(:,1),es22_norm(:,2),linspace(6,16,nq)+5-ii);

    gs12_inter(ii,:) = interp1(gs12_norm(:,1),gs12_norm(:,2),linspace(6,16,nq)+5-ii);

    gs22_inter(ii,:) = interp1(gs22_norm(:,1),gs22_norm(:,2),linspace(6,16,nq)+5-ii);

    ezz_inter(ii,:) = interp1(ezz_norm(:,1),ezz_norm(:,2),linspace(6,16,nq)+5-ii);

    ezz2_inter(ii,:) = interp1(ezz2_norm(:,1),ezz2_norm(:,2),linspace(6,16,nq)+5-ii);


    sum_es1 = sum_es1 + map1d(ii) * es1_inter(ii,:) + map1d(ii+10) * es12_inter(ii,:);

    sum_gs1 = sum_gs1 + map1d(ii) * gs1_inter(ii,:) + map1d(ii+10) * gs12_inter(ii,:);

    sum_es2 = sum_es2 + map1d(ii) * es2_inter(ii,:) + map1d(ii+10) * es22_inter(ii,:);

    sum_gs2 = sum_gs2 + map1d(ii) * gs2_inter(ii,:) + map1d(ii+10) * gs22_inter(ii,:);

    sum_ezz = sum_ezz + map1d(ii) * ezz_inter(ii,:) + map1d(ii+10) * ezz2_inter(ii,:);



end

Delta = sqrt(46^2 + 4*sum_es1.^2 + 4*sum_es2.^2)*1/2 - 46/2 - sqrt(255^2 + 4*sum_gs1.^2 + 4*sum_gs2.^2)*1/2 + 255/2 - sum_ezz;


figure(1)

plot((xx-11)/2,Delta,'linewidth',2)


grid on

xlabel('x/a')
ylabel('\Delta (GHz)')

set(gca,'Fontsize',15)

xticks([-2 -1 0 1 2])

set(gca,'ylim',[-10 70])

legend('Frequency shift')


%%


es1TwoD_norm(:,1) = es1TwoD(:,1)/244e-9;
es1TwoD_norm(:,2) = es1TwoD(:,2)/244e-9;
es1TwoD_norm(:,3) = es1TwoD(:,3)/1e9;

es2TwoD_norm(:,1) = es2TwoD(:,1)/244e-9;
es2TwoD_norm(:,2) = es2TwoD(:,2)/244e-9;
es2TwoD_norm(:,3) = es2TwoD(:,3)/1e9;

gs1TwoD_norm(:,1) = gs1TwoD(:,1)/244e-9;
gs1TwoD_norm(:,2) = gs1TwoD(:,2)/244e-9;
gs1TwoD_norm(:,3) = gs1TwoD(:,3)/1e9;

gs2TwoD_norm(:,1) = gs2TwoD(:,1)/244e-9;
gs2TwoD_norm(:,2) = gs2TwoD(:,2)/244e-9;
gs2TwoD_norm(:,3) = gs2TwoD(:,3)/1e9;

ezzTwoD_norm(:,1) = ezzTwoD(:,1)/244e-9;
ezzTwoD_norm(:,2) = ezzTwoD(:,2)/244e-9;
ezzTwoD_norm(:,3) = ezzTwoD(:,3)/1e9;



es1TwoD2_norm(:,1) = es1TwoD2(:,1)/244e-9;
es1TwoD2_norm(:,2) = es1TwoD2(:,2)/244e-9;
es1TwoD2_norm(:,3) = es1TwoD2(:,3)/1e9;

es2TwoD2_norm(:,1) = es2TwoD2(:,1)/244e-9;
es2TwoD2_norm(:,2) = es2TwoD2(:,2)/244e-9;
es2TwoD2_norm(:,3) = es2TwoD2(:,3)/1e9;

gs1TwoD2_norm(:,1) = gs1TwoD2(:,1)/244e-9;
gs1TwoD2_norm(:,2) = gs1TwoD2(:,2)/244e-9;
gs1TwoD2_norm(:,3) = gs1TwoD2(:,3)/1e9;

gs2TwoD2_norm(:,1) = gs2TwoD2(:,1)/244e-9;
gs2TwoD2_norm(:,2) = gs2TwoD2(:,2)/244e-9;
gs2TwoD2_norm(:,3) = gs2TwoD2(:,3)/1e9;

ezzTwoD2_norm(:,1) = ezzTwoD2(:,1)/244e-9;
ezzTwoD2_norm(:,2) = ezzTwoD2(:,2)/244e-9;
ezzTwoD2_norm(:,3) = ezzTwoD2(:,3)/1e9;

% [xq,yq] = meshgrid(-6:1e-2:6, -3:1e-2:3);
% vq = griddata(beta2D_norm(:,1),beta2D_norm(:,2),beta2D_norm(:,3),xq,yq);

sum1 = zeros(601,1001);
sum2 = zeros(601,1001);
sum3 = zeros(601,1001);
sum4 = zeros(601,1001);
sum5 = zeros(601,1001);



for ii = 1:10

    ii


    [xq,yq] = meshgrid((-5:1e-2:5)+5-ii, (-3:1e-2:3));

    vq1 = griddata(es1TwoD_norm(:,1),es1TwoD_norm(:,2),es1TwoD_norm(:,3),xq,yq);

    vq2 = griddata(es2TwoD_norm(:,1),es2TwoD_norm(:,2),es2TwoD_norm(:,3),xq,yq);

    vq3 = griddata(gs1TwoD_norm(:,1),gs1TwoD_norm(:,2),gs1TwoD_norm(:,3),xq,yq);

    vq4 = griddata(gs2TwoD_norm(:,1),gs2TwoD_norm(:,2),gs2TwoD_norm(:,3),xq,yq);

    vq5 = griddata(ezzTwoD_norm(:,1),ezzTwoD_norm(:,2),ezzTwoD_norm(:,3),xq,yq);

    vq6 = griddata(es1TwoD2_norm(:,1),es1TwoD2_norm(:,2),es1TwoD2_norm(:,3),xq,yq);

    vq7 = griddata(es2TwoD2_norm(:,1),es2TwoD2_norm(:,2),es2TwoD2_norm(:,3),xq,yq);

    vq8 = griddata(gs1TwoD2_norm(:,1),gs1TwoD2_norm(:,2),gs1TwoD2_norm(:,3),xq,yq);

    vq9 = griddata(gs2TwoD2_norm(:,1),gs2TwoD2_norm(:,2),gs2TwoD2_norm(:,3),xq,yq);

    vq10 = griddata(ezzTwoD2_norm(:,1),ezzTwoD2_norm(:,2),ezzTwoD2_norm(:,3),xq,yq);


    sum1 = sum1 + map1d(ii) * vq1 + map1d(ii+10) * vq6;
    sum2 = sum2 + map1d(ii) * vq2 + map1d(ii+10) * vq7;
    sum3 = sum3 + map1d(ii) * vq3 + map1d(ii+10) * vq8;
    sum4 = sum4 + map1d(ii) * vq4 + map1d(ii+10) * vq9;
    sum5 = sum5 + map1d(ii) * vq5 + map1d(ii+10) * vq10;


end

%%

[xq yq] = meshgrid((-5:1e-2:5)/2, (-0.84:1e-2:0.84));

Delta = sqrt(46^2 + 4*sum1.^2 + 4*sum2.^2)*1/2 - 46/2 - sqrt(255^2 + 4*sum3.^2 + 4*sum4.^2)*1/2 + 255/2 - sum5;

pcolor(xq,yq,abs(Delta(217:385,:)))

shading interp

colormap hot

colorbar

set(gca, 'ColorScale', 'log')

% caxis([-10 70])


set(gca,'xlim',[-2.5 2.5])

%%


figure(2)

semilogy(xx-11,abs(sum_gs1),'linewidth',2)

hold on

semilogy(xx-11,abs(sum_gs2),'linewidth',2)

%%

figure(10)

dcafteroptimize_norm(:,1) = dcafteroptimize(:,1)/244e-9;
dcafteroptimize_norm(:,2) = dcafteroptimize(:,2)/1e9;

plot(dcafteroptimize_norm(:,1)-10,dcafteroptimize_norm(:,2))

set(gca,'xlim',[-4.5 4.5])


%%

figure(1)

plot((gs1afterCE(1:3:2797,1)/244e-9-4)/2,gs1afterCE(1:3:2797,2)/1e9,'linewidth',2)

hold on

plot((es1afterCE(1:3:2797,1)/244e-9-4)/2,es1afterCE(1:3:2797,2)/1e9,'linewidth',2)

plot((gs2afterCE(1:3:2797,1)/244e-9-4)/2,gs2afterCE(1:3:2797,2)/1e9,'linewidth',2)




plot((es2afterCE(1:3:2797,1)/244e-9-4)/2,es2afterCE(1:3:2797,2)/1e9,'linewidth',2)


plot((gs1(1:3:2797,1)/244e-9-4)/2,gs1(1:3:2797,2)/1e9*(-1.7),'linewidth',2)


set(gca,'xlim',[0 5])
set(gca,'ylim',[-25 15])

xlabel('x/a')
ylabel('Frequency (GHz)')

legend('\Delta_{gs1}','\Delta_{es1}','\Delta_{gs2}','\Delta_{es2}','\Delta_{ZPL,zz}')

set(gca,'Fontsize',15)


grid on



figure(2)

Delta_gs = sqrt(255^2 + 4*gs1afterCE.^2 + 4*gs2afterCE.^2) - 255;

Delta_es = sqrt(46^2 + 4*es1afterCE.^2 + 4*es2afterCE.^2) - 46;

Delta_ZPL = -1.7*gs1;

Delta = -Delta_gs + Delta_es + Delta_ZPL;



plot((gs1afterCE(1:3:2797,1)/244e-9-4)/2,Delta(1:3:2797,2)/1e9,'linewidth',2)




xlabel('x/a')
ylabel('Frequency (GHz)')

legend('\Delta\nu')

set(gca,'Fontsize',15)

set(gca,'xlim',[0 5])
set(gca,'ylim',[10^-1.5 10^2.5])


grid on

% hold on
% 
% plot((es1afterCE(1:3:2797,1)/244e-9-4)/2,es1afterCE(1:3:2797,2)/1e9,'linewidth',2)
% 
% plot((gs2afterCE(1:3:2797,1)/244e-9-4)/2,gs2afterCE(1:3:2797,2)/1e9,'linewidth',2)
% 
% 
% 
% 
% plot((es2afterCE(1:3:2797,1)/244e-9-4)/2,es2afterCE(1:3:2797,2)/1e9,'linewidth',2)
% 
% 
% plot((gs1(1:3:2797,1)/244e-9-4)/2,gs1(1:3:2797,2)/1e9*(-1.7),'linewidth',2)

