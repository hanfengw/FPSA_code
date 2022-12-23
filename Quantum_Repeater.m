% Hanfeng Wang 12.23.2022
% For "Field Programmable Spin Arrays for Scalable Quantum Repeaters" paper

% This code is for the FPSA and MZI rate calculation. The detailed method
% is mentioned in the supplementary material of this paper. 

%% fpsa
clear

distance = 1;

N = 24:24:4000;

eta = exp(-0.041*distance);

loss = 0.0008; % 4e-3 dB



for ii = 1:length(N)

    p = zeros(1,N(ii)/2);

    for jj = 1:N(ii)/2
    
        p(jj) = 2*0.01*eta*0.33*0.83*0.25*((1-loss)^jj)/2;
        
    end

    A = zeros(length(p),length(p));

    p = 1 - (1-p).^250;

    q = p./(1-p);

    for kk = 1:length(p)
        
        A(kk,1) = sum(q(1:kk));

    end

    for i = 2:length(p)
        
        for j = 2:length(p)
            
            A(i,j) = A(i-1,j) + q(i)*A(i-1,j-1);

        end

    end

    norm = 1;

    for i = 1:length(p)

        norm = norm * (1-p(i));

    end

    A = A*norm;

    p2 = 0;
    
    for i = 1:length(p)

        p2 = p2 + 2*A(end,i)*i*sum(A(end,i+1:end)) + i*A(end,i)^2;
       %p2 = p2 + A(end,i)*i;

    end

    NN(ii) = p2;


end




for ii = 1:length(N)

    p1(ii) = (0.83*0.25*(1-loss)^(N(ii)/2))^2/2;

    T(ii) = 40e-6 + 11e-6/p1(ii) + 5.5e-6*250 + 1e3/3e8*250;

end

rate = NN./T;

plot(N,rate,'linewidth',2)

hold on

%% mzi

distance = 1;

N = 24:24:4000;

eta = exp(-0.041*distance);



for ii = 1:length(N)

    
    pmzi = 2*0.01*eta*0.33*0.83*0.25*0.912^log2(N(ii))/2;

    pmzi = 1 - (1-pmzi).^250;
        
    p2 = N(ii)/2*pmzi;

    NN(ii) = p2;


end



for ii = 1:length(N)

    p1(ii) = (0.83*0.25*0.912^log2(N(ii)))^2/2;

    T(ii) = 40e-6 + 11e-6/p1(ii) + 5.5e-6*250 + 1e3/3e8*250;

end

rate2 = NN./T;

plot(N,rate2,'linewidth',2)

hold on

set(gca,'ylim',[0 2.5e4])

set(gca,'FontSize',15)

legend('FPSA','MZI tree')

xlabel('Number of qubits')
ylabel('Entanglement rate (ebits/s)')


%% N devices

clear

distance = 1;

N = 24:24:4400;

% N = 4:5000:50000;

eta = exp(-0.041*distance);

loss = 0.0008; % 4e-3 dB

ndevlist = 1:1:10;



for kkk = 1:length(ndevlist)

    kkk

    ndev = ndevlist(kkk);

    N_per = floor(N/ndev);

for ii = 1:length(N)


    p = zeros(1,floor(N_per(ii)/2));

    for jj = 1:floor(N_per(ii)/2)
    
        p(jj) = 2*0.01*eta*0.33*0.83*0.25*((1-loss)^jj)/2*0.912^log2(ndev);
        
    end

    pall = repmat(p,1,ndev);

    pall = 1 - (1-pall).^250;

    norm1 = 1;

    for i = 1:floor(length(pall)/8*3)

        norm1 = norm1 * (1-pall(i));

    end

    norm2 = 1;

    for i = floor(length(pall)/8*3) + 1:length(pall)

        norm2 = norm2 * (1-pall(i));

    end



    A = zeros(length(pall),length(pall));

    q = pall./(1-pall);

    for kk = 1:length(pall)
        
        A(kk,1) = sum(q(1:kk))*norm1;

    end

    for i = 2:length(pall)
        
        for j = 2:length(pall)
            
            A(i,j) = A(i-1,j) + q(i)*A(i-1,j-1);

        end

    end



    A = A*norm2;

    p2 = 0;
    
    for i = 1:length(pall)

       p2 = p2 + 2*A(end,i)*i*sum(A(end,i+1:end)) + i*A(end,i)^2;


    end

    NN(ii) = p2;

   


end




for ii = 1:length(N)

    p1(ii) = (0.83*0.25*(1-loss)^(N(ii)/2/ndev)*0.912^log2(ndev))^2/2;

    T(ii) = 40e-6 + 11e-6/p1(ii) + 5.5e-6*250 + 1e3/3e8*250;

end

rate(kkk,:) = NN./T;

end

figure(1)

plot(N,rate(1,:),'linewidth',2)

hold on

plot(N,rate(2,:),'linewidth',2)

% plot(N,rate(3,:),'linewidth',2)

plot(N,max(rate),'--','linewidth',2)

set(gca,'Xlim',[0 4000])

set(gca,'Ylim',[0 12e4])

set(gca,'Fontsize',15)

yticks([0 4e4 8e4 12e4])

legend('N_{dev} = 1','N_{dev} = 2','Envelope')



%% n devices

clear

distancelist = 1:1:20;

% distancelist = 10;

N = logspace(2,4.8,20);




loss = 0.0008; % 4e-3 dB

ndevlist = 1:1:40;



% ndevlist = [40];

for xx = 1:length(distancelist)

    xx

    distance = distancelist(xx);

    eta = exp(-0.041.*distance);

    for kkk = 1:length(ndevlist)


        ndev = ndevlist(kkk);

        N_per = floor(N/ndev);

        for ii = 1:length(N)

            p = zeros(1,floor(N_per(ii)/2));

            for jj = 1:floor(N_per(ii)/2)

                p(jj) = 2*0.01*eta*0.33*0.83*0.25*((1-loss)^jj)/2*0.912^log2(ndev);

            end

            pall = repmat(p,1,ndev);

            pall = 1 - (1-pall).^250;

            norm = ones(1,4);

            i0 = 5e10;

            i1 = 5e10;

            i2 = 5e10;

            i3 = 5e10;

            for i0 = 1:length(pall)

                

                if norm(1) < 1e-300

                    break

                end

                norm(1) = norm(1) * (1-pall(i0));

            end

            if i0 < length(pall)

                for i1 = i0:length(pall)

                    

                    if norm(2) < 1e-300

                        break

                    end

                    norm(2) = norm(2) * (1-pall(i1));

                end

            end

            if i1 < length(pall)

                for i2 = i1:length(pall)

                    

                    if norm(3) < 1e-300

                        break

                    end

                    norm(3) = norm(3) * (1-pall(i2));

                end

            end

            if i2 < length(pall)

                for i3 = i2:length(pall)

                    norm(4) = norm(4) * (1-pall(i3));

                end

            end

            A = zeros(length(pall),length(pall));

            q = pall./(1-pall);

            for kk = 1:length(pall)

                A(kk,1) = sum(q(1:kk))*norm(1);

            end

            normi = 2;

            for i = 2:length(pall)

                for j = 2:length(pall)

                    A(i,j) = A(i-1,j) + q(i)*A(i-1,j-1);

                    if A(i,j) > 1e300
                        
                        A = A*norm(normi);

                        normi = normi+1;
                        
                    end

                end

            end

            if sum(A(end,:))>1
                
                A = A*norm(normi);

            end
% 
%             A = A*norm2;

            p2 = 0;

            for i = 1:length(pall)

                % p2 = p2 + 2*A(end,i)*i*sum(A(end,i:end));
                p2 = p2 + 2*A(end,i)*i*sum(A(end,i+1:end)) + i*A(end,i)^2;

            end

            NN(ii) = p2;

        end




        for ii = 1:length(N)

            p1(ii) = (0.83*0.25*(1-loss)^(N(ii)/2/ndev)*0.912^log2(ndev))^2/2;

            T(ii) = 40e-6 + 11e-6/p1(ii) + 5.5e-6*250 + 1e3/3e8*250;

        end

        rate(kkk,:) = NN./T;



    end

    Finalrate(xx,:) = max(rate);


end

%%
h = pcolor(round(N),distancelist,Finalrate);

hold on
shading interp
set(h, 'EdgeColor', 'none');

set(gca, 'XScale', 'log')

set(gca, 'Xlim', [1e2 1e5])

set(gca, 'ColorScale', 'log')

colormap hot

colorbar

caxis([0.9e3 3e5])


distancelist1 = [1:0.1:40];

area(3000*distancelist1,distancelist1,'LineStyle','--','LineWidth',1.1)

colororder([0.7 0.7 0.7]*8.5/7)

% set(gca, 'YScale', 'log')


% plot(N,max(rate),'--','linewidth',2)
% 
% set(gca,'Xlim',[0 4000])
% 
% set(gca,'Ylim',[0 12e4])
% 
% set(gca,'Fontsize',15)
% 
% yticks([0 4e4 8e4 12e4])
% 
% legend('N_{dev} = 1','N_{dev} = 2','Envelope')





