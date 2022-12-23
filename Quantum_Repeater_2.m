%%
% Hanfeng Wang 12.23.2022
% For "Field Programmable Spin Arrays for Scalable Quantum Repeaters" paper
% Added rate per added qubit

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


%% mzi

distance = 1;



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




%% 

rateenvelop = max(rate);

figure(2)

for ii = 2:length(rateenvelop)-1
    
    slope(ii) = (rate(1,ii+1)-rate(1,ii-1))/(N(ii+1)-N(ii-1));

    slope2(ii) = (rate2(ii+1)-rate2(ii-1))/(N(ii+1)-N(ii-1));

    slope3(ii) = (rateenvelop(ii+1)-rateenvelop(ii-1))/(N(ii+1)-N(ii-1));


end

plot(N(2:end-1),slope(2:end),'linewidth',2)

hold on

plot(N(2:end-1),slope2(2:end),'linewidth',2)



plot(N(2:end-1),slope3(2:end),'.','MarkerSize',10)



grid on

set(gca,'Fontsize',15)

xlabel('Number of qubits')
ylabel('\Delta\Gamma/qubit (ebit/s)')

legend('FPSA','MZI tree','Hybrid structure')


set(gca,'xlim',[100 4000])

xticks([1000 2000 3000 4000])
