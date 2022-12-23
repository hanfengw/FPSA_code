%%
% Hanfeng Wang 12.23.2022
% For "Field Programmable Spin Arrays for Scalable Quantum Repeaters" paper
% Superradiance

%% Fig.S4e

figure(1)
Gammasp = 100e6;

t = linspace(0,80,100)*1e-9

semilogy(t,Gammasp*exp(-Gammasp*t)/1e9)

hold on

semilogy(t,2*Gammasp*exp(-2*Gammasp*t)/1e9)

set(gca,'Fontsize',15)

set(gca,'ylim',[1e-5 1e0])

figure(2)

F = exp(Gammasp*t)./(2+exp(Gammasp*t))

semilogy(t,1-F)

%% Fig.S4f

pdet = logspace(-3,0,100);

R1 = 2*0.01*pdet;

T1 = 40e-6 + 5.5e-6./R1;

R2 = pdet.^2/2;

T2 = 40e-6 + 11e-6./R2;

loglog(pdet,1./T1)

hold on

loglog(pdet,1./T2)

set(gca,'ylim',[1e-2 2e4])


F = 0.99;

t = log(1.98/0.01)/100e6;

p = exp(-100e6*t);

R3 = pdet*p/2;

T3 = 40e-6 + 5.5e-6./R3;

loglog(pdet,1./T3)


%%

clear

pdetlist = logspace(-3,0,20);
p = 0.005;

for iii = 1:length(pdetlist)

    N = 0;
    T = 0;
    pdet = pdetlist(iii);

    for ii = 1:1000000

        index = 0;

        R = rand(1);

        if R<p*pdet/2
            index = 1;
        end

        if R>=pdet/2
            index = 0;
        end

        if R>=p*pdet/2 && R<pdet/2

            R2 = rand(1);

            if R2<=pdet
                index = 2;
            else
                index = 3;
            end

        end

        if index == 1
            T = T+5.5e-6;
            N = N+1;
        end

        if index == 3
            T = T+11e-6;
        end

        if index == 0
            T = T+5.5e-6;
        end

        if index == 2
            T = T+11e-6;
            N = N+1;
        end

    end

    pf(iii) = N./(T+40e-6*N)

end

loglog(pdetlist,pf)

hold on

