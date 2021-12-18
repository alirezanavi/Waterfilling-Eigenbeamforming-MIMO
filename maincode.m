% code written by Alireza Navi
clc;
clear;
close all
Nr = [5 3]; % to obtain figure 3.4 and 3.5 of the book set Nr(1) = 4 and Nr(2) = 3
Nt = Nr(2:-1:1);
snrdb = -20:.1:20;
snr = 10.^(snrdb/10);
Niterr = 1000; % number of monte carlo iterations( since it runs the code 1000 times for faster result set the No. of iterations to 1 and for better results leave it be! :) )
c1 = zeros(length(snr),2);
c2 = zeros(length(snr),2);
c3 = zeros(length(snr),2);
c4 = zeros(length(snr),2);
mu = zeros(1,length(snr));
p_opt = zeros(length(snr),3);
for q =1:Niterr
h1 = 1/sqrt(2)*(randn(Nr(1),Nt(1))+1j*randn(Nr(1),Nt(1)));
lambda = sort(eig(h1*h1'),"descend")';
for k = 1:2
if k ~=1
    h1 = h1.';
end
%% equal power capacity(or ONLY CSIR)
for j = 1:length(snr)
c1(j,k) = c1(j,k) + sum(log2(1+1/Nt(k)*snr(j)*lambda(1:rank(h1))));
%% Eigenbeamforming
r = rank(h1);
    niter=0;
    while 1
    mu(j) = 1/(r-niter)*(1+1/snr(j)*sum(1./lambda(1:(r-niter))));
    p_opt(j,1:(r-niter)) = mu(j)- 1./(snr(j)*lambda(1:(r-niter)));
    if abs(p_opt(j,:))== p_opt(j,:)
    break
    else
    p_opt(j,end-niter:end) = 0;        
    niter = niter+1;
    end
    end
c2(j,k) = c2(j,k)+ sum(log2(1+snr(j)*p_opt(j,:).*lambda(1:rank(h1))));
%% Single Mode
c3(j,k) = c3(j,k)+ log2(1+snr(j)*1*lambda(1));
%% SISO
c4(j,k) = c4(j,k)+ log2(1+snr(j)*1);
end
end
end
for k=1:2
plot(snrdb,c1(:,k)/Niterr,LineWidth=.7,Color='k');
hold on;
grid on;
plot(snrdb,c2(:,k)/Niterr,'-.',LineWidth=.7,Color='k');
plot(snrdb,c3(:,k)/Niterr,'--',LineWidth=.7,Color='k');
plot(snrdb,c4(:,k)/Niterr,':',LineWidth=1.5,Color='k');
legend('Equal power','Eigenbeamforming','Single mode Beamforming','SISO');
xlabel('SNR(dB)');
ylabel('bits per channel use');
if k == 1
title('Nt=3 ,Nr = 5')
figure;
else
    title('Nt=5 ,Nr = 3')
end
end