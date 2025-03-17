
%% This file is to explore the influence of the varying parameters.

%% Default Setting Parameters

nb_reciever=5;                        % Number of receivers
ns = 256;                             % Number of points in RIR.
SNR = 20;                             % Signal Noise Ratio
Nx_set = ns+100;
r0 = 0.2;
center = [2 1.5 2];                   % Center of the reciever circle


%% Change of the Nb_Sensors
avg_nb = 100;
ns=256;
Nx =ns+100;
SNR = 20;
r0= 0.2;
max_n = 10;

NB_recievers = (2:1:10);
numNBr = length(NB_recievers);


nbr_lasso = zeros(numNBr,avg_nb);
nbr_tik = zeros(numNBr,avg_nb);
nbr_aot = zeros(numNBr,avg_nb);
nbr_bcot = zeros(numNBr,avg_nb);
nbr_lc2 = zeros(numNBr,avg_nb);

for seed =(1:avg_nb)
    for idx = 1:numNBr
        current_nbr = NB_recievers(idx);
        disp('Number of Sensor')
        disp(current_nbr);
        nmse_nbr = rir_estimation_func(seed,ns, Nx, current_nbr, SNR,r0);
        nbr_lasso(idx,seed) = nmse_nbr.lasso;
        nbr_tik(idx,seed) = nmse_nbr.tik;
        nbr_aot(idx,seed) = nmse_nbr.aot;
        nbr_bcot(idx,seed) = nmse_nbr.bcot;
        nbr_lc2(idx,seed)= nmse_nbr.L2;
    end
end

%% Plot NBr
e_lasso_nbr = 20*log10(mean(nbr_lasso,2));
e_tik_nbr = 20*log10(mean(nbr_tik,2));
e_aot_nbr = 20*log10(mean(nbr_aot,2));
e_bcot_nbr = 20*log10(mean(nbr_bcot,2));
e_lc2_nbr = 20*log10(mean(nbr_lc2,2));

figure;
plot(NB_recievers,e_lasso_nbr,'-s','LineWidth',1.5);
hold on
plot(NB_recievers,e_tik_nbr,'-o','LineWidth',1.5);
hold on
plot(NB_recievers,e_lc2_nbr,"-v",'LineWidth',1.5,'Color',"#77AC30");
hold on
plot(NB_recievers,e_aot_nbr,"-Hexagram",'LineWidth',1.5,'Color',"#EDB120");
hold on
plot(NB_recievers,e_bcot_nbr,"-pentagram",'LineWidth',1.5,'Color',"#7E2F8E");
hold off
grid("on");

legend({'Lasso','Tikhonov','$\ell_2$','Adjacent OMT','Barycenter OMT'},'Interpreter','latex','Fontsize',12);
ylabel('NMSE/dB','Interpreter','latex');
xlabel('Number of sensor','Interpreter','latex');
ylim([-5,-1.8]);
% xlim([-7,inf])


%% Change of the r0
avg_nb = 100;
ns=256;
SNR = 20;
nb_reciever = 5;
Nx =ns+100;
r0s= (0.05:0.05:0.5);
numr0 = length(r0s);

r0_lasso = zeros(numr0,avg_nb);
r0_tik = zeros(numr0,avg_nb);
r0_aot = zeros(numr0,avg_nb);
r0_bcot = zeros(numr0,avg_nb);
r0_lc2 = zeros(numr0,avg_nb);

for seed =(1:avg_nb)
    for idx = 1:numr0
        current_r0 = r0s(idx);
        disp('r0')
        disp(current_r0);
        nmse_r0 = rir_estimation_func(seed,ns, Nx, nb_reciever,SNR,current_r0);
        r0_lasso(idx,seed) = nmse_r0.lasso;
        r0_tik(idx,seed) = nmse_r0.tik;
        r0_aot(idx,seed) = nmse_r0.aot;
        r0_bcot(idx,seed) = nmse_r0.bcot;
        r0_lc2(idx,seed) = nmse_r0.L2;
    end
end


%% Plot r0
e_lasso_r0 = 20*log10(mean(r0_lasso,2));
e_tik_r0 = 20*log10(mean(r0_tik,2));
e_aot_r0 = 20*log10(mean(r0_aot,2));
e_bcot_r0 = 20*log10(mean(r0_bcot,2));
e_lc2_r0 = 20*log10(mean(r0_lc2,2));
r0s= (0.05:0.05:0.5);

figure;
plot(r0s,e_lasso_r0,'-s','LineWidth',1.5);
hold on
plot(r0s,e_tik_r0,'-o','LineWidth',1.5);
hold on
plot(r0s,e_lc2_r0,"-v",'LineWidth',1.5,'Color',"#77AC30");
hold on
plot(r0s,e_aot_r0 ,"-Hexagram",'LineWidth',1.5,'Color',"#EDB120");
hold on
plot(r0s,e_bcot_r0,"-pentagram",'LineWidth',1.5,'Color',"#7E2F8E");
hold off

grid("on");
legend({'Lasso','Tikhonov','$\ell_2$','Adjacent OMT','Barycenter OMT'},'Interpreter','latex','Fontsize',12);
ylabel('NMSE/dB','Interpreter','latex');
xlabel('$r_0/m$','Interpreter','latex');
% xlim([-inf,inf]);
% ylim([-5.5,inf]);


%% Changing the SNR
% close all
avg_nb = 100;
ns=256;
nb_reciever = 5;
Nx_set = ns+100;
r0= 0.2;
SNRs = 0:5:40;
numSNR = length(SNRs);

snr_lasso = zeros(numSNR,avg_nb);
snr_tik = zeros(numSNR,avg_nb);
snr_aot = zeros(numSNR,avg_nb);
snr_bcot = zeros(numSNR,avg_nb);
snr_lc2 = zeros(numSNR,avg_nb);


for seed =(11:20)
    for idx = 1:numSNR
        current_snr = SNRs(idx);
        disp('SNR')
        disp(current_snr);
        nmse_snr = rir_estimation_func(seed,ns, Nx_set, nb_reciever,current_snr,r0);
        snr_lasso(idx,seed) = nmse_snr.lasso;
        snr_tik(idx,seed) = nmse_snr.tik;
        snr_aot(idx,seed) = nmse_snr.aot;
        snr_bcot(idx,seed) = nmse_snr.bcot;
        snr_lc2(idx,seed) =  nmse_snr.L2;
    end
end

%% Plot 
SNRs = 0:5:40;
e_lasso_snr = 20*log10(mean(snr_lasso,2));
e_tik_snr = 20*log10(mean(snr_tik,2));
e_aot_snr = 20*log10(mean(snr_aot,2));
e_bcot_snr = 20*log10(mean(snr_bcot,2));
e_lc2_snr = 20*log10(mean(snr_lc2,2));

figure;
plot(SNRs,e_lasso_snr,'-s','LineWidth',1.5);
hold on
plot(SNRs,e_tik_snr,'-o','LineWidth',1.5);
hold on
plot(SNRs,e_lc2_snr,"-v",'LineWidth',1.5,'Color',"#77AC30");
hold on
plot(SNRs,e_aot_snr,"-Hexagram",'LineWidth',1.5,'Color',"#EDB120");
hold on
plot(SNRs,e_bcot_snr,"-pentagram",'LineWidth',1.5,'Color',"#7E2F8E");
hold off

grid("on");
ylabel('NMSE/dB','Interpreter','latex');
xlabel('SNR/dB','Interpreter','latex');
legend({'Lasso','Tikhonov','$\ell_2$','Adjacent OMT','Barycenter OMT',},'Fontsize',12,'Interpreter','latex');
% xlim([0,40]);
% ylim([-inf,inf]);


%% Change of Nx_set
avg_nb = 50;
ns=256;
SNR = 20;
nb_reciever = 5;
r0= 0.2;
Nx_sets =ns+(0:50:ns+150);
numNx = length(Nx_sets);

nx_lasso = zeros(numNx,avg_nb);
nx_tik = zeros(numNx,avg_nb);
nx_aot = zeros(numNx,avg_nb);
nx_bcot = zeros(numNx,avg_nb);
nx_lc2 = zeros(numNx,avg_nb);

for seed =(1:avg_nb)
    for idx = 1:numNx
        current_nx = Nx_sets(idx);
        disp('Nx')
        disp(current_nx);
        nmse_nx = rir_estimation_func(seed,ns, current_nx, nb_reciever,SNR,r0);
        nx_lasso(idx,seed) = nmse_nx.lasso;
        nx_tik(idx,seed) =  nmse_nx.tik;
        nx_aot(idx,seed) =  nmse_nx.aot;
        nx_bcot(idx,seed) = nmse_nx.bcot;
        nx_lc2(idx,seed) = nmse_nx.L2;
    end
end


%% Plot Nx

e_lasso_nx = 20*log10(mean(nx_lasso,2));
e_tik_nx = 20*log10(mean(nx_tik,2));
e_aot_nx = 20*log10(mean(nx_aot,2));
e_bcot_nx = 20*log10(mean(nx_bcot,2));
e_lc2_nx = 20*log10(mean(nx_lc2,2));

ratio_nx = (Nx_sets-ns)/ns;

figure;
plot(ratio_nx,e_lasso_nx,'-s','LineWidth',1.5);
hold on
plot(ratio_nx,e_tik_nx,'-o','LineWidth',1.5);
hold on
plot(ratio_nx,e_lc2_nx,"-v",'LineWidth',1.5,'Color',"#77AC30");
hold on
plot(ratio_nx,e_aot_nx,"-Hexagram",'LineWidth',1.5,'Color',"#EDB120");
hold on
plot(ratio_nx,e_bcot_nx,"-pentagram",'LineWidth',1.5,'Color',"#7E2F8E");
hold off
grid("on");

legend({'Lasso','Tikhonov','$\ell_2$','Adjacent OMT','Barycenter OMT',},'Fontsize',12,'Interpreter','latex');
xlabel('$(N_x-N_h)/N_h$','Interpreter','latex');
ylabel('NMSE/dB','Interpreter','latex');
% xlim([0.05,inf]);
% ylim([-8.5,0]);

