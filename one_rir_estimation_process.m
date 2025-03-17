clear all
close all
clc
%% This file is an entire process of the RIR estimation with all regularization methods.
% The scenarios has been set but can be customized.
% Current scenario is: SNR = 20dB, n_receiver = 5, Nx = ns + 100, r0=0.2m
% This file also plot the input signal waveform
% Room-Source-Array setup visualization
% One original RIR and estimated RIRs.
% 3D plot of all RIRs and estimated RIRs.
% NMSE are also calculated.
%% Function to calculate the NMSE
function NMSE = NMSE(H,H_hat,n)
    NMSE = 0;
    for i = 1:n
        nmse_i = sum_square(H(:,i)-H_hat(:,i))/sum_square(H(:,i));
        NMSE = NMSE + nmse_i;
    end
    NMSE = NMSE/n;
end

%% Experiment Parameters
% sigma2_t = 0 ;                      % Noise ratio, Not in Ehabets model
c = 343;                              % Speed of Sound
Room = [5 4 6];                       % Room Dimensions
source = [2 3.5 2];                   % Source Coordinates
beta = 0.5;                           % Reflection coefficient of walls.  
ns = 256;                             % number of samples
Nh = ns;
mtype ='omnidirectional';
SNR = 20;                             % Signal Noise Rati
order = 4;                            % Reflection order
%% Measurement Point Cloud (Without Error: if points larger than Room.)
center = [2 1.5 2];                   % Center of the reciever circle
r0=0.2;                               % Radius of the point cloud
n = 5;                                % Number of points in the cloud
seed = 1;                             % Determine the random point group 
meas_pts=meas_cloud_u(center,r0,seed,n,'circle');  % Coordinates list (n*3)

%% Visualize the Room, Source and Meas_Clds
Room_Visualization(Room,source,meas_pts)

%% Voice input signal 
[x, fs] = audioread('voicevibrato.wav');  % Voice Signal
x = x/abs(max(x));
downsample = 6;
x = decimate(x,downsample);
fs = fs/downsample;
start_t = 1.1;                            % Start Time 
Nx_set = Nh+100;  
N_start = ceil(start_t*fs);
N_end = N_start+Nx_set;
x = x(N_start:N_end);
Nx = length(x);                           % Identify the length of input signal
Xs = zeros(Nx-Nh,Nh);
for k = 1:Nx-Nh
    Xs(k,:) = flip(x(k+1:k+Nh));
end

%% Identify the properties of the input signal 

tx = (0:Nx-1)/fs;              % Identify the time axis correspoing to x
figure;
plot(tx, x)
title('Input Signal');
 
%% Calculate one reference RIR (By default, the time axes are always the same).
h_ref = rir_generator(c, fs, meas_pts(1,:), source, Room, beta, ns, mtype, order);  % Create a referenc RIR.
Nh = length(h_ref);                       % Length of IR
H = zeros(Nh,n);                          % Create an empty matrix to carry RIRs received by different microphones.
for i = 1:n
    h = rir_generator(c, fs, meas_pts(i,:), source, Room, beta, ns, mtype, order);   % One RIR. 
    H(:,i)=h(1:Nh)/sum(abs(h(1:Nh)));                                                % Fill to the matrix, with the normalization.
end



%% Calculate the Output signal. y=h*x+v 
Y = zeros(Nx,n);                          % Empty container of output signal y
for i = 1:n                               
    sigma = sqrt(sumsqr(filter(H(:,i),1,x))/(Nx)/power(10,SNR/10)); % Calculate the variance of the noise.
    a = sigma*randn(Nx,1);                % Noise
    Y(:,i)=filter(H(:,i),1,x)+a;          % yk=hk*x+vk 
end                                                                  
Y = Y(Nh+1:end,:);                        % Only take the latest Nx-Nh values to do the estimation.                                                           
t_est=(0:Nh-1)/fs;                        % The time axis of RIR.

%% Verification
cvx_solver mosek
cvx_begin quiet
    variable H_est(Nh,n)
    obj_est= 0;
    for i = 1:n
        obj_est= obj_est + sum_square(Y(:,i)-Xs*H_est(:,i));
    end
    minimize obj_est
cvx_end

disp('Opt Value:')
disp(obj_est)

%% Lasso Regularization  
N_lambda =10;                             % Number of tested regularization coefficients.                                       
lambdas = logspace(-6,0,N_lambda);        % Distribution of the regularization coefficients.
H_lassos = zeros(Nh,n,N_lambda);          % Empty matrix for all estimated RIRs.
NMSE_lasso = zeros(N_lambda,1);           % Empty vector of the NMSE values for different lambda.

for l = (1:N_lambda)  
    lambda = lambdas(l);
    cvx_begin quiet
    cvx_solver mosek
        variable H_lasso(Nh,n)            % Empty estimated matrix of all RIRs for one lambda.
        obj_lasso= 0;                     % Optimized value
        for i = 1:n                       
            obj_lasso= obj_lasso+sum_square(Y(:,i)-Xs*H_lasso(:,i));       % Term |y-Xh|.
        end
        obj_lasso = obj_lasso + lambda*norms(H_lasso(:),1);                % Regularization term.
        minimize(obj_lasso)                                                % Minimize part
    cvx_end

    H_lassos(:,:,l) = H_lasso;            % Record the estimated RIR for on lambda 
    nmse_lasso = NMSE(H,H_lasso,n);       % Calculated NMSE
    NMSE_lasso(l,1) = nmse_lasso;         % Record the NMSE value for the lambda
end

[min_lasso_nmse,I_lasso_nmse] = min(NMSE_lasso); % Find the optimum lambda who causes the lowest NMSE
H_lasso = H_lassos(:,:,I_lasso_nmse);            % Related estimated RIRs.
fprintf("Lambda within interval: %d \n " , I_lasso_nmse>1 && I_lasso_nmse<N_lambda);  % Represent if the optimum lambda is considered or out of expectation.
disp('Opt Lass NMSE Value:')
disp(min_lasso_nmse)

%% Tikhonov
N_lambda =10;
lambdas = logspace(-3,3,N_lambda);
H_tiks = zeros(Nh,n,N_lambda);
NMSE_tik = zeros(N_lambda,1);
for l = (1:N_lambda)
    lambda = lambdas(l);
    cvx_begin quiet
    cvx_solver mosek
        variable H_tik(Nh,n)
        obj_tik= 0;
        for i = 1:n
            obj_tik= obj_tik+sum_square(Y(:,i)-Xs*H_tik(:,i))+...
            lambda*sum_square(H_tik(:,i));
        end

        minimize(obj_tik)
    cvx_end

    H_tiks(:,:,l) = H_tik;
    nmse_tik = NMSE(H,H_tik,n);
    NMSE_tik(l,1) = nmse_tik;
end

[min_tik_nmse,I_tik_nmse] = min(NMSE_tik);
H_tik = H_tiks(:,:,I_tik_nmse);
fprintf("Lambda within interval: %d \n " , I_tik_nmse>1 && I_tik_nmse<N_lambda)
disp('Opt Tikhonov NMSE Value:')
disp(min_tik_nmse)

%% Adjacent OT
N_lambda =10;
[T1,T2]=meshgrid(t_est,t_est);
C=(T1-T2).^2*fs+power(10,-3);                  % Cost Function
lambdas = logspace(-3,3,N_lambda);
H_aots = zeros(Nh,n,N_lambda);
NMSE_aot = zeros(N_lambda,1);
for l = (1:N_lambda)
    lambda = lambdas(l);
    cvx_begin quiet
    cvx_solver mosek
        variable H_aot_p(Nh,n) nonnegative
        variable H_aot_n(Nh,n) nonnegative
        variable M_p(Nh, Nh,n-1) nonnegative
        variable M_n(Nh, Nh,n-1) nonnegative
        obj_aot= 0;
        for i = 1:n
            obj_aot= obj_aot+sum_square(Y(:,i)-Xs*(H_aot_p(:,i)-H_aot_n(:,i)));
        end
        for k = 2:n
            obj_aot = obj_aot+ lambda*sum(sum(C.*M_p(:,:,k-1)))...
                      + lambda*sum(sum(C.*M_n(:,:,k-1)));
            subject to 
                sum(M_p(:,:,k-1),2) == H_aot_p(:,k);
                sum(M_p(:,:,k-1),1)' == H_aot_p(:,k-1);
                sum(M_n(:,:,k-1),2) == H_aot_n(:,k);
                sum(M_n(:,:,k-1),1)' == H_aot_n(:,k-1);
        end
        minimize obj_aot 
        H_aot = H_aot_p-H_aot_n;
    cvx_end

    H_aots(:,:,l) = H_aot;
    nmse_aot = NMSE(H,H_aot,n);
    NMSE_aot(l,1) = nmse_aot;
end

[min_aot_nmse,I_aot_nmse] = min(NMSE_aot);
H_aot = H_aots(:,:,I_aot_nmse);
fprintf("Lambda within interval: %d \n " , I_aot_nmse>1 && I_aot_nmse<N_lambda)
disp('Opt AOT NMSE Value:')
disp(min_aot_nmse)


%% BaryCenter OT
N_lambda =10;
lambdas = logspace(-3,3,N_lambda);
H_bcots = zeros(Nh,n,N_lambda);
NMSE_bcot = zeros(N_lambda,1);
for l = (1:N_lambda)
    lambda = lambdas(l);
    cvx_begin quiet
    cvx_solver mosek
        variable H_bcot_p(Nh,n) nonnegative
        variable H_bcot_n(Nh,n) nonnegative
        variable M_bcot_p(Nh, Nh,n) nonnegative
        variable M_bcot_n(Nh, Nh,n) nonnegative
        variable h_bary_p(Nh,1) nonnegative
        variable h_bary_n(Nh,1) nonnegative
        obj_bcot= 0;
        for i = 1:n
            obj_bcot= obj_bcot+sum_square(Y(:,i)- ...
                        Xs*(H_bcot_p(:,i)-H_bcot_n(:,i)))+...
                        lambda*sum(sum(C.*M_bcot_p(:,:,i)))+...
                        lambda*sum(sum(C.*M_bcot_n(:,:,i)));
            subject to 
                sum(M_bcot_p(:,:,i),2) == H_bcot_p(:,i);
                sum(M_bcot_p(:,:,i),1)' == h_bary_p;
                sum(M_bcot_n(:,:,i),2) == H_bcot_n(:,i);
                sum(M_bcot_n(:,:,i),1)' == h_bary_n;
        end
        minimize(obj_bcot)
        H_bcot = H_bcot_p-H_bcot_n;
    cvx_end

    H_bcots(:,:,l) = H_bcot;
    nmse_bcot = NMSE(H,H_bcot,n);
    NMSE_bcot(l,1) = nmse_bcot;
end

[min_bcot_nmse,I_bcot_nmse] = min(NMSE_bcot);
H_bcot = H_bcots(:,:,I_bcot_nmse);
fprintf("Lambda within interval: %d \n " , I_bcot_nmse>1 && I_bcot_nmse<N_lambda)
disp('Opt BaryCenter OT NMSE Value:')
disp(min_bcot_nmse)

%% L2 (2 regularization coefficients, double coefficients-searching loop)

N_lambda=10;
lambdas = logspace(-3,3,N_lambda);
H_lc2s= zeros(Nh,n,N_lambda);
NMSE_lc2 = zeros(N_lambda,1);

N_mu = 2;
mus = logspace(-3,3,N_mu);
H_lc2s_m= zeros(Nh,n,N_mu);
NMSE_lc2_m = zeros(N_mu,1);

for l = (1:N_lambda)
    lambda = lambdas(l);
    for m =(1:N_mu)
        mu = mus(m);
        cvx_begin  quiet
            variable H_L2(Nh, n)
            variable H_l2(Nh,1)
            obj_l2 = 0;
            re = 0;
            for k = 1:n
                obj_l2 = obj_l2 + sum_square(Y(:,k) - Xs*H_L2(:,k)); 
            end
            for k = 1:n
                re = re+ lambda*(sum_square(H_L2(:, k) - H_l2)+mu*sum_square(H_L2(:, k)));
            end
    
            subject to
                H_l2 == sum(H_L2,2)/n;
            minimize(obj_l2 + re)
        cvx_end
        
        H_lc2s_m(:,:,m) = H_L2;
        nmse_lc2_m = NMSE(H,H_L2,n);
        NMSE_lc2_m(m,1) = nmse_lc2_m;
        
    end
    [~,I_Lc2_nmse_m] = min(NMSE_lc2_m);

    H_lc2s(:,:,l) = H_lc2s_m(:,:,I_Lc2_nmse_m);
    nmse_lc2 =  NMSE_lc2_m(I_Lc2_nmse_m,1);
    NMSE_lc2(l,1) = nmse_lc2;
end

[min_Lc2_nmse,I_Lc2_nmse] = min(NMSE_lc2);
H_L2 = H_lc2s(:,:,I_Lc2_nmse);
fprintf("Lambda within interval: %d \n " , I_Lc2_nmse>1 && I_Lc2_nmse<N_lambda)
disp('Opt L2 NMSE Value:')
disp(min_Lc2_nmse)



%% Representation of results of the example results.
idx_ex = 2;
figure();
plot(H(:,idx_ex),linewidth=1);
hold on
plot(H_lasso(:,idx_ex),linewidth=1);
hold on
plot(H_aot(:,idx_ex),linewidth=1);
hold on
plot(H_bcot(:,idx_ex),'--',linewidth=1);
hold on
plot(H_tik(:,idx_ex),linewidth=1,color='green');
hold off

legend('Standard','Lasso','AOT','BCOT','tik')



%% Plot and Compare the Results.
low_cut = 5e-3;
hlog = log10(max(abs(H),low_cut));
hlog_lasso = log10(max(abs(H_lasso),low_cut));
hlog_tik = log10(max(abs(H_tik),low_cut));
hlog_aot = log10(max(abs(H_aot),low_cut));
hlog_bcot = log10(max(abs(H_bcot),low_cut));
hlog_lc2 = log10(max(abs(H_L2),low_cut));


h_max = max([max(hlog(:)),max(hlog_lasso(:)),max(hlog_tik(:)),max(hlog_aot(:)),max(hlog_bcot(:)),max(hlog_lc2(:))]);
h_min = log10(low_cut);

[Tap,Nb_meas] = meshgrid(t_est,(1:n));
figure();
subplot(6,1,1)
surf(Tap,Nb_meas,hlog','EdgeColor', 'none');
title('Standard IRs');

subplot(6,1,2)
surf(Tap,Nb_meas,hlog_lasso','EdgeColor', 'none');
title('Lasso');

subplot(6,1,3)
surf(Tap,Nb_meas,hlog_tik','EdgeColor', 'none');
title('Tikhonov');

subplot(6,1,4)
surf(Tap,Nb_meas,hlog_aot','EdgeColor', 'none');
title('OT');

subplot(6,1,5)
surf(Tap,Nb_meas,hlog_bcot','EdgeColor', 'none');
title('BCOT');

subplot(6,1,6)
surf(Tap,Nb_meas,hlog_lc2','EdgeColor', 'none');
title('L2');

for j=(1:6)
    subplot(6,1,j)
    view(2); % View as 2D
    ylabel('Idx of IR');
    xlabel('Time/s');
    colormap([ parula(256)]);    % Add white for zero-mass areas 
    clim([h_min, h_max]);
    colorbar;                            % Add a colorbar
end




%% Evaluation with the NMSE signal
L = 20;
nmse_lasso = NMSE(H,H_lasso,n);
nmse_tik = NMSE(H,H_tik,n);
nmse_aot = NMSE(H,H_aot,n);
nmse_bcot = NMSE(H,H_bcot,n);
nmse_lc2 = NMSE(H,H_L2,n);

