clear; clc; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is for "parfor" and CoreNum = 1 is like using simple "for"
CoreNum = 12;
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(CoreNum);
else
    disp('matlab pool already started');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = rng(100);       % Seed the RNG for repeatability

% You may use a finite #cycles in practice or large N 
maxcycles = 1e15;

EbNo_vec = 1:0.5:3.5; % we have loop for different SNR values
K = 64;
N = 128;
M = 16;  % M can be selected from [1,2,4, 8, 16, 32,64] for testing.
R = K/N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fer = zeros(1,length(EbNo_vec)); % Frame error rate for each SNR
FE = zeros(1,length(EbNo_vec));  % #of frames in errors
Nruns = zeros(1,length(EbNo_vec)); % #of actual runs
total_num_of_ps_update = zeros(1, length(EbNo_vec));
maxRun = 1e9;

ClustNum = [5000 5000 5000 5e3*(ones(1,length(EbNo_vec)-3))]; % adjust this based on N and K
maxFE = 100; % maximum number of frame errors, you may use 300 for smoother curves

fprintf('-------------------------------------\n');
for EbNo_count = 1:length(EbNo_vec)
    tic;
    EbNo_dB = EbNo_vec(EbNo_count);
    EbNo = 10^(EbNo_dB/10);
    sigma = 1/sqrt(2*R*EbNo);
    Nblkerrs = 0;
    tot_num_of_ps_upd = 0;

    fprintf('[%02d:%02d] Starting! SNR = %.2f\n',0, 0, EbNo_dB);
    for i = 1:maxRun/ClustNum(EbNo_count)
        
        parfor j = 1:ClustNum(EbNo_count)
            %Generate random message
            msg = randi([0 1], 1, K);
            % encoding
            obj = FPSU_SC(K, N, randi([0 1], 1, K));
            [U, x] = obj.polar_enc();
            % BPSK modulation
            modulated = 1 - 2*x;
            % AWGN
            r = modulated + randn(1,N)*sigma;
            % Decoding
            [U_hat, X_hat, num] = obj.polar_SCD(r, M);
            tot_num_of_ps_upd = tot_num_of_ps_upd + num;
            Nblkerrs = Nblkerrs + any(X_hat~=x);
        end
        if Nblkerrs >= maxFE
            break;
        end
        t = toc;
        
        elapsed_m = t/60;
        elapsed_h = floor(elapsed_m/60);
        elapsed_m = floor(elapsed_m - elapsed_h*60);
        fprintf(2,'[%02d:%02d] EbNo = %.2f, Frame = %d, FE = %d\n',elapsed_h,elapsed_m,EbNo_dB,i*ClustNum(EbNo_count),Nblkerrs);
    end
    t = toc;
    elapsed_m = t/60;
    elapsed_h = floor(elapsed_m/60);
    elapsed_m = floor(elapsed_m - elapsed_h*60);
    fprintf(2,'[%02d:%02d] SNR = %.1f, Frame = %d, FE = %d\n',elapsed_h,elapsed_m,EbNo_dB,i*ClustNum(EbNo_count),Nblkerrs);
    
    temp = (i*ClustNum(EbNo_count));
    Nruns(EbNo_count) = temp;
    Fer(EbNo_count) = Nblkerrs/temp;
    FE(EbNo_count) = Nblkerrs;
    total_num_of_ps_update(EbNo_count) = tot_num_of_ps_upd;
    fprintf('-------------------------------------\n');
end

rng(s);     % Restore RNG

% Creating sim_state
res.trials = Nruns;
res.frame_errors = FE;
res.FER = FE./Nruns;
res.average_num_of_ps_update = total_num_of_ps_update ./ Nruns;

res.SNR = EbNo_vec;
res.K = K;
res.delta = delta;

res.max_cycles_allowed = maxcycles;
res.max_trials = maxRun;
res.max_frame_errors = maxFE;

filename = sprintf('Polar_FPSU_K%d_N%d_SNR%0.1f-%0.1f_M%d.mat',K,N,EbNo_vec(1),EbNo_vec(end), M);
save(filename, 'res');

figure;
semilogy(res.SNR, res.FER, '-s', 'LineWidth', 1.5);
grid on;
xlabel('$E_b/N_0$ (dB)', 'Interpreter', 'latex');
ylabel('BLER', 'Interpreter', 'latex');
