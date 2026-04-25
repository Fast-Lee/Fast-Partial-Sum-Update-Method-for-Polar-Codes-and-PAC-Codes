%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2026.3, Xun Li
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that:
% the source code retains the above copyright notice, 
% and te redistribtuion condition.
% Freely distributed for educational and research purposes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
EbNo_vec = 1:0.5:4; % we have loop for different SNR values
K = 64;
N = 128;
R = K/N;
poly = 3211;
polyb = dec2bin(base2dec(num2str(poly), 8)) - '0';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fer = zeros(1,length(EbNo_vec)); % Frame error rate for each SNR
FE = zeros(1,length(EbNo_vec));  % #of frames in errors
Nruns = zeros(1,length(EbNo_vec)); % #of actual runs
maxRun = 1e9;

ClustNum = [500 500 1e3 5e3*(ones(1,length(EbNo_vec)-3))]; % adjust this based on N and K
maxFE = 500; % maximum number of frame errors, you may use 300 for smoother curves

fprintf('-------------------------------------\n');
for EbNo_count = 1:length(EbNo_vec)
    tic;
    EbNo_dB = EbNo_vec(EbNo_count);
    EbNo = 10^(EbNo_dB/10);
    sigma = 1/sqrt(2*R*EbNo);
    Nblkerrs = 0;
    
    fprintf('[%02d:%02d] Starting! SNR = %.2f\n',0, 0, EbNo_dB);
    for i = 1:maxRun/ClustNum(EbNo_count)
        
        parfor j = 1:ClustNum(EbNo_count)
            %Generate random message
            msg = randi([0 1], 1, K);
            % Polar transformation
            pac = FPSU_PAC(N, K, polyb);
            x = pac.encode(msg)';
            % BPSK modulation
            modulated = 1 - 2*x;
            % AWGN
            r = modulated + randn(1, N)*sigma;
            % Decoding
            % d_esti = pac.SC_decoder(r')';
            M = 1;  % select M = [1,2 4,8,16, 32] for test.
            [d_esti] = pac.SC_decoder(r, M);
            Nblkerrs = Nblkerrs + any(d_esti~=msg);
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
    fprintf('-------------------------------------\n');
end

rng(s);     % Restore RNG
% Creating sim_state
res.trials = Nruns;
res.frame_errors = FE;
res.FER = FE./Nruns;

res.SNR = EbNo_vec;
res.K = K;

res.max_cycles_allowed = maxcycles;
res.max_trials = maxRun;
res.max_frame_errors = maxFE;

filename = sprintf('PAC_FPSU_K%d_N%d_SNR%0.1f-%0.1f.mat',K,N,EbNo_vec(1),EbNo_vec(end));
save(filename, 'res');

figure;
semilogy(res.SNR, res.FER, '-s', 'LineWidth', 1.5);
grid on;
xlabel('$E_b/N_0$ (dB)', 'Interpreter', 'latex');
ylabel('BLER', 'Interpreter', 'latex');
