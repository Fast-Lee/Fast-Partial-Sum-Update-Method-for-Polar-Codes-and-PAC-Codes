%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2026.3, XXXXXX & XXXXXX
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that:
% the source code retains the above copyright notice, 
% and te redistribtuion condition.
% Freely distributed for educational and research purposes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef FPSU_SCL

    properties
        N         
        K  
        rate               
        rate_profiling      
        QW_i
        B
        GN                  % polar generator matrix
    end
    
    methods
        function obj = FPSU_SCL(N, K, QW_i)
            n = ceil(log2(N));
            N = 2^n;
            obj.N = N;
            obj.K = K;
            obj.rate = K / N;
            obj.QW_i = QW_i;
            obj.rate_profiling = PW_rate_profiling(obj);
            obj.GN = FPSU_SCL.get_GN(obj.N);                 
        end
        
        % polarization-weight rate profile
        function info_indices = PW_rate_profiling(obj)
            valid = obj.QW_i(obj.QW_i < obj.N);
            info_indices = sort(valid(end-obj.K+1:end) + 1);
        end
        
        % polar encoder
        function x = encode(obj, d)
            if (length(d) ~= obj.K)
                error('The length of the input d is not equal to k.')
            end
            u = zeros(1, obj.N);
            u(obj.rate_profiling) = d;
            % Polar Encoding
            x = mod(u*obj.GN, 2)';
        end

        % FPSU polar SCL decoding
        function [d_esti, PM_best, U_list, total_number_of_ps_update] = FPSU_Polar_SCL(obj, llr, L, M)
            total_number_of_ps_update = 0;
            B_RAM = zeros(1, obj.N);
            B_RAM(obj.rate_profiling) = 1;
            % Same LLR RAM depth rule as the verified FPSU_SC
            DEPTH = FPSU_SCL.get_llr_ram_depth(obj.N, M);
            % Path states
            LLR_Ram = zeros(DEPTH, 2*M, L);
            PS_RAM  = zeros(obj.N/M/2, M, L);
            U_path  = zeros(obj.K, L);
            % Only path 1 is active at the beginning
            PM = inf(1, L);
            PM(1) = 0;
            active_path = false(1, L);
            active_path(1) = true;
            % Initialize channel llr to path 1 with the same butterfly mapping
            for i = 0:obj.N-1
                if i < obj.N/2
                    addr = i * 2;
                else
                    addr = (i - obj.N/2) * 2 + 1;
                end
                LLR_Ram(floor(addr/(2*M)) + 1, mod(addr, 2*M) + 1, 1) = llr(i + 1);
            end

            S_MAX = log2(obj.N);
            STAGE = log2(M) + 1;
            info_cnt = 1;
            G = [1 zeros(1, M - 1)];

            for i = 0:obj.N-1
                %==========================================================
                % 1) LLR update for every active path
                %==========================================================
                for l_index = 1:L
                    if ~active_path(l_index)
                        continue;
                    end
                    LLR_Ram(:,:,l_index) = obj.update_llr_path(LLR_Ram(:,:,l_index), PS_RAM(:,:,l_index),i, S_MAX, STAGE, DEPTH, M);
                end

                %==========================================================
                % 2) Bit decision / path split / path pruning
                %==========================================================
                chosen_bit = zeros(1, L);
                if B_RAM(i + 1) == 0  % Frozen bit, all active paths take 0
                    for l_index = 1:L
                        if ~active_path(l_index)
                            continue;
                        end
                        root_llr = LLR_Ram(DEPTH, 1, l_index);
                        PM(l_index) = FPSU_SCL.update_pm(PM(l_index), root_llr, 0);
                        chosen_bit(l_index) = 0;
                    end
                else   % Information bit, split every active path
                    active_idx = find(active_path);
                    num_cand = 2 * numel(active_idx);
                    cand_PM  = inf(1, num_cand); 
                    cand_src = zeros(1, num_cand);
                    cand_bit = zeros(1, num_cand);
                    c = 1;
                    for k = 1:numel(active_idx)
                        src = active_idx(k);
                        root_llr = LLR_Ram(DEPTH, 1, src);
                        cand_PM(c)   = FPSU_SCL.update_pm(PM(src), root_llr, 0);
                        cand_src(c)  = src;
                        cand_bit(c)  = 0;
                        c = c + 1;
                        cand_PM(c)   = FPSU_SCL.update_pm(PM(src), root_llr, 1);
                        cand_src(c)  = src;
                        cand_bit(c)  = 1;
                        c = c + 1;
                    end

                    [~, order] = sort(cand_PM, 'ascend');
                    keep_num = min(L, num_cand);

                    old_LLR_Ram = LLR_Ram;
                    old_PS_RAM  = PS_RAM;
                    old_U_path  = U_path;
                    LLR_Ram = zeros(size(old_LLR_Ram));
                    PS_RAM  = zeros(size(old_PS_RAM));
                    U_path  = zeros(size(old_U_path));
                    PM = inf(1, L);
                    active_path = false(1, L);
                    chosen_bit = zeros(1, L);
                    for l_index = 1:keep_num
                        idx = order(l_index);
                        src = cand_src(idx);
                        active_path(l_index) = true;
                        PM(l_index) = cand_PM(idx);
                        chosen_bit(l_index) = cand_bit(idx);
                        LLR_Ram(:, :, l_index) = old_LLR_Ram(:, :, src);
                        PS_RAM(:, :, l_index)  = old_PS_RAM(:, :, src);
                        U_path(:, l_index)     = old_U_path(:, src);
                        U_path(info_cnt, l_index) = cand_bit(idx);
                    end
                    info_cnt = info_cnt + 1;
                end

                %==========================================================
                % 3) Starting tree level for decoding the next bit
                %==========================================================
                S_MAX = sum(bitget(bitxor(i, i+1), 1:32));

                %==========================================================
                % 4) Partial-sum update for every active path 
                %==========================================================
                for l_index = 1:L
                    if ~active_path(l_index)
                        continue;
                    end
                    [PS_RAM(:,:,l_index), num] = obj.FPSU_path(PS_RAM(:,:,l_index), chosen_bit(l_index), G, i, M);
                    total_number_of_ps_update = total_number_of_ps_update + num;
                end
                G = xor([0 G(1:end - 1)], G);
            end

            [PM_best, best_path] = min(PM);
            d_esti = U_path(:, best_path).';
            U_list = U_path(:, active_path);
        end
    
        function LLR_Ram_path = update_llr_path(obj, LLR_Ram_path, PS_RAM_path, i, S_MAX, STAGE, DEPTH, M)
            for s = S_MAX:-1:1
                if s < STAGE
                    Cs = 1;
                    rd_start_addr = DEPTH - 1 - s;
                    wr_start_addr = DEPTH - s;
                else
                    Cs = floor(2^(s-1) / M);
                    rd_start_addr = floor((obj.N - 2^s) / M);
                    wr_start_addr = floor((obj.N - 2^(s-1)) / M);
                end
        
                fg = bitget(i, s);
        
                for n = 0:Cs-1
                    % read LLR RAM
                    LLR_rd = LLR_Ram_path(rd_start_addr + 1, :);
                    rd_start_addr = rd_start_addr + 1;
        
                    % read partial sums if g function is needed
                    if fg == 1
                        PS_Rd_Addr = mod(floor((i - 2^(s-1)) / M) + n, obj.N/M/2);
                        beta = PS_RAM_path(PS_Rd_Addr + 1, :);
                        beta = FPSU_SCL.f_g_bits_selector(beta, i, s);
                    else
                        beta = zeros(1, M);
                    end
        
                    % execute M parallel f/g operations
                    for x = 1:M
                        if x <= 2^(s-1)
                            wr_data = FPSU_SCL.f_g(LLR_rd(2*x - 1), LLR_rd(2*x), fg, beta(x));
                        else
                            wr_data = 0;
                        end
        
                        source_addr = n * M + x - 1;
                        if source_addr < 2^(s-2)
                            trans_addr = source_addr * 2;
                        else
                            trans_addr = (source_addr - 2^(s-2)) * 2 + 1;
                        end
        
                        Base_addr   = floor(trans_addr / (2*M)) + wr_start_addr;
                        Offset_addr = mod(trans_addr, 2*M) + 1;
                        LLR_Ram_path(Base_addr + 1, Offset_addr) = wr_data;
                    end
                end
            end
        end
    
        function [PS_RAM_path, total_number_of_ps_update] = FPSU_path(obj, PS_RAM_path, bit, G, i, M)
            total_number_of_ps_update = 0;
            bit_and = and(repmat(bit, 1, M), G);
            if i == obj.N/2
                PS_RAM_path(:,:) = 0;
                PS_RAM_path(1,:) = bit_and;
                total_number_of_ps_update = total_number_of_ps_update + M;
                return;
            end
            if bit == 0
                return;
            end
            j0 = mod(floor(i / M), obj.N/M/2);
            j = j0;
            while true
                PS_RAM_path(j + 1, :) = xor(PS_RAM_path(j + 1, :), bit_and);
                total_number_of_ps_update = total_number_of_ps_update + M;
                if j == 0
                    break;
                end
                j = bitand(j - 1, j0);
            end
        end
    
    end

    methods(Static)

        function PM_new = update_pm(PM_old, llr, u)
            hard_bit = double(llr < 0);
            PM_new = PM_old + double(u ~= hard_bit) * abs(llr);
        end

        function LLR = f_g(L1, L2, fg, PSN_bit)
            if fg == 0
                LLR =  sign(L1)*sign(L2)*min(abs(L1), abs(L2));
            else
                if PSN_bit == 0
                    LLR = L1 + L2;
                else
                    LLR = L2 - L1;
                end
            end
        end

        function selected_bits = f_g_bits_selector(f_g_bits, i, s)
            M = numel(f_g_bits);
            selected_bits = f_g_bits;
            if s <= 0
                return;
            end
            group_len = 2^s;
            half_len  = 2^(s-1);
            if group_len > M
                return;
            end
            base = floor(mod(i, M) / group_len) * group_len;
            selected_bits(1:half_len) = f_g_bits(base + (1:half_len));
        end

        function depth = get_llr_ram_depth(N, M)
            depth = 0;
            t = 1;
            while t <= N
                if t <= M
                    depth = depth + 1;
                else
                    depth = depth + t/(2*M);
                end
                t = t * 2;
            end
        end
        
        function PM = calc_PM(PM, llr, u)
            if(u ~= 0.5*(1 - sign(llr)))
                PM = PM + abs(llr);
            end
        end
        
        function P = get_GN(N)
            F = [1, 0 ; 1, 1];
            P = zeros(N, N);
            P(1 : 2, 1 : 2) = F;
            for i = 2 : log2(N)
                P(1 : 2^i, 1 : 2^i) = kron(P(1 : 2^(i - 1), 1 : 2^(i - 1)), F);
            end
        end

    end
end

