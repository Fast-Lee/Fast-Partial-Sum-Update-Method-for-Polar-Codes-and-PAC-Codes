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

classdef FPSU_PAC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        N                   % code length
        K                   % info length
        gen                 % convolutional generator polynomial
        rate                
        rate_profiling      
        B_RAM               
        conv_depth          
        GN                  % Polar generator matrix
        convolution_matrix  % upper-triangular Toeplitz matrix
        llr
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        function obj = FPSU_PAC(N, k, g)
            n = ceil(log2(N));
            N = 2^n;
            obj.N = N;
            obj.K = k;
            obj.gen = g;
            obj.rate = k / N;
            obj.conv_depth = length(g);
            obj.rate_profiling = RM_rate_profiling(obj);
            obj.B_RAM = setB(obj);

            obj.GN = FPSU_PAC.get_GN(obj.N);                    % Polar generator matrix
            g_zp = [obj.gen, zeros(1, obj.N-obj.conv_depth)];
            obj.convolution_matrix = triu(toeplitz(g_zp));      % upper-triangular Toeplitz matrix
        end
        
        function info_indices = RM_rate_profiling(obj)
            %RM_rate_profiling 
            Channel_indices = (0:obj.N - 1)';
            bitStr = dec2bin(Channel_indices);
            bit = abs(bitStr) - 48;       % char -> int
            RM_score = sum(bit, 2);
            [~, sorted_indices] = sort(RM_score, 'ascend');
            info_indices = sort(sorted_indices(end-(obj.K)+1:end), 'ascend');
        end

        function B_RAM = setB(obj)
            B_RAM = zeros(1, obj.N);
            for i = 1: obj.N
                if ismember(i, obj.rate_profiling)
                    B_RAM(i) = 1;
                else
                    B_RAM(i) = 0;
                end
            end
        end
        
        function x = encode(obj, d)
            if (length(d) ~= obj.K)
                error('The length of the input d is not equal to k.')
            end
            
            v = zeros(1, obj.N);
            v(obj.rate_profiling) = d;
            % convolutional encoder
            u = mod(v*obj.convolution_matrix, 2);
            % Polar Encoding
            x = mod(u*obj.GN, 2)';
        end
        
        function [d_esti] = SC_decoder(obj, LLR, M)
            % calculate LLR RAM depth.
            DEPTH = 0;
            i = 1;
            while i <= obj.N
                if i <= M
                    DEPTH = DEPTH + 1;
                else
                    DEPTH = DEPTH + i/(2*M);
                end
                i = i * 2;
            end

            % address butterfly permutation
            LLR_Ram = zeros(DEPTH, 2*M);
            for i = 0: obj.N-1
                if (i < obj.N/2)
                    addr = i * 2;
                else
                    addr = (i - obj.N/2) * 2 + 1;
                end
                LLR_Ram(floor(addr/(2 * M))+1, mod(addr, (2 * M))+1) = LLR(i + 1); 
            end
            
            S_MAX = log2(obj.N);            % Maximum depth (or height) of the binary tree
            STAGE = log2(M) + 1;            
            info_cnt = 1;                   % Information bit counter
            d_esti = zeros(1, obj.K);       % decoding result
            G = [1 zeros(1, M - 1)];        % First row of Pascal's triangle
            PS_RAM = zeros(obj.N/M/2, M);   % Partial-sum RAM
            curr_state = zeros(1, length(obj.gen)-1); % Stored sequence of the convolutional code
            
            % decoding start
            for i = 0: obj.N-1
                for s = S_MAX: -1: 1        
                    if s < STAGE
                        Cs = 1;
                        rd_start_addr = DEPTH-1 - s;
                        wr_start_addr = DEPTH - s;
                    else
                        Cs = floor(2^(s-1) / M);
                        rd_start_addr = floor((obj.N-2^s)/M);
                        wr_start_addr = floor((obj.N-2^(s-1))/M);
                    end
                    
                    % Determine whether uses the f or g function based on the index i and stage s 
                    fg = FPSU_PAC.fg_selector(i, s);

                    for n = 0 : Cs-1
                        % read LLR RAM
                        LLR_rd = LLR_Ram(rd_start_addr+1 , :);
                        rd_start_addr = rd_start_addr + 1;
                        
                        if fg == 1
                            PS_Rd_Addr = mod(floor((i-2^(s-1))/M) + n, obj.N/M/2);
                            beta = PS_RAM(PS_Rd_Addr + 1, :);
                            if M == 1
                                beta = PS_RAM(PS_Rd_Addr + 1, :);
                            elseif M == 2
                                beta = FPSU_PAC.f_g_bits_selector2(beta, i, s);
                            elseif M == 4
                                beta = FPSU_PAC.f_g_bits_selector4(beta, i, s);
                            elseif M == 8
                                beta = FPSU_PAC.f_g_bits_selector8(beta, i, s);
                            elseif M == 16
                                beta = FPSU_PAC.f_g_bits_selector16(beta, i, s);
                            elseif M == 32
                                beta = FPSU_PAC.f_g_bits_selector32(beta, i, s);
                            end
                        else
                            beta = zeros(1, M);
                        end
            
                        % Execute M f- or g-functions per operation
                        for x = 1: M
                            if x <= 2^(s-1)
                                wr_data = FPSU_PAC.f_g(LLR_rd(2*x - 1), LLR_rd(2*x), fg, beta(x));
                            else
                                wr_data = 0;
                            end
                            source_addr = n * M + x - 1;
                            if(source_addr < 2^(s-2))
                                trans_addr = source_addr * 2;
                            else
                                trans_addr = (source_addr-2^(s-2)) * 2 + 1;
                            end
                            
                            Base_addr = floor(trans_addr/(2*M)) + wr_start_addr; 
                            Offset_addr = mod(trans_addr, 2*M) + 1;            
                            LLR_Ram(Base_addr+1, Offset_addr) = wr_data ;
                        end
                    end
                end
                
                % get vi and ui
                if obj.B_RAM(i+1) == 0        % frozen-bit
                    [u_hat, curr_state] = FPSU_PAC.conv1b(0, curr_state , obj.gen);
                elseif obj.B_RAM(i+1) == 1    % infomation-bit
                    PM_pair = zeros(1, 2);
                    curr_state_temp  = curr_state ;
                    [u_left , curr_state] = FPSU_PAC.conv1b(0, curr_state , obj.gen);
                    if (LLR_Ram(DEPTH, 1) < 0) ~= u_left
                        PM_pair(1) = abs(LLR_Ram(DEPTH, 1));
                    end
                    [u_right , curr_state_temp] = FPSU_PAC.conv1b(1, curr_state_temp , obj.gen);
                    if (LLR_Ram(DEPTH, 1) < 0) ~= u_right
                        PM_pair(2) = abs(LLR_Ram(DEPTH, 1));
                    end
                    
                    if PM_pair(1) > PM_pair(2)
                        d_esti(info_cnt) = 1;
                        u_hat = u_right;
                        curr_state  = curr_state_temp ;
                    else 
                        d_esti(info_cnt) = 0;
                        u_hat = u_left;
                    end
                    info_cnt = info_cnt + 1;
                end
                
                % Starting tree level for decoding the next bit
                S_MAX = FPSU_PAC.count_one(i);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % update PS RAM
                bit_and = and(repmat(u_hat, 1, M), G);
  
                if i == obj.N/2
                    PS_RAM = zeros(obj.N/M/2, M);
                    PS_RAM(1, :) = bit_and;
                else
                    j = mod(floor(i / M), obj.N/M/2);
                    while (u_hat ~= 0)
                        Ram_Val = PS_RAM(j + 1, :);
                        PS_RAM(j + 1, :) = xor(Ram_Val, bit_and);
                        if j == 0
                            break;
                        end
                        j = bitand(j - 1, mod(floor(i / M), obj.N/M/2));
                    end
                end
                G = xor([0 G(1: end - 1)], G);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end

    
    end    
        

    methods(Static)

        function [u, next_state] = conv1b(v, curr_state, g)
            u = v*g(1);
            for i = 2: length(g)
                u = xor(u, g(i)*curr_state(i-1));
            end
            
            next_state=[v curr_state(1:end-1)];
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


        % Determine whether stage s requires an f or g operation
        function fg = fg_selector(i, s)
            lst = zeros(1, 10);
            cnt = 1;
            while(i > 0)
                lst(cnt) = mod(i, 2);
                i = floor(i / 2);
                cnt = cnt + 1;
            end
            fg = lst(s);
        end


        % f-function and g-function
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
        
        function S_max = count_one(i)
            x = i;
            y = i + 1;
            
            lst_x = zeros(1, 12);
            lst_y = zeros(1, 12);
            
            cnt = 1;
            while(x > 0)
                lst_x(cnt) = mod(x , 2);
                cnt = cnt + 1;
                x = floor(x / 2);
            end
            
            cnt = 1;
            while(y > 0)
                lst_y(cnt) = mod(y , 2);
                cnt = cnt + 1;
                y = floor(y / 2);
            end
            
            S_max = 0;
            xor_ = xor(lst_x, lst_y);
            for i = 1: 11
               if xor_(i) == 1
                   S_max = S_max + 1;
               end
            end
        end

        % Output the required partial sums for the current stage from 4 bits
        function selected_bits = f_g_bits_selector2(f_g_bits, i, s)
            selected_bits = zeros(1, 2);
            selected_bits(1) = f_g_bits(1);
            selected_bits(2) = f_g_bits(2);
        end

        % Output the required partial sums for the current stage from 4 bits
        function selected_bits = f_g_bits_selector4(f_g_bits, i, s)
            selected_bits = zeros(1, 4);
            % s=1
            if mod(i, 4) == 1
                b0_1 = f_g_bits(0+1);
            elseif mod(i, 4) == 3
                b0_1 = f_g_bits(2+1);
            end
           
            if s == 1
                selected_bits(1) = b0_1;
            else
                selected_bits(1) = f_g_bits(1);
            end
            
            selected_bits(2) = f_g_bits(2);
            selected_bits(3) = f_g_bits(3);
            selected_bits(4) = f_g_bits(4);

        end


        % Output the required partial sums for the current stage from 8 bits
        function selected_bits = f_g_bits_selector8(f_g_bits, i, s)
            
            selected_bits = zeros(1, 8);
            
            %%%%%%%%%%%%%%% bit 1 %%%%%%%%%%%%%%%
            % s=1
            if mod(i, 8) == 1
                b0_1 = f_g_bits(0+1);
            elseif mod(i, 8) == 3
                b0_1 = f_g_bits(2+1);
            elseif mod(i, 8) == 5
                b0_1 = f_g_bits(4+1);
            elseif mod(i, 8) == 7
                b0_1 = f_g_bits(6+1);
            end
            
            % s=2
            if mod(i, 8) == 2
                b0_2 = f_g_bits(0+1);
            elseif mod(i, 8) == 6
                b0_2 = f_g_bits(4+1);
            end
            
            if s == 1
                selected_bits(0+1) = b0_1;
            elseif s == 2
                selected_bits(0+1) = b0_2;
            else
                selected_bits(0+1) = f_g_bits(0+1);
            end
            
            %%%%%%%%%%%%%%% bit 2 %%%%%%%%%%%%%%%
            % s=2
            if mod(i, 8) == 2
                b1_2 = f_g_bits(1+1);
            elseif mod(i, 8) == 6
                b1_2 = f_g_bits(5+1);
            end

            if s == 2
                selected_bits(1+1) = b1_2;
            else
                selected_bits(1+1) = f_g_bits(1+1);
            end
            
            %%%%%%%%%%%%%%% others %%%%%%%%%%%%%%%
            selected_bits(2 + 1) = f_g_bits(2+1);
            selected_bits(3 + 1) = f_g_bits(3+1);
            selected_bits(4 + 1) = f_g_bits(4+1);
            selected_bits(5 + 1) = f_g_bits(5+1);
            selected_bits(6 + 1) = f_g_bits(6+1);
            selected_bits(7 + 1) = f_g_bits(7+1);
        end
        

        % Output the required partial sums for the current stage from 16 bits
        function selected_bits = f_g_bits_selector16(f_g_bits, i, s)
            
            selected_bits = zeros(1, 16);
            
            %%%%%%%%%%%%%%% bit 1 %%%%%%%%%%%%%%%
            % s=1
            if mod(i, 16) == 1
                b0_1 = f_g_bits(0+1);
            elseif mod(i, 16) == 3
                b0_1 = f_g_bits(2+1);
            elseif mod(i, 16) == 5
                b0_1 = f_g_bits(4+1);
            elseif mod(i, 16) == 7
                b0_1 = f_g_bits(6+1);
            elseif mod(i, 16) == 9
                b0_1 = f_g_bits(8+1);
            elseif mod(i, 16) == 11
                b0_1 = f_g_bits(10+1);
            elseif mod(i, 16) == 13
                b0_1 = f_g_bits(12+1);
            else
                b0_1 = f_g_bits(14+1);
            end
            
            % s=2
            if mod(i, 16) == 2
                b0_2 = f_g_bits(0+1);
            elseif mod(i, 16) == 6
                b0_2 = f_g_bits(4+1);
            elseif mod(i, 16) == 10
                b0_2 = f_g_bits(8+1);
            else
                b0_2 = f_g_bits(12+1);
            end
            
            % s=3
            if mod(i, 16) == 4
                b0_3 = f_g_bits(0+1);
            else
                b0_3 = f_g_bits(8+1);
            end
            
            if s == 1
                selected_bits(0+1) = b0_1;
            elseif s == 2
                selected_bits(0+1) = b0_2;
            elseif s == 3
                selected_bits(0+1) = b0_3;
            else
                selected_bits(0+1) = f_g_bits(0+1);
            end
            
            %%%%%%%%%%%%%%% bit 2 %%%%%%%%%%%%%%%
            % s=2
            if mod(i, 16) == 2
                b1_2 = f_g_bits(1+1);
            elseif mod(i, 16) == 6
                b1_2 = f_g_bits(5+1);
            elseif mod(i, 16) == 10
                b1_2 = f_g_bits(9+1);
            else
                b1_2 = f_g_bits(13+1);
            end
            
            % s=3
            if mod(i, 16) == 4
                b1_3 = f_g_bits(1+1);
            else
                b1_3 = f_g_bits(9+1);
            end

            if s == 2
                selected_bits(1+1) = b1_2;
            elseif s == 3
                selected_bits(1+1) = b1_3;
            else
                selected_bits(1+1) = f_g_bits(1+1);
            end
            
            %%%%%%%%%%%%%%% bit 3 %%%%%%%%%%%%%%%
            % s=3
            if mod(i, 16) == 4
                b2_3 = f_g_bits(2+1);
            else
                b2_3 = f_g_bits(10+1);
            end

            if s == 3
                selected_bits(2+1) = b2_3;
            else
                selected_bits(2+1) = f_g_bits(2+1);
            end
            
            %%%%%%%%%%%%%%% bit 4 %%%%%%%%%%%%%%%
            % s=3
            if mod(i, 16) == 4
                b3_3 = f_g_bits(3+1);
            else
                b3_3 = f_g_bits(11+1);
            end

            if s == 3
                selected_bits(3+1) = b3_3;
            else
                selected_bits(3+1) = f_g_bits(3+1);
            end
            
            %%%%%%%%%%%%%%% others %%%%%%%%%%%%%%%
            selected_bits(4 + 1) = f_g_bits(4+1);
            selected_bits(5 + 1) = f_g_bits(5+1);
            selected_bits(6 + 1) = f_g_bits(6+1);
            selected_bits(7 + 1) = f_g_bits(7+1);
            selected_bits(8 + 1) = f_g_bits(8+1);
            selected_bits(9 + 1) = f_g_bits(9+1);
            selected_bits(10 + 1) = f_g_bits(10+1);
            selected_bits(11 + 1) = f_g_bits(11+1);
            selected_bits(12 + 1) = f_g_bits(12+1);
            selected_bits(13 + 1) = f_g_bits(13+1);
            selected_bits(14 + 1) = f_g_bits(14+1);
            selected_bits(15 + 1) = f_g_bits(15+1); 
        end
        
        % Output the required partial sums for the current stage from 32 bits
        function selected_bits = f_g_bits_selector32(f_g_bits, i, s)
            selected_bits = zeros(1, 32);
            
            %%%%%%%%%%%%%%% bit 1 %%%%%%%%%%%%%%%
            % s=1
            if mod(i, 32) == 1
                b0_1 = f_g_bits(0+1);
            elseif mod(i, 32) == 3
                b0_1 = f_g_bits(2+1);
            elseif mod(i, 32) == 5
                b0_1 = f_g_bits(4+1);
            elseif mod(i, 32) == 7
                b0_1 = f_g_bits(6+1);
            elseif mod(i, 32) == 9
                b0_1 = f_g_bits(8+1);
            elseif mod(i, 32) == 11
                b0_1 = f_g_bits(10+1);
            elseif mod(i, 32) == 13
                b0_1 = f_g_bits(12+1);
            elseif mod(i, 32) == 15
                b0_1 = f_g_bits(14+1);
            elseif mod(i, 32) == 17
                b0_1 = f_g_bits(16+1);
            elseif mod(i, 32) == 19
                b0_1 = f_g_bits(18+1);
            elseif mod(i, 32) == 21
                b0_1 = f_g_bits(20+1);
            elseif mod(i, 32) == 23
                b0_1 = f_g_bits(22+1);
            elseif mod(i, 32) == 25
                b0_1 = f_g_bits(24+1);
            elseif mod(i, 32) == 27
                b0_1 = f_g_bits(26+1);
            elseif mod(i, 32) == 29
                b0_1 = f_g_bits(28+1);
            else%if mod(i, 32) == 31
                b0_1 = f_g_bits(30+1);
            end
            
            % s=2
            if mod(i, 32) == 2
                b0_2 = f_g_bits(0+1);
            elseif mod(i, 32) == 6
                b0_2 = f_g_bits(4+1);
            elseif mod(i, 32) == 10
                b0_2 = f_g_bits(8+1);
            elseif mod(i, 32) == 14
                b0_2 = f_g_bits(12+1);
            elseif mod(i, 32) == 18
                b0_2 = f_g_bits(16+1);
            elseif mod(i, 32) == 22
                b0_2 = f_g_bits(20+1);
            elseif mod(i, 32) == 26
                b0_2 = f_g_bits(24+1);
            else%if mod(i, 32) == 30
                b0_2 = f_g_bits(28+1);
            end
            
            % s=3
            if mod(i, 32) == 4
                b0_3 = f_g_bits(0+1);
            elseif mod(i, 32) == 12
                b0_3 = f_g_bits(8+1);
            elseif mod(i, 32) == 20
                b0_3 = f_g_bits(16+1);
            else%if mod(i, 32) == 28
                b0_3 = f_g_bits(24+1);
            end
            
            
            % s=4
            if mod(i, 32) == 8
                b0_4 = f_g_bits(0+1);
            else%if mod(i, 32) == 24
                b0_4 = f_g_bits(16+1);
            end

            if s == 1
                selected_bits(0+1) = b0_1;
            elseif s == 2
                selected_bits(0+1) = b0_2;
            elseif s == 3
                selected_bits(0+1) = b0_3;
            elseif s == 4
                selected_bits(0+1) = b0_4;
            else
                selected_bits(0+1) = f_g_bits(0+1);
            end
            
            %%%%%%%%%%%%%%% bit 2 %%%%%%%%%%%%%%%
            % s=2
            if mod(i, 32) == 2
                b1_2 = f_g_bits(1+1);
            elseif mod(i, 32) == 6
                b1_2 = f_g_bits(5+1);
            elseif mod(i, 32) == 10
                b1_2 = f_g_bits(9+1);
            elseif mod(i, 32) == 14
                b1_2 = f_g_bits(13+1);
            elseif mod(i, 32) == 18
                b1_2 = f_g_bits(17+1);
            elseif mod(i, 32) == 22
                b1_2 = f_g_bits(21+1);
            elseif mod(i, 32) == 26
                b1_2 = f_g_bits(25+1);
            else%if mod(i, 32) == 30
                b1_2 = f_g_bits(29+1);
            end

            % s=3
            if mod(i, 32) == 4
                b1_3 = f_g_bits(1+1);
            elseif mod(i, 32) == 12
                b1_3 = f_g_bits(9+1);
            elseif mod(i, 32) == 20
                b1_3 = f_g_bits(17+1);
            else%if mod(i, 32) == 28
                b1_3 = f_g_bits(25+1);
            end
            
            % s=4
            if mod(i, 32) == 8
                b1_4 = f_g_bits(1+1);
            else%if mod(i, 32) == 24
                b1_4 = f_g_bits(17+1);
            end

            if s == 2
                selected_bits(1+1) = b1_2;
            elseif s == 3
                selected_bits(1+1) = b1_3;
            elseif s == 4
                selected_bits(1+1) = b1_4;
            else
                selected_bits(1+1) = f_g_bits(1+1);
            end
            
            %%%%%%%%%%%%%%% bit 3 %%%%%%%%%%%%%%%
            % s=3
            if mod(i, 32) == 4
                b2_3 = f_g_bits(2+1);
            elseif mod(i, 32) == 12
                b2_3 = f_g_bits(10+1);
            elseif mod(i, 32) == 20
                b2_3 = f_g_bits(18+1);
            else%if mod(i, 32) == 28
                b2_3 = f_g_bits(26+1);
            end
            
            % s=4
            if mod(i, 32) == 8
                b2_4 = f_g_bits(2+1);
            else%if mod(i, 32) == 24
                b2_4 = f_g_bits(18+1);
            end
            
            if s == 3
                selected_bits(2+1) = b2_3;
            elseif s == 4
                selected_bits(2+1) = b2_4;
            else
                selected_bits(2+1) = f_g_bits(2+1);
            end
            
            %%%%%%%%%%%%%%% bit 4 %%%%%%%%%%%%%%%
            % s=3
            if mod(i, 32) == 4
                b3_3 = f_g_bits(3+1);
            elseif mod(i, 32) == 12
                b3_3 = f_g_bits(11+1);
            elseif mod(i, 32) == 20
                b3_3 = f_g_bits(19+1);
            else%if mod(i, 32) == 28
                b3_3 = f_g_bits(27+1);
            end
            
            % s=3
            if mod(i, 32) == 8
                b3_4 = f_g_bits(3+1);
            else%if mod(i, 32) == 24
                b3_4 = f_g_bits(19+1);
            end

            if s == 3
                selected_bits(3+1) = b3_3;
            elseif s == 4
                selected_bits(3+1) = b3_4;
            else
                selected_bits(3+1) = f_g_bits(3+1);
            end
            
            %%%%%%%%%%%%%%% bit 5 %%%%%%%%%%%%%%%
            % s=4
            if mod(i, 32) == 8
                b4_4 = f_g_bits(4+1);
            else%if mod(i, 32) == 24
                b4_4 = f_g_bits(20+1);
            end

            if s == 4
                selected_bits(4+1) = b4_4;
            else
                selected_bits(4+1) = f_g_bits(4+1);
            end
            
            %%%%%%%%%%%%%%% bit 6 %%%%%%%%%%%%%%%
            % s=4
            if mod(i, 32) == 8
                b5_4 = f_g_bits(5+1);
            else%if mod(i, 32) == 24
                b5_4 = f_g_bits(21+1);
            end
            
            if s == 4
                selected_bits(5+1) = b5_4;
            else
                selected_bits(5+1) = f_g_bits(5+1);
            end
            
            %%%%%%%%%%%%%%% bit 7 %%%%%%%%%%%%%%%
            % s=4
            if mod(i, 32) == 8
                b6_4 = f_g_bits(6+1);
            else%if mod(i, 32) == 24
                b6_4 = f_g_bits(22+1);
            end

            if s == 4
                selected_bits(6+1) = b6_4;
            else
                selected_bits(6+1) = f_g_bits(6+1);
            end
            
            %%%%%%%%%%%%%%% bit 8 %%%%%%%%%%%%%%%
            % s=4
            if mod(i, 32) == 8
                b7_4 = f_g_bits(7+1);
            else%if mod(i, 32) == 24
                b7_4 = f_g_bits(23+1);
            end

            if s == 4
                selected_bits(7+1) = b7_4;
            else
                selected_bits(7+1) = f_g_bits(7+1);
            end
            
            %%%%%%%%%%%%%%% other bits %%%%%%%%%%%%%%%
            selected_bits(8+1) = f_g_bits(8+1);
            selected_bits(9+1) = f_g_bits(9+1);
            selected_bits(10+1) = f_g_bits(10+1);
            selected_bits(11+1) = f_g_bits(11+1);
            selected_bits(12+1) = f_g_bits(12+1);
            selected_bits(13+1) = f_g_bits(13+1);
            selected_bits(14+1) = f_g_bits(14+1);
            selected_bits(15+1) = f_g_bits(15+1);
            selected_bits(16+1) = f_g_bits(16+1);
            selected_bits(17+1) = f_g_bits(17+1);
            selected_bits(18+1) = f_g_bits(18+1);
            selected_bits(19+1) = f_g_bits(19+1);
            selected_bits(20+1) = f_g_bits(20+1);
            selected_bits(21+1) = f_g_bits(21+1);
            selected_bits(22+1) = f_g_bits(22+1);
            selected_bits(23+1) = f_g_bits(23+1);
            selected_bits(24+1) = f_g_bits(24+1);
            selected_bits(25+1) = f_g_bits(25+1);
            selected_bits(26+1) = f_g_bits(26+1);
            selected_bits(27+1) = f_g_bits(27+1);
            selected_bits(28+1) = f_g_bits(28+1);
            selected_bits(29+1) = f_g_bits(29+1);
            selected_bits(30+1) = f_g_bits(30+1);
            selected_bits(31+1) = f_g_bits(31+1);
        end
    end

end
