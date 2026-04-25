
% This project primarily implements the proposed fast partial-sum update (FPSU)
% algorithm and compares the decoding performance and total number of partial-sum
% updates of the FPSU Polar SCD and FPSU PAC SCD with conventional Polar and 
% PAC codes, thereby demonstrating the advantages of our method.

% The project mainly consists of 5 folders, organized as follows:

%% - Polar
%     -- Polar_SC
%         --- Polar_Conventional_SC.m     : Conventional Polar codes using (3) as the partial-sum update strategy.
%         --- BLER.m                      : Block-error rate of Polar_Conventional_SC.
%         --- Average_Num_of_PS_Update.m  : Simulate the impact of different R, values of M, and N, and Eb/N0 on the number of partial-sum bits update.

%     -- Fast_Polar_Decoders
%         --- Fast_Polar.m                : Fast Polar decoders in [10].
%         --- BLER.m                      : Block-error rate of Fast Polar decoders in [10].
%         --- Average_Num_of_PS_Update.m  : Simulate the impact of different R, values of M, and N, and Eb/N0 on the number of partial-sum bits update.

%     -- Fast_SCD
%         --- Fast_Polar_SC2.m            : Fast Polar SCD in [11].
%         --- BLER.m                      : Block-error rate of Fast Polar SCD in [11].
%         --- Average_Num_of_PS_Update.m  : Simulate the impact of different R, values of M, and N, and Eb/N0 on the number of partial-sum bits update.
 
%     -- FPSU_Polar_SC
%         --- FPSU_SC.m                   : Proposed fast partial-sum update (FPSU) algorithm used in Polar SCD.
%         --- BLER.m                      : Block-error rate of FPSU Polar SCD.
%         --- Average_Num_of_PS_Update.m  : Simulate the impact of different R, values of M, and N, and Eb/N0 on the number of partial-sum bits update.

%     -- FPSU_Polar_SCL
%         --- FPSU_SCL.m                  : Proposed fast partial-sum update (FPSU) algorithm used in Polar SCL.
%         --- BLER.m                      : Block-error rate of FPSU Polar SCL.
%         --- Average_Num_of_PS_Update.m  : Simulate the impact of different R, values of M, and N, and Eb/N0 on the number of partial-sum bits update.

%     -- Polar_SCL
%         --- PolarSCLRef.m               : Polar SCL algorithm for reference.
%         --- BLER.m                      : Block-error rate of FPSU Polar SCL.
%         --- Average_Num_of_PS_Update.m  : Simulate the impact of different R, values of M, and N, and Eb/N0 on the number of partial-sum bits update.


%% - PAC
%     -- PACCodes
%         --- PACCodes.m                  : PAC SCD and list decoders in [2].
%         --- BLER.m                      : Block-error rate of PAC SCD and list decoders in [2].
%         --- Average_Num_of_PS_Update.m  : Simulate the impact of different R, values of M, and N, and Eb/N0 on the number of partial-sum bits update.
% 
%     -- FPSU_PAC
%         --- FPSU_PAC.m                  : Proposed fast partial-sum update (FPSU) algorithm used in PAC SCD.
%         --- BLER.m                      : Block-error rate of FPSU PAC SCD.
%         --- Average_Num_of_PS_Update.m  : Simulate the impact of different R, values of M, and N, and Eb/N0 on the number of partial-sum bits update.
 
%% - BLER_FPSU_SC
%         --- Draw_Fig3.m                 : Comparasion of Block error rate

%% - BLER_FPSU_SCL
%         --- Draw_Fig.m                  : Comparasion of Block error rate

%% - ANPSU
%         --- Complexity_Analysis.m       : Comparison of Partial-Sum Bit Updates With and Without Algorithm 1 at Code Rate R = 1
%         --- Polar_ANPSU.m               : Reduction Ratio of FPSU SCD Compared With other algorithms for Polar Codes
%         --- PolarSCL_ANPSU.m            : Reduction Ratio of FPSU SCL Compared With Polar SCL
%         --- PAC_ANPSU.m                 : Reduction Ratio of FPSU SCD Compared With other algorithms for PAC Codes
