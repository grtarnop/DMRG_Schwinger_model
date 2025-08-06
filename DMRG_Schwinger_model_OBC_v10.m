% Dempsey et all convention
% v8. 09/28/2024. (new version with updated functions)
% v9. 10/05/2024. (improved guessing eigs functions for sweeps)
% v10. 10/19/2024. (saves MPS)
% This is a finite DMRG main script with spin projection 
% for the Schwinger model (Dempsey et all convention)
% In this code 1 DMRG site = 2 spin-1/2 sites
clear all
clc
format longG

% The Schwinger model parameters
x = 50;       % gauge coupling x = 1/(ea)^2
me = 0.2;    % m/e
shift = 1;   % mass shift 
mu = 2*sqrt(x)*(me) - 1/4*shift;  % mass term mu = 2m/(e^2 a)
y = 1;        % parameter to regulate Coulomb energy term
alpha =  1;
theta = pi*alpha;   % theta angle

% half of the total DMRG sites in iDMRG procedure 
% (4*N qubits in the whole chain)
N = 250;        
m_inf = 100;   % number of states kept under truncation for infinite DMRG
m_fin = 400;    % number of states kept under truncation for finite DMRG

% If dynm = 1 then m is chosen such that trunc_error < trun_tol
% If dynm = 0 then m is fixed to the initial value
dynm_inf = 0;   % switch parameter for infinite DMRG part
dynm_fin = 1;   % switch parameter for finite DMRG part
% truncation error tolerance (used when m is dynamical)
% should be bigger than 5*1e-15 
trun_tol = 1e-14; 

eig_tol = 1e-10; % eigenvalues accuracy
maxit = 1;       % maximum number of iterations for fDMRG
kryl = 10;       % dimension of the Krylov space

lvs = 1;        % number of lowest energy levels
s =  0;         % total spin projection of the superblock
sweeps = 2;     % number of sweeps
lev = 1;        % energy level (ground state is lev = 1)

% Initialize lists of N empty cells for all left and right objects
% Left blocks and their quantum sectors dimensions:
HL = cell(N, 1); QL = cell(N, 1); SL = cell(N, 1); 
SpL = cell(N, 1); mL = cell(N, 1);  
% Truncated left blocks and their quantum sectors dimensions:
HLb = cell(N, 1); QLb = cell(N, 1); SLb = cell(N, 1); 
SpLb = cell(N, 1); mLb = cell(N, 1);
% Unitary matrices which represent truncation and ledgers:
UL = cell(N, 1); TL = cell(N, 1);
% Right blocks and their quantum sectors dimensions:
HR = cell(N, 1); QR = cell(N, 1); SR = cell(N, 1) ; 
SpR = cell(N, 1); mR = cell(N, 1);
% Rruncated right blocks and their quantum sectors dimensions:
HRb = cell(N, 1);  QRb = cell(N, 1); SRb = cell(N, 1);  
SpRb = cell(N, 1); mRb = cell(N, 1); 
% Unitary matrices which represent truncation and ledgers:
UR = cell(N, 1); TR = cell(N, 1);

% Initialize a list of N empty cells for wave functions in LR form
psiLR = cell(N, 1);

% Initialize a list of N empty cells for indices of quantum combinations 
idx_qc = cell(N, 1);
idx_qcL = cell(N, 1);
idx_qcR = cell(N, 1);

% Build first DMRG site for the Schwinger left and right blocks
[HLb{1}, QLb{1}, SLb{1}, SpLb{1}, mLb{1}, ...
 HRb{1}, QRb{1}, SRb{1}, SpRb{1}, mRb{1}] = ...
                  build_first_DMRG_site_Schwinger_blocks(x, mu, y, theta);

                                   
startTime = datetime('now');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Start infinite DMRG loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['****************************************'])
disp(['********* Infinite DMRG begins *********'])
disp(['****************************************'])

% List of energies for infinite DMRG part and sweeps
list_of_energies = [];  
% Set initial value of the maximum bond dimension of the iDMRG
m_max = 4;

for n = 2:N
disp(['========================================'])    
disp(['Number of DMRG sites of the whole chain = ', num2str(2*n)])    
disp(['========================================'])  

%%%%  Build all operators as q-q cells for n DMRG sites %%%%%%%%%%%%%%%%%%%
tic
[HL{n}, QL{n}, SL{n}, SpL{n}, mL{n}, TL{n}] = ...
          build_Schwinger_blocks_left(HLb{n - 1}, QLb{n-1}, SLb{n - 1}, ...
                       SpLb{n - 1}, mLb{n - 1}, n, N, s,  x, mu, y, theta);  
toc

tic
[HR{n}, QR{n}, SR{n}, SpR{n}, mR{n}, TR{n}] = ...
            build_Schwinger_blocks_right(HRb{n - 1}, QRb{n-1}, SRb{n - 1}, ...
                       SpRb{n - 1}, mRb{n - 1}, n, N, s, x, mu, y, theta);  
toc
%%%% Find ground eigenstate(s) and eigenenergy(ies) of the superblock %%%%%    
% Start DMRG procedure only if the current size of the chain can have 
% required total spin projection s
if n >= abs(s)
tic
[psiLR{n}, e0, idx_qc{n}] = eigs_Schwinger_superblock(HL{n}, QL{n}, SpL{n},...
                                    mL{n}, HR{n}, SR{n}, SpR{n}, mR{n}, ...
                                              [n, n], s, x, y, theta, ...
                                              lvs, kryl, eig_tol, maxit);   
disp(['Energies of the target states = '])
E0 = (diag(e0) - mu*N*2)/x
list_of_energies = [list_of_energies, [toc; 2*n; E0]];
toc

%%%%%%%%%% Truncate left and right blocks %%%%%%%%%
tic
[HLb{n}, QLb{n}, SLb{n}, SpLb{n}, UL{n}, mLb{n}, ...
 HRb{n}, QRb{n}, SRb{n}, SpRb{n}, UR{n}, mRb{n}, m_current] = ...
                 truncate_Schwinger_blocks_iDMRG(HL{n}, QL{n}, SL{n}, ...
                     SpL{n}, mL{n}, HR{n}, QR{n}, SR{n}, SpR{n}, mR{n}, ...   
                            psiLR{n}, [n, n], s, m_inf, dynm_inf, trun_tol);
% Update maximum bond dimension m_max
if m_current > m_max
    m_max = m_current;
end
toc
else
% Assign trunctated operators to be equal to the non-truncated
% Left block operators
HLb{n} = HL{n}; QLb{n} = QL{n}; SLb{n} = SL{n}; SpLb{n} = SpL{n};
% Left operators quantum sectors dimensions
mLb{n} = mL{n};
% Right block operators
HRb{n} = HR{n}; QRb{n} = QR{n}; SRb{n} = SR{n}; SpRb{n} = SpR{n};  
% Right operators quantum sectors dimensions
mRb{n} = mR{n}; 
end

end                                        
endTime = datetime('now');
elapsedTime = seconds(endTime - startTime)


% Create separately left and right indices of quantum combinations  
idx_qcL = idx_qc;
idx_qcR = idx_qc;
% Save MPS data after infinite DMRG into a speparate file
MPS_iDMRG_file_name_save = strcat('MPS_Schwinger_OBC_N', num2str(4*N),'_x', num2str(x), ...
   '_me', num2str(me),'_shf', num2str(shift),'_th', num2str(alpha), ...
   '_s', num2str(s), '_swps',num2str(0), '_m', num2str(m_max), ...
   '_lev', num2str(lev),'.mat');
save(MPS_iDMRG_file_name_save, 'UL','mL','mLb', 'UR','mR','mRb', ...
                           'psiLR', 'idx_qcL', 'idx_qcR'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Start finite DMRG loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Arrays of computation times
build_times = [];
eigs_times = [];
trunc_times = [];
sweep_times = [];

% Array of all energies for finite DMRG part
list_of_all_energies = [];  
% Set initial value of the maximum bond dimension of the fDMRG
m_max = 4;


disp(['**************************************'])
disp(['********* Finite DMRG begins *********'])
disp(['**************************************'])   
for sweep = 1:sweeps
startTime = datetime('now');

disp(['Sweep number = ', num2str(sweep)])       
%%%%%% Half sweep from left to right   ------------->   %%%%%%%%%%%%%%%%%%%
    for n = (N + 1):(2*N - 2)
        disp(['----> Sweep ',num2str(sweep), ...
              ' from left to right ----> Chain size:  (', ...
              num2str(n),',', num2str(2*N - n),')' ])
        tic       
        [HL{n}, QL{n}, SL{n}, SpL{n}, mL{n}, TL{n}] = ...
          build_Schwinger_blocks_left(HLb{n - 1}, QLb{n - 1}, SLb{n - 1}, ...
                       SpLb{n - 1}, mLb{n - 1}, n, N, s, x, mu, y, theta); 
        toc
        build_times = [build_times, [n ; 2*N - n; toc ] ];

        tic
%%%%%%% Find ground eigenstate and eigenenergy of the superblock %%%%%%%%%%
        % Use proper combination of quantum numbers for the first quarter
        % of the first sweep
        if sweep == 1
            idx_qc_old = idx_qc{2*N - n};
        else
            idx_qc_old = idx_qc{n};
        end
        [psiLR{n}, e0, idx_qc{n}] = ...
             eigs_Schwinger_superblock_guess_ltor(HL{n}, QL{n}, SpL{n}, mL{n},...
                 HR{2*N - n}, SR{2*N - n}, SpR{2*N - n}, mR{2*N - n}, ...
                 UL{n - 1}, TL{n}, UR{2*N - n}, TR{2*N - n + 1}, {psiLR{n-1}{1,:}},...
                 idx_qc{n - 1}, idx_qc_old, mRb{2*N - n}, [n, 2*N - n], ...
                                s, x, y, theta, lvs, kryl, eig_tol, maxit);     
        disp(['Energies of the target states = '])
        E0 = (diag(e0) - mu*N*2)/x
        list_of_all_energies = horzcat(list_of_all_energies, [n; 2*N-n; E0]);
        toc
        eigs_times = [eigs_times, [n ; 2*N - n; toc ]];
        
%%%%%%% Truncate left and right blocks %%%%%%%%%  
       tic
       [HLb{n}, QLb{n}, SLb{n}, SpLb{n}, UL{n}, mLb{n}, ...
        HRb{2*N - n}, QRb{2*N - n}, SRb{2*N - n}, SpRb{2*N - n},...
                                    UR{2*N - n}, mRb{2*N - n}, m_current] = ...
                truncate_Schwinger_blocks_fDMRG(HL{n}, QL{n}, SL{n},...
                            SpL{n}, mL{n}, HR{2*N - n}, QR{2*N - n}, ...
                            SR{2*N - n}, SpR{2*N - n}, mR{2*N - n}, ...
                                    psiLR{n}, idx_qc{n}, [n, 2*N - n], ...
                                        s, m_fin, dynm_fin, trun_tol); 
       % Update maximum bond dimension m_max
       if m_current > m_max
           m_max = m_current;
       end                            
       toc
       trunc_times = [trunc_times,  [n ; 2*N - n; m_current; toc ]];
    end
    

%%%%%%  Half sweep from right to left  <------------- %%%%%%%%%%%%%%%%%%%%%
    for n = 3:(2*N - 2)
        disp(['<---- Sweep ', num2str(sweep), ...
             ' from right to left <---- Chain size: (', ...
            num2str(2*N - n),',', num2str(n),')' ])
        tic
        [HR{n}, QR{n}, SR{n}, SpR{n}, mR{n}, TR{n}] = ...
            build_Schwinger_blocks_right(HRb{n - 1}, QRb{n-1}, SRb{n - 1}, ...
                       SpRb{n - 1}, mRb{n - 1}, n, N, s, x, mu, y, theta); 
        toc
        build_times = [build_times, [2*N - n ; n; toc]];

%%%%%%%%%% Find ground eigenstate and eigenenergy of the superblock %%%%%%%
        tic
        [psiLR{2*N - n}, e0, idx_qc{2*N - n}] = ...
                  eigs_Schwinger_superblock_guess_rtol(HL{2*N - n}, ...
                        QL{2*N - n}, SpL{2*N - n}, mL{2*N - n}, ...
                        HR{n}, SR{n}, SpR{n}, mR{n}, ...
                        UL{2*N - n}, TL{2*N - n + 1}, UR{n - 1}, TR{n}, ...
                        {psiLR{2*N - n + 1}{1,:}}, idx_qc{2*N - n + 1}, ...
                        idx_qc{2*N - n}, mLb{2*N - n}, ...
                        [2*N - n, n], s,  x, y, theta, ...
                        lvs, kryl, eig_tol, maxit);
        disp(['Energies of the target states = '])
        E0 = (diag(e0) - mu*N*2)/x
        list_of_all_energies = horzcat(list_of_all_energies, [2*N-n; n; E0]);
        toc
        eigs_times = [eigs_times,  [2*N - n ; n; toc]];
  %%%%%%%%%% Truncate left and right blocks %%%%%%%%%  
    tic
    [HLb{2*N - n}, QLb{2*N - n}, SLb{2*N - n}, SpLb{2*N - n}, ...
        UL{2*N - n}, mLb{2*N - n}, HRb{n}, QRb{n}, SRb{n}, SpRb{n},...
                                    UR{n}, mRb{n}, m_current] = ...
                truncate_Schwinger_blocks_fDMRG(HL{2*N - n}, QL{2*N - n}, ...
                            SL{2*N - n}, SpL{2*N - n}, mL{2*N - n}, ...
                                 HR{n}, QR{n}, SR{n}, SpR{n}, mR{n}, ...
                                  psiLR{2*N - n}, idx_qc{2*N - n}, [2*N - n, n], ...
                                        s, m_fin, dynm_fin, trun_tol);  
    % Update maximum bond dimension m_max
    if m_current > m_max
        m_max = m_current;
    end    
    toc   
    trunc_times = [trunc_times,  [2*N-n ; n; m_current; toc]];
    end


%%%%%% Sweep from left to right to the center ------------->   %%%%%%%%%%
    % Check if it is a final sweep
    if sweep == sweeps
        final_site = 2*N - 2;
         % Save UR MPS data
        MPS_fDMRG_file_name_save = strcat('MPS_Schwinger_OBC_N', ...
                num2str(4*N),'_x', num2str(x), '_me', num2str(me), ...
                '_shf', num2str(shift),'_th', num2str(alpha),...
                '_s', num2str(s), '_swps',num2str(sweeps), ...
                '_m', num2str(m_max), '_lev', num2str(lev),'.mat');
        save(MPS_fDMRG_file_name_save, 'UR','mR','mRb');        
    else
        final_site = N;
    end
    for n = 3:final_site
        disp(['----> Sweep ',num2str(sweep), ...
              ' from left to right ----> Chain size:  (', ...
              num2str(n),',', num2str(2*N - n),')' ])
        tic       
        [HL{n}, QL{n}, SL{n}, SpL{n}, mL{n}, TL{n}] = ...
          build_Schwinger_blocks_left(HLb{n - 1}, QLb{n-1}, SLb{n - 1}, ...
                       SpLb{n - 1}, mLb{n - 1}, n, N, s, x, mu, y, theta); 
        toc
        build_times = [build_times, [n ; 2*N - n; toc ] ];

        tic
%%%%%%% Find ground eigenstate and eigenenergy of the superblock %%%%%%%
        [psiLR{n}, e0, idx_qc{n}] = ...
             eigs_Schwinger_superblock_guess_ltor(HL{n}, QL{n}, SpL{n}, mL{n},...
                 HR{2*N - n}, SR{2*N - n}, SpR{2*N - n}, mR{2*N - n}, ...
                 UL{n - 1}, TL{n}, UR{2*N - n}, TR{2*N - n + 1}, {psiLR{n - 1}{1,:}},...
                 idx_qc{n - 1}, idx_qc{n}, mRb{2*N - n}, [n, 2*N - n], ...
                                s, x, y, theta, lvs, kryl, eig_tol, maxit);     
        disp(['Energies of the target states = '])
        E0 = (diag(e0) - mu*N*2)/x
        list_of_all_energies = horzcat(list_of_all_energies, [n; 2*N-n; E0]);
        toc
        eigs_times = [eigs_times, [n ; 2*N - n; toc ]];
%%%%%%% Truncate left and right blocks %%%%%%%%%  
       tic
       [HLb{n}, QLb{n}, SLb{n}, SpLb{n}, UL{n}, mLb{n}, ...
        HRb{2*N - n}, QRb{2*N - n}, SRb{2*N - n}, SpRb{2*N - n},...
                                    UR{2*N - n}, mRb{2*N - n}, m_current] = ...
                truncate_Schwinger_blocks_fDMRG(HL{n}, QL{n}, SL{n},...
                            SpL{n}, mL{n}, HR{2*N - n}, QR{2*N - n}, ...
                            SR{2*N - n}, SpR{2*N - n}, mR{2*N - n}, ...
                                    psiLR{n}, idx_qc{n}, [n, 2*N - n], ...
                                        s, m_fin, dynm_fin, trun_tol); 
       % Update maximum bond dimension m_max
       if m_current > m_max
           m_max = m_current;
       end    
       toc
       trunc_times = [trunc_times,  [n ; 2*N-n; m_current; toc ]];
    end 

    
% Add new energy to the list and stop timer and store the elapsed time    
list_of_energies = [list_of_energies, [toc; 2*n; E0]];


endTime = datetime('now');
elapsedTime = seconds(endTime - startTime)
sweep_times = [sweep_times, [sweep ; elapsedTime]];
%}
end 

if sweeps > 0 
    % Create separately left and right indices of quantum combinations  
    idx_qcL = idx_qc;
    for n = 2:(2*N - 2)
        idx_qcR{n} = idx_qc{2*N - n};
    end  
    % Save remaining MPS data after finite DMRG into the same MPS file
    save(MPS_fDMRG_file_name_save, 'UL','mL','mLb', 'psiLR', ...
                            'idx_qcL', 'idx_qcR', '-append'); 
end




%{
if sweeps > 0 
% Plot of times 
% Create figure and plot
figure;
yyaxis left
plot(build_times(3,:),'.-r')
hold on
plot(eigs_times(3,:),'.-b')
plot(trunc_times(4,:),'.-y')
yyaxis right
plot(trunc_times(3,:),'.-g')
hold off
end
%}

    
