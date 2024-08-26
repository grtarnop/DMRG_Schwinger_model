% v6. 08/15/2024 
% finite DMRG with spin projection for the Schwinger model with OBC (Dempsey convention)
clear all
format long
x = 100;        % gauge coupling x = 1/(ea)^2
me = 0.333333;  % m/e
mu = 2*sqrt(x)*(me) -1/4*1  ;       % mass term mu = 2m/(e^2 a)
y = 1;        % parameter to regulate Coulomb energy term
alpha = 1;
theta = pi*alpha;   % theta angle

N = 75;        % final number of DMRG sites (double-sites) in the left and right blocks (Length of spin-1/2 chain 4*N )
m_inf = 50;    % number of states kept under truncation for infinite DMRG
m_fin = 50;    % number of states kept under truncation for finite DMRG

% if dynm = 1 then m is chosen such that trunc_error < trun_tol
% if dynm = 0 then m is fixed to the initial value
dynm_inf = 0;   % switch parameter for infinite DMRG part
dynm_fin = 1;   % switch parameter for finite DMRG part
trun_tol = 1e-14; % truncation error tolerance (used when m is dynamical)
eig_tol = 1e-10;    % eigenvalue accuracy tolerance 
kryl = 8;       % dimension of the Krylov space

lvs = 1;      % number of lowest energy levels
s = 0;        % Total spin projection Sz 
sweeps = 4;   % number of DMRG sweeps



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% First DMRG site preparation %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic single-site operators 
I2 = sparse([1 0; 0, 1]);      % Identity matrix
sp = sparse([0 1; 0 0]);       % S^+ matrix
sz = sparse([1 0; 0 -1]);      % Sz Pauli matrix
% basic double-site operators 
H2 = (1/2)*kron((1 + theta/pi)*(I2 + sz) + theta^2/(2*pi^2)*I2, I2);   % Hamiltonian operator HC for 2 sites           
Q2 = (1/2)*(kron(sz, I2) + kron(I2, sz));   % charge operator Q
S2 = kron(sz, I2);                          % operator S             
M2 = kron(I2 + sz, I2) + kron(I2, I2 - sz); % mass operator for 2 sites
HXY2 = kron(sp, sp') + kron(sp', sp);       % XY Hamiltonian for 2 sites
SpL2 = kron(I2, sp);
SpR2 = kron(sp, I2);


%%% Construct spin-spin sectors of the basic double-site operators %%%%%%
v = cell(3, 1);
for q = 1:3
    v{q} = subspace_fixed_Sz(2, -1 + q - 1);
end

I4ss = cell(3, 3);      % identity for two spin-1/2 sites
H2ss = cell(3, 3);      % H2 
Q2ss = cell(3, 3);      % Q2
S2ss = cell(3, 3);      % S2
M2ss = cell(3, 3);      % M2   
HXY2ss = cell(3, 3);    % HXY2
SpL2ss = cell(3, 3);    % SpL2
SpR2ss = cell(3, 3);    % SpR2
for q1 = 1:3            % runs over all 3 different spin sectors 
    for q2 = 1:3        % runs over all 3 different spin sectors 
        I4ss{q1, q2} = v{q1}'*kron(I2, I2)*v{q2};
        H2ss{q1, q2} = v{q1}'*H2*v{q2};
        Q2ss{q1, q2} = v{q1}'*Q2*v{q2};
        S2ss{q1, q2} = v{q1}'*S2*v{q2};
        M2ss{q1, q2} = v{q1}'*M2*v{q2};
        HXY2ss{q1, q2} = v{q1}'*HXY2*v{q2};
        SpL2ss{q1, q2} = v{q1}'*SpL2*v{q2};
        SpR2ss{q1, q2} = v{q1}'*SpR2*v{q2};
    end
end

% initialize a list of N empty matrices 
% for HL, HR,  QL, QR, SL, SR, SpL, SpR, IdL, IdR operators
% left blocks
HL = cell(N,1); QL = cell(N,1); SL = cell(N,1); 
SpL = cell(N,1);IdL = cell(N, 1);
% right blocks
HR = cell(N,1); QR = cell(N,1); SR = cell(N,1); 
SpR = cell(N,1); IdR = cell(N, 1);

% initialize a list of N empty cells 
% for HLss, SLss ...,HRss, SRss ... operators
% left blocks
HLss = cell(N, 1); QLss = cell(N, 1); SLss = cell(N, 1); 
SpLss = cell(N,1); IdLss = cell(N, 1); UL = cell(N, 1);
% right blocks
HRss = cell(N, 1); QRss = cell(N, 1); SRss = cell(N, 1); 
SpRss = cell(N,1); IdRss = cell(N, 1); UR = cell(N, 1);

% initialize a list of N empty cells for spin combinations 
% and space dimensions ledgers TL and TR
TL = cell(N, 1);
TR = cell(N, 1);

% Values for the initial double-site operators
% Left blocks
HL{1} = x*(kron(sp, sp') + kron(sp', sp)) ...
        + (mu/2)*(kron(I2 + sz, I2) + kron(I2, I2 - sz)) + y*H2;
QL{1} = Q2;
SL{1} = S2;
SpL{1} = kron(I2, sp);
IdL{1} = kron(I2, I2);
% Right blocks
HR{1} = x*(kron(sp, sp') + kron(sp', sp)) ... 
        + (mu/2)*(kron(I2 + sz, I2) + kron(I2, I2 - sz)) + y*H2;

QR{1} = Q2;
SR{1} = S2;
SpR{1} = kron(sp, I2);
IdR{1} = kron(I2, I2);


%%% Construct all Hamiltonians in spin-spin sectors for double-site %%%%%%%
nd = 1; % current double-site chain size
% Set of vectors with a given Sz value  
v = cell(2*nd + 1, 1);
for q = 1:(2*nd + 1)
    v{q} = subspace_fixed_Sz(2*nd, -nd + q - 1);
end
% Left blocks 
HLss{nd} = cell(2*nd + 1, 2*nd + 1); 
QLss{nd} = cell(2*nd + 1, 2*nd + 1); 
SLss{nd} = cell(2*nd + 1, 2*nd + 1); 
SpLss{nd} = cell(2*nd + 1, 2*nd + 1); 
IdLss{nd} = cell(2*nd + 1, 2*nd + 1);
% Right blocks
HRss{nd} = cell(2*nd + 1, 2*nd + 1); 
QRss{nd} = cell(2*nd + 1, 2*nd + 1); 
SRss{nd} = cell(2*nd + 1, 2*nd + 1);
SpRss{nd} = cell(2*nd + 1, 2*nd + 1); 
IdRss{nd} = cell(2*nd + 1, 2*nd + 1);

for q1 = 1:(2*nd + 1)       % runs over all 2*nd + 1 different spin sectors 
    for q2 = 1:(2*nd + 1)   % runs over all 2*nd + 1 different spin sectors 
        % Left blocks
        HLss{nd}{q1, q2} = v{q1}'*HL{nd}*v{q2}; 
        QLss{nd}{q1, q2} = v{q1}'*QL{nd}*v{q2};
        SLss{nd}{q1, q2} = v{q1}'*SL{nd}*v{q2}; 
        SpLss{nd}{q1, q2} = v{q1}'*SpL{nd}*v{q2};
        IdLss{nd}{q1, q2} = v{q1}'*IdL{nd}*v{q2};
        % Right blocks
        HRss{nd}{q1, q2} = v{q1}'*HR{nd}*v{q2}; 
        QRss{nd}{q1, q2} = v{q1}'*QR{nd}*v{q2};
        SRss{nd}{q1, q2} = v{q1}'*SR{nd}*v{q2}; 
        SpRss{nd}{q1, q2} = v{q1}'*SpR{nd}*v{q2};
        IdRss{nd}{q1, q2} = v{q1}'*IdR{nd}*v{q2};
    end
end





%tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Start infinite DMRG loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['****************************************'])
disp(['********* Infinite DMRG begins *********'])
disp(['****************************************'])

% Start a parallel pool if it doesn't exist
if isempty(gcp('nocreate'))
    parpool;
end

% List of times and all energies 
list_of_times = [];  

tic
for nd = 2:N

disp(['========================================'])    
disp(['Number of DMRG sites of the whole chain = ', num2str(2*nd)])    
disp(['========================================']) 

%%%%  Build all Hamiltonians in spin-spin sectors for nd double-sites %%%%%
[HLss{nd}, QLss{nd}, SLss{nd}, SpLss{nd}, IdLss{nd}, TL{nd}] = ...
    build_Schwinger_blocks_left_ss(HLss{nd - 1}, QLss{nd - 1}, SLss{nd - 1}, ...
              SpLss{nd - 1}, IdLss{nd - 1}, H2ss, Q2ss, S2ss, M2ss, ...
              HXY2ss, SpL2ss, SpR2ss, I4ss, 2*nd, x, mu, y, theta, s, 2*N); 

[HRss{nd}, QRss{nd}, SRss{nd}, SpRss{nd}, IdRss{nd}, TR{nd}] = ...
    build_Schwinger_blocks_right_ss(HRss{nd - 1}, QRss{nd - 1}, SRss{nd - 1}, ...
              SpRss{nd - 1}, IdRss{nd - 1}, H2ss, Q2ss, S2ss, M2ss, ...
              HXY2ss, SpL2ss, SpR2ss, I4ss, 2*nd, x, mu, y, theta, s, 2*N);
toc              
list_of_times = [list_of_times, [toc; nd; nd; 0]];

%%%% Find ground eigenstate and eigenenergy of the superblock %%%%%%%%%%%%%
    if 2*nd >= abs(s)
        tic
        [psiLR, e0] = eigs_Schwinger_superblock_ss(HLss{nd}, QLss{nd}, ...
                          SpLss{nd}, IdLss{nd}, HRss{nd}, SRss{nd},...
                          SpRss{nd}, IdRss{nd}, [2*nd, 2*nd], s, x, y, ...
                                             theta, lvs, kryl, eig_tol);
                                     
        disp(['Energies of the target states = '])
        E0 = (diag(e0) - mu*N*2)/x
        toc
        list_of_times = [list_of_times, [toc; nd; nd; E0]];

        %%%%%%%%%% Truncate left and right blocks %%%%%%%%%  
        tic
        [HLss{nd}, QLss{nd}, SLss{nd}, SpLss{nd}, IdLss{nd}, UL{nd}, ...
         HRss{nd}, QRss{nd}, SRss{nd}, SpRss{nd}, IdRss{nd}, UR{nd}] = ...
            truncate_Schwinger_blocks_ss(HLss{nd}, QLss{nd}, SLss{nd}, ...
                                         SpLss{nd}, IdLss{nd}, HRss{nd}, ...
                                         QRss{nd}, SRss{nd}, SpRss{nd}, ...
                                         IdRss{nd},psiLR, [2*nd, 2*nd], ...   
                                         s, 2*N, m_inf, dynm_inf, trun_tol);
        toc
        list_of_times = [list_of_times, [toc; nd; nd; E0]];
    end                          


end
%toc

% List of energies for infinite DMRG part and sweeps
list_of_energies = [];  
% Add new energy to the list and stop timer and store the elapsed time
list_of_energies = [list_of_energies, [2*nd; E0]];





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Start finite DMRG loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['**************************'])
disp(['*** Finite DMRG begins ***'])
disp(['**************************'])

for sweep = 1:sweeps
%tic 
disp(['Sweep number = ', num2str(sweep)])    
    %%%%%% Half sweep from left to right   ------------->   %%%%%%%%%%%%%%%
    for nd = (N + 1):(2*N - 2)
        disp(['----> Sweep ',num2str(sweep),' from left to right ----> Chain size:  (', ...
                              num2str(2*nd),',', num2str(4*N - 2*nd),')' ])
        tic
        [HLss{nd}, QLss{nd}, SLss{nd}, SpLss{nd}, IdLss{nd}, TL{nd}] = ...
          build_Schwinger_blocks_left_ss(HLss{nd - 1}, QLss{nd - 1}, ...
                                        SLss{nd - 1}, SpLss{nd - 1}, ...
                                        IdLss{nd - 1}, H2ss, Q2ss, S2ss, ...
                                        M2ss, HXY2ss, SpL2ss, SpR2ss, ...
                                        I4ss, 2*nd, x, mu, y, theta, s, 2*N);
        
        [HRss{2*N - nd}, QRss{2*N - nd}, SRss{2*N - nd}, ...
         SpRss{2*N - nd}, IdRss{2*N - nd}, TR{2*N - nd}] = ...
          build_Schwinger_blocks_right_ss(HRss{2*N - nd - 1}, QRss{2*N -nd - 1}, ...
                                        SRss{2*N - nd - 1}, SpRss{2*N - nd - 1},...
                                        IdRss{2*N - nd - 1}, H2ss, Q2ss, ...
                                        S2ss, M2ss, HXY2ss, SpL2ss, SpR2ss, ...
                                        I4ss, 4*N - 2*nd, x, mu, y, theta, s, 2*N);
        toc    
        list_of_times = [list_of_times, [toc; nd; 2*N - nd; E0]];
        %%%% Find ground eigenstate and eigenenergy of the superblock %%%%%%%%%
        tic
        [psiLR, e0] = eigs_Schwinger_superblock_guess_ltor_ss(HLss{nd}, ...
                            QLss{nd}, SpLss{nd}, IdLss{nd}, UL{nd - 1}, ...
                            HRss{2*N - nd}, SRss{2*N - nd}, SpRss{2*N - nd},...
                            IdRss{2*N - nd}, UR{2*N - nd}, TR{2*N - nd + 1}, ...
                            {psiLR{1,:}}, [2*nd, 4*N - 2*nd], s, x, y, theta, ...
                            lvs, kryl, eig_tol, sweep);
        
        disp('Lowest energies = ')
        E0 = (diag(e0) - mu*N*2)/x
        toc 
        list_of_times = [list_of_times, [toc; nd; 2*N - nd; E0]]; 

        %%%% Truncate left and right blocks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        tic
        [HLss{nd}, QLss{nd}, SLss{nd}, SpLss{nd}, IdLss{nd}, UL{nd}, ...
                    HRss{2*N - nd}, QRss{2*N - nd}, SRss{2*N - nd}, ...
                    SpRss{2*N - nd}, IdRss{2*N - nd}, UR{2*N - nd}] = ...
                            truncate_Schwinger_blocks_ss(HLss{nd}, QLss{nd}, ...
                                            SLss{nd}, SpLss{nd}, IdLss{nd}, ...
                                            HRss{2*N - nd}, QRss{2*N - nd}, ...
                                            SRss{2*N - nd}, SpRss{2*N - nd}, ...
                                            IdRss{2*N - nd}, psiLR, [2*nd, 4*N - 2*nd], ...   
                                            s, 2*N, m_fin, dynm_fin, trun_tol);
        toc
        list_of_times = [list_of_times, [toc; nd; 2*N - nd; E0]]; 
    end


    %%%%%%  Half sweep from right to left  <------------- %%%%%%%%%%%%%%%%%
    for nd = 3:(2*N - 2)
        disp(['<---- Sweep ', num2str(sweep),' from right to left <---- Chain size: (', ...
                                num2str(4*N - 2*nd),',', num2str(2*nd),')' ])
        tic 
        [HLss{2*N - nd}, QLss{2*N - nd}, SLss{2*N - nd}, ...
         SpLss{2*N - nd}, IdLss{2*N - nd}, TL{2*N - nd}] = ...
             build_Schwinger_blocks_left_ss(HLss{2*N - nd - 1}, QLss{2*N - nd - 1}, ...
                                SLss{2*N - nd - 1}, SpLss{2*N - nd - 1}, ... 
                                IdLss{2*N - nd - 1}, H2ss, Q2ss, S2ss, M2ss, ...
                                HXY2ss, SpL2ss, SpR2ss, I4ss, 4*N - 2*nd, ...
                                x, mu, y, theta,  s, 2*N);
        [HRss{nd}, QRss{nd}, SRss{nd}, SpRss{nd}, IdRss{nd}, TR{nd}] = ...
             build_Schwinger_blocks_right_ss(HRss{nd - 1}, QRss{nd - 1}, ...
                                            SRss{nd - 1}, SpRss{nd - 1}, ... 
                                            IdRss{nd - 1}, H2ss, Q2ss, S2ss, ...
                                            M2ss, HXY2ss, SpL2ss, SpR2ss, I4ss, ...
                                            2*nd, x, mu, y, theta, s, 2*N);
        toc
        list_of_times = [list_of_times, [toc; 2*N - nd; nd; E0]]; 

        %%%%% Find ground eigenstate and eigenenergy of the superblock %%%%
        tic
        [psiLR, e0] = eigs_Schwinger_superblock_guess_rtol_ss(HLss{2*N - nd},  ...
                           QLss{2*N - nd}, SpLss{2*N - nd}, IdLss{2*N - nd}, ...
                           UL{2*N - nd}, TL{2*N - nd + 1}, HRss{nd}, ...
                           SRss{nd}, SpRss{nd}, IdRss{nd}, UR{nd - 1}, ...
                           {psiLR{1,:}}, [4*N - 2*nd, 2*nd], s, x, y, ...
                           theta, lvs, kryl, eig_tol, sweep);
        disp('Lowest energies = ')
        E0 = (diag(e0) - mu*N*2)/x
        toc
        list_of_times = [list_of_times, [toc; 2*N - nd; nd; E0]]; 

      %%%%%%%%%% Truncate left and right blocks %%%%%%%%%%%%%%%%%%%%%%%%%%%    
      tic
      [HLss{2*N - nd}, QLss{2*N - nd}, SLss{2*N - nd}, SpLss{2*N - nd}, ...
       IdLss{2*N - nd}, UL{2*N - nd}, HRss{nd}, QRss{nd}, SRss{nd},  ...
       SpRss{nd}, IdRss{nd}, UR{nd}] = ...
          truncate_Schwinger_blocks_ss(HLss{2*N - nd}, QLss{2*N - nd}, ...
               SLss{2*N - nd}, SpLss{2*N - nd}, IdLss{2*N - nd}, ...
               HRss{nd}, QRss{nd}, SRss{nd}, SpRss{nd}, IdRss{nd}, ...   
               psiLR, [4*N - 2*nd, 2*nd], s, 2*N, m_fin, dynm_fin, trun_tol);
      toc
      list_of_times = [list_of_times, [toc; 2*N - nd; nd; E0]]; 
    end
  
    %%%%%% Sweep from left to right to the center ------------->   %%%%%%%
    for nd = 3:N
    disp(['----> Sweep ',num2str(sweep),' from left to right ----> Chain size:  (', ...
                            num2str(2*nd),',', num2str(4*N - 2*nd),')' ])
    tic  
        [HLss{nd}, QLss{nd}, SLss{nd}, SpLss{nd}, IdLss{nd}, TL{nd}] = ...
          build_Schwinger_blocks_left_ss(HLss{nd - 1}, QLss{nd - 1}, ...
                                        SLss{nd - 1}, SpLss{nd - 1}, ...
                                        IdLss{nd - 1}, H2ss, Q2ss, S2ss, ...
                                        M2ss, HXY2ss, SpL2ss, SpR2ss, ...
                                        I4ss, 2*nd, x, mu, y, theta, s, 2*N);
        
        [HRss{2*N - nd}, QRss{2*N - nd}, SRss{2*N - nd}, ...
         SpRss{2*N - nd}, IdRss{2*N - nd}, TR{2*N - nd}] = ...
          build_Schwinger_blocks_right_ss(HRss{2*N - nd - 1}, QRss{2*N -nd - 1}, ...
                                        SRss{2*N - nd - 1}, SpRss{2*N - nd - 1},...
                                        IdRss{2*N - nd - 1}, H2ss, Q2ss, ...
                                        S2ss, M2ss, HXY2ss, SpL2ss, SpR2ss, ...
                                        I4ss, 4*N - 2*nd, x, mu, y, theta, s, 2*N);
    toc
    list_of_times = [list_of_times, [toc; nd; 2*N - nd; E0]]; 
    %%%%%% Find ground eigenstate and eigenenergy of the superblock %%%%%%%
        tic
        [psiLR, e0] = eigs_Schwinger_superblock_guess_ltor_ss(HLss{nd}, ...
                            QLss{nd}, SpLss{nd}, IdLss{nd}, UL{nd - 1}, ...
                            HRss{2*N - nd}, SRss{2*N - nd}, SpRss{2*N - nd},...
                            IdRss{2*N - nd}, UR{2*N - nd}, TR{2*N - nd + 1}, ...
                            {psiLR{1,:}}, [2*nd, 4*N - 2*nd], s, x, y, theta, ...
                            lvs, kryl, eig_tol, sweep);
        
        disp('Lowest energies = ')
        E0 = (diag(e0) - mu*N*2)/x
        toc 
        list_of_times = [list_of_times, [toc; nd; 2*N - nd; E0]]; 
     

    %%%% Truncate left and right blocks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    tic
    [HLss{nd}, QLss{nd}, SLss{nd}, SpLss{nd}, IdLss{nd}, UL{nd}, ...
                HRss{2*N - nd}, QRss{2*N - nd}, SRss{2*N - nd}, ...
                SpRss{2*N - nd}, IdRss{2*N - nd}, UR{2*N - nd}] = ...
                        truncate_Schwinger_blocks_ss(HLss{nd}, QLss{nd}, ...
                                        SLss{nd}, SpLss{nd}, IdLss{nd}, ...
                                        HRss{2*N - nd}, QRss{2*N - nd}, ...
                                        SRss{2*N - nd}, SpRss{2*N - nd}, ...
                                        IdRss{2*N - nd}, psiLR, [2*nd, 4*N - 2*nd], ...   
                                        s, 2*N, m_fin, dynm_fin, trun_tol);
    toc
    list_of_times = [list_of_times, [toc; nd; 2*N - nd; E0]]; 
    end   

%toc
% add new energy to the list and stop timer and store the elapsed time    
list_of_energies = [list_of_energies, [2*nd; diag(E0)]];
end 






