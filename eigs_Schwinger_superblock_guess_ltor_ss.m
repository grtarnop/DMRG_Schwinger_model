% v1. 08/16/2024
% This function finds the ground eigenstate (in LR form) and energy of 
% the superblock for the Schwinger model in a given spin sector s
% and uses guessing algorithm for Left to Right step of the sweep
% Input: Left blocks: HLss, QLss, SLss, SpLss, IdLs
%        Right blocks: HRss, QRss, SRss, SpRss, IdRss
%        Number of chain sites of the left and right blocks: N = [N1, N2] 
%        Spin sector of the superblock: s
%        Couplings: x and y
%        Theta angle: theta
%        Number of levels: lvs
%        Dimension of the Krylov space: kryl
%        Accuracy of the eigenstates: eig_tol
%        Number of the sweep: sweep


function [ psiLR, e0 ] = eigs_Schwinger_superblock_guess_ltor_ss(HLss, QLss, SpLss, IdLss, UL, ...
                                     HRss, SRss, SpRss, IdRss, UR, TR, psi, ...
                                     N, s,  x, y, theta, lvs, kryl, eig_tol, sweep)
                                 
                                 
sct = spin_combinations(N, s);     % combinations of spin sectors
k = size(sct, 1);             % number of combinations of the LR spin sectors
% for each LR combination of the spin sectors find dimension of left and 
% right blocks  
for j = 1:k
    % dimensions of the left block for spin combination j
    dsl(j) = size(HLss{sct(j, 1) + N(1)/2 + 1, sct(j, 1) + N(1)/2 + 1}, 1);  
    % dimensions of the right block for spin combination j
    dsr(j) = size(HRss{sct(j, 2) + N(2)/2 + 1, sct(j, 2) + N(2)/2 + 1}, 1); 
end
ds = sum(dsl.*dsr);          % total dimension of the superblock sector s

%tic 
OPTS.disp = 0;   % Diagnostic information display level
if ds > kryl
    %OPTS.p = kryl;     % number of basis vectors 
end
OPTS.tol = eig_tol;  % Accuracy of the eigenstates

% Number of eigenenergies
if lvs > ds
    levels = ds;
else
    levels = lvs;      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Eigenstate guessing part: from left to right  %%%%%%%%%%%%%%%%%%%%
psi0 = psi;
k0 = size(psi, 2);  % number of spin combinations j0 = 1,...,k0 

% Step 1: Contraction of psi0 with UL
for j0= 1:k0
    psi0s1{j0} = UL{j0}'*psi0{j0};
end


% Step 2: Reshaping psi0s1
% Part a:
% all previous combinations which give total projection s
sct0 = spin_combinations([N(1) - 2,  N(2) + 2], s); 
for j0 = 1:k0
    p = 1;
    lrj = (N(2) + 2)/2 + 1 + sct0(j0, 2);    % lambda_{r,j}
    kij = size(TR{lrj}, 1);                  % k_{ij}
    for ij = 1:kij
        m1 = size(psi0s1{j0}, 1);         % m^{s_{l, j}}
        m2 = TR{lrj}(ij, 3);              % m^{sigma_{l+1, ij}}
        m3 = TR{lrj}(ij, 4);              % m^{lambda_{r-1,ij}}
        psi0s2a{j0}{ij} = permute(reshape(psi0s1{j0}(:, p:p + m2*m3 - 1), [m1, m3, m2]), [1, 3, 2]);
        p = p + m2*m3;
    end
end


% Part b:
for j = 1:k
    p = 1;
    for j0 = 1:k0
        lrj = (N(2) + 2)/2 + 1 + sct0(j0, 2);    % lambda_{r,j}
        kij = size(TR{lrj}, 1);                  % k_{ij}
        for ij = 1:kij
            m1 = size(psi0s2a{j0}{ij}, 1);    % m^{s_{l, j}}
            m2 = TR{lrj}(ij, 3);              % m^{sigma_{l+1, ij}}
            m3 = TR{lrj}(ij, 4);              % m^{lambda_{r-1,ij}}
            sgj = TR{lrj}(ij, 1);             % sigma_{l+1, ij}
            if sct0(j0, 1) + sgj == sct(j, 1)
                psi0s2b{j}(p:p + m1*m2 - 1,:) = reshape( permute(psi0s2a{j0}{ij}, [2, 1, 3]),[m1*m2, m3]);
                p = p + m1*m2;
            end
        end
    end
end


% Step 3: Contraction psi0s2b with UR
% adjust UR for the first quarter of the first sweep (some UR{j} may not exists before):
if sweep == 1 && (N(1) > N(2))
    sct_old = spin_combinations([N(2),  N(2)], s);
    if isempty(sct_old)
        for j = 1:k
            UR_new{j} = eye(dsr(j), size(psi0s2b{j}, 2));
        end
    else
        for j = 1:k
            j_old = find(sct_old(:, 2) == sct(j, 2));
            if isempty(j_old)
                UR_new{j} = zeros(dsr(j), size(psi0s2b{j}, 2));
            else
                UR_new{j} = UR{j_old};    
            end
        end
    end
UR = UR_new;    
end

for j = 1:k
    psi0s3{j} =  psi0s2b{j}*transpose(UR{j});
end


% represent the output vector psi0guess as a big column of size ds
psi0guess = zeros(ds, 1);
p = 0;
for j = 1:k
    psi0guess(p + 1: p + dsl(j)*dsr(j)) = ...
        reshape( permute(psi0s3{j}, [2, 1]), [dsl(j)*dsr(j), 1] );                              
    p = p + dsl(j)*dsr(j);
end
OPTS.v0 = psi0guess;   % eigenstate initial guess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% define a function handle which implements multiplication by the Hamiltonian 
Afun = @(psiIn)apply_Schwinger_superblock_ss(psiIn, HLss, QLss, SpLss, IdLss, ...
                                          HRss, SRss, SpRss, IdRss, ...
                                          N, s,  x, y, theta);
% compute eigenstates and eigenenergies
[psi0, e0] = eigs(Afun, ds, levels, 'SA',  OPTS);   
%toc

% represent the wave function(s) psi0 in LR form as a cells for each 
% energy level and for each combination of the spin projections   
psiLR = cell(levels, k);
for l = 1:levels
    p = 0;
    for j = 1:k
        psiLR{l, j} = permute( reshape( psi0(p + 1: p + dsl(j)*dsr(j), l), ...
                                    [dsr(j), dsl(j)]), [2, 1]);
        p = p + dsl(j)*dsr(j);
    end
end


end




% v1. 11/11/2023
% v2. 11/17/2023
% v3. 05/14/2024  (included theta angle)
% v4. 08/16/2024  (changes in spin_combinations convention)
% apply Schwinger superblock projected on a given sector s
% input: spin-spin left and right blocks for the Schwinger model
%        N = [N1, N2]: number of chain sites of the left and right blocks
%        s: superblock spin sector 
function psiOut = apply_Schwinger_superblock_ss(psiIn, HLss, QLss, SpLss, IdLss, ...
                                    HRss, SRss, SpRss, IdRss, N, s, x, y, theta)
sct = spin_combinations(N, s); % combinations of spin sectors
k = size(sct, 1);         % number of the combinations of the spin sectors
% for each LR combination of the spin sectors find dimension of left and 
% right blocks 
for j = 1:k
    % dimensions of the left block for spin combination j
    dsl(j) = size(HLss{sct(j, 1) + N(1)/2 + 1, sct(j, 1) + N(1)/2 + 1}, 1);
    % dimensions of the right block for spin combination j
    dsr(j) = size(HRss{sct(j, 2) + N(2)/2 + 1, sct(j, 2) + N(2)/2 + 1}, 1);
end
ds = sum(dsl.*dsr);   % total dimension of the superblock in the sector s

% To deal with the combinations of the spin sectors separately we create 
% k cells for psi, which are represented as dsl x dsr matrices
% we reshape psiIn into matrices for each combination of spin sectors j
psi = cell(k, 1);
p = 0;
for j = 1:k
    psi{j} = permute( reshape( psiIn(p + 1: p + dsl(j)*dsr(j)), ...
                                [dsr(j), dsl(j)]), [2, 1]);
    p = p + dsl(j)*dsr(j);
end

% Implement multiplication by the superblock
% Perform multiplication by left and right blocks
%tic
w = cell(k, 1);
for j1 = 1:k
    w{j1} = sparse(dsl(j1), dsr(j1)); % initialize empty sparse arrays
    for j2 = 1:k
        % positive interger indices for spins projections
        q1L = sct(j1, 1) + N(1)/2 + 1;  
        q2L = sct(j2, 1) + N(1)/2 + 1;
        q1R = sct(j1, 2) + N(2)/2 + 1;
        q2R = sct(j2, 2) + N(2)/2 + 1;
        w{j1} = w{j1} + HLss{q1L, q2L}*psi{j2}*transpose(IdRss{q1R, q2R});
        w{j1} = w{j1} + IdLss{q1L, q2L}*psi{j2}*transpose(HRss{q1R, q2R});     
        w{j1} = w{j1} + y*QLss{q1L, q2L}*psi{j2}*transpose(SRss{q1R, q2R}); 
        w{j1} = w{j1} + y*N(2)*(QLss{q1L,q2L}.^2 + (1/2 + theta/pi)*QLss{q1L,q2L})...
                        *psi{j2}*transpose(IdRss{q1R,q2R});
        w{j1} = w{j1} + y*(theta^2/(2*pi)^2)*IdLss{q1L, q2L}*psi{j2}*transpose(IdRss{q1R,q2R});
        w{j1} = w{j1} + x*SpLss{q1L, q2L}*psi{j2}*transpose(SpRss{q2R, q1R}');
        w{j1} = w{j1} + x*SpLss{q2L, q1L}'*psi{j2}*transpose(SpRss{q1R, q2R}); 

    end
end
%toc

% represent the output vector psiOut as a big column of size ds
psiOut = zeros(ds, 1);
p = 0;
for j = 1:k
    psiOut(p + 1: p + dsl(j)*dsr(j)) = ...
        reshape( permute(w{j}, [2, 1]), [dsl(j)*dsr(j), 1] );                              
    p = p + dsl(j)*dsr(j);
end


end
