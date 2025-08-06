% v1. 09/29/2024. (Dempsey et all convention)
% v2. 10/06/2024. (improved reshaping algorithm)
% This function finds the ground eigenstate (in LR form) and energy of 
% the superblock for the Schwinger model in a given total spin projection s
% and uses guessing algorithm for Left to Right step of the sweep
% ans uses Dempsey et all convention
% Input: 
%        HL, QL, SpL:   left blocks
%        mL: dimensions of left blocks 
%        HR, SR, SpR:  right blocks
%        mR: dimensions of right blocks 
%        UL(NL - 1) and UR(NR): unitary matrices for NL - 1 and NR
%        TL(NL), TR(NR + 1): spins and dimensions ledgers 
%        N = [NL, NR]: number of DMRG sites of the left and right blocks 
%        psi0: previous wave function in LR form for [NL - 1, NR + 1] chain 
%        idx_qc0: indies of quantum combinations for psi0
%        idx_qc_old: indies of quantum combinations for psi_old
%        mRb: truncated dimensions of right blocks for psi_old
%        s: spin projection of the superblock
%        x, y, theta: Schwinger model parameters  
%        lvs: number of energy levels
%        kryl: dimension of the Krylov space
%        eig_tol: accuracy of the eigenvalue 
%        maxit: maximum number of iterations
% Output:
%       psiLR: eigenfunctions in LR format
%       e0: diag matrix of eigenenergies 
%       idx_qc: list of indices of quantum combinations 
function [ psiLR, e0, idx_qc ] = ...
              eigs_Schwinger_superblock_guess_ltor(HL, QL, SpL, mL, ...
                                     HR, SR, SpR, mR, UL, TL, UR, TR, ...
                                     psi0, idx_qc0, idx_qc_old, mRb, ...
                                            N, s, x, y, theta, ...
                                            lvs, kryl, eig_tol, maxit)
% All combinations of L and R indices iL and iR which correspond to the qL
% and qR quantum numbers, which give total quantum sector q (qL + qR = q + 1)
idx_qc = [];
for iL = 1:size(mL, 2)   % grows in ascending order
    for iR = size(mR, 2):-1:1   % grows in descending order
        if mL(2, iL) + mR(2, iR) == sum(2*N)/2 + s + 2
            idx_qc = [idx_qc ; [iL, iR, mL(2, iL), mR(2, iR), mL(1, iL), mR(1, iR)]];
        end
    end
end
k = size(idx_qc, 1);  % number of combinations of the L and R quantum numbers
% For each LR combination of indices get dimensions of L and R blocks  
mLc = mL(1, idx_qc(:, 1));
mRc = mR(1, idx_qc(:, 2));
mSBs = sum(mLc.*mRc);    % total dimension of the SuperBlock with spin s

levels = 1;      % Basic number of energy levels
if mSBs > 1
    % check that Krylov dimension is not larger than the superblock size ds
    if kryl < mSBs
        OPTS.p = kryl;   % (Dimension of the Krylov spase) number of basis vectors
    end   
    % check that number of eigenvalues is smaller than the superblock size mSBs
    if lvs < mSBs
        levels = lvs;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Eigenstate guessing part: from left to right  %%%%%%%%%%%%%%%%%%%%
% Get number of the previous quantum combinations j0 = 1,...,k0 
k0 = size(idx_qc0, 1);  

% Step 1: Contraction of psi0 with UL(NL - 1)
psi0p = cell(1, k0);
for j0 = 1:k0
    psi0p{j0} = UL{j0}'*psi0{j0};
end


% Step 2: Reshaping psi0p into psi0p2
psi0p2 = cell(1, k);
for j = 1:size(idx_qc, 1)
    i = idx_qc(j, 1);
    q_rm1j = idx_qc(j, 4);
    p = 1;
    for ij = 1:size(TL{i}, 1)        
        q_lij = TL{i}(ij, 1);
        mLb_lij = TL{i}(ij, 3);
        mS_ij = TL{i}(ij, 4);
        j0 = find(idx_qc0(:, 3) == q_lij);
        q_rj0 = idx_qc0(j0, 4);
        i0 = idx_qc0(j0, 2);
        ij0 = find(TR{i0}(:, 2) == q_rm1j);
        if ~isempty(ij0)
            mRb_rm1j = TR{i0}(ij0, 4);
            p0 = 1 + TR{i0}(1:(ij0-1), 3)'*TR{i0}(1:(ij0-1), 4);
            psi0p2{j}(p:p + mLb_lij*mS_ij - 1, :) = ...
                  permute(reshape( permute(...
                    psi0p{j0}(:, p0:p0 + mS_ij*mRb_rm1j - 1), [2, 1]), ...
                    [mRb_rm1j, mLb_lij*mS_ij] ), [2, 1]);
            p = p + mLb_lij*mS_ij;
        else
            psi0p2{j}(p:p + mLb_lij*mS_ij - 1,:) = zeros(mLb_lij*mS_ij, 0);
            p = p + mLb_lij*mS_ij;
        end
    end
end


% Step 3: Contraction psi0p2 with UR
% Some quantum combinations labeled by j may not even exist before
% and so some UR{j}
if isempty(idx_qc_old)
    disp('There are no previous combinations thus UR')
    for j = 1:k
        % Check that psi0p2{j) exists 
        if ~isempty(psi0p2{j})
            % Multiply by UR^T.
            psi0p3{j} =  psi0p2{j}*speye(size(psi0p2{j}, 2), mRc(j));
        else
            psi0p3{j} =  sparse(mLc(j), mRc(j));
        end
    end  
else
    for j = 1:k
        if ~isempty(psi0p2{j})
            % Find old index j_old of q.c. which corresponds to the new one j 
            % by comparing quantum charges
            j_old = find(idx_qc_old(:, 4) == idx_qc(j, 4));
            if isempty(j_old)
                %disp('no j_old')
                psi0p3{j} =  psi0p2{j}*speye(size(psi0p2{j}, 2), mRc(j));
            else
                % Multiply by UR^T.
                psi0p3{j} =  psi0p2{j}*transpose(UR{j_old});
            end
        else
            psi0p3{j} =  zeros(mLc(j), mRc(j));
        end
    end
end




% represent the output vector psi0guess as a big column of size ds
psi0guess = zeros(mSBs, 1);
p = 0;
for j = 1:k
    psi0guess(p + 1: p + mLc(j)*mRc(j)) = ...
        reshape( permute(psi0p3{j}, [2, 1]), [mLc(j)*mRc(j), 1] );                              
    p = p + mLc(j)*mRc(j);
end
OPTS.v0 = psi0guess;   % eigenstate initial guess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}


% check that Krylov dimension is not larger than the superblock size ds
if kryl < mSBs
   OPTS.p = kryl;   % (Dimension of the Krylov spase) number of basis vectors
end   

OPTS.disp = 0;   % Diagnostic information display level
OPTS.tol = eig_tol;   % Eigenvalue accuracy 
OPTS.maxit = maxit;        % Maximum number of iterations



% Define a function handle which implements multiplication by the Hamiltonian 
Afun = @(psiIn)apply_Schwinger_superblock(psiIn, HL, QL, SpL, mLc, ...
                             HR, SR, SpR, mRc, idx_qc, N, x, y, theta);
% Compute eigenstates and eigenenergies
[psi0, e0] = eigs(Afun, mSBs, levels, 'SA', OPTS);   
%toc

% Represent the wave function(s) psi0 in LR form as a cells for each 
% energy level and for each combination of the spin projections   
psiLR = cell(levels, k);
for l = 1:levels
    p = 0;
    for j = 1:k
        psiLR{l, j} = permute( reshape( psi0(p + 1: p + mLc(j)*mRc(j), l), ...
                                    [mRc(j), mLc(j)]), [2, 1]);
        p = p + mLc(j)*mRc(j);
    end
end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v1. 09/29/2024. 
% apply Schwinger superblock projected on a given total spin projection s
% it uses Dempsey et all convention
% Input: 
%     psiIn: Initial wave vector in form a big column
%     HL, QL, SpL: left blocks 
%     mLc:  array of dimensions of L blocks needed for total spin s
%     HR, SR, SpR: right blocks 
%     mRc:  array of dimensions of R blocks needed for total spin s
%     idx_qc: combinations of L and R indices and quantum numbers for spin s
%     N = [NL, NR]: number of DMRG sites of the left and right blocks
%     x, y, theta: Schwinger model parameters  
% Output: 
%        psiOut: wave vector obtained from psiIn in form a big column        
function psiOut = apply_Schwinger_superblock(psiIn, HL, QL, SpL, mLc,...
                                            HR, SR, SpR, mRc, idx_qc,...
                                            N, x, y, theta)
                              
k = size(idx_qc, 1);  % number of combinations of the L and R quantum numbers
mSBs = sum(mLc.*mRc);  % total dimension of the SuperBlock with spin s

% To deal with the combinations of the spins separately we create 
% k cells for psiLR, which are represented as mLc(j) x mRc(j) matrices
% Reshape psiIn into matrices for each combination j of quantum numbers
psiLR = cell(k, 1);
p = 0;
for j = 1:k
    psiLR{j} = permute( reshape( psiIn(p + 1: p + mLc(j)*mRc(j)), ...
                                [mRc(j), mLc(j)]), [2, 1]);
    p = p + mLc(j)*mRc(j);
end

% Implement multiplication by the superblock
% Perform multiplication by left and right blocks (Dempsey convention)
w = cell(k, 1);
for j1 = 1:k
    w{j1} = sparse(mLc(j1), mRc(j1)); % initialize empty sparse arrays
    for j2 = 1:k
        % Left and Right indices of the combinations
        iL1 = idx_qc(j1, 1); 
        iR1 = idx_qc(j1, 2);       
        iL2 = idx_qc(j2, 1); 
        iR2 = idx_qc(j2, 2); 
        % Left and Right quantum numbers of the combinations
        qL1 = idx_qc(j1, 3);
        qR1 = idx_qc(j1, 4);
        qL2 = idx_qc(j2, 3);
        qR2 = idx_qc(j2, 4);

        if (iL1 == iL2) && (iR1 == iR2)
        w{j1} = w{j1} + full(HL{iL1, iL2})*psiLR{j2};
        w{j1} = w{j1} + psiLR{j2}*full(transpose(HR{iR1, iR2}));
        w{j1} = w{j1} + y*QL{iL1, iL2}*psiLR{j2}*transpose(SR{iR1, iR2}); 
        w{j1} = w{j1} + y*2*N(2)*(QL{iL1, iL2}.^2 + ...
                        (1/2 + theta/pi)*QL{iL1,iL2})*psiLR{j2};
        w{j1} = w{j1} + y*(theta^2/(2*pi)^2)*psiLR{j2};
        end
        
        if (qL1 - qL2 == 1) && (qR2 - qR1 == 1)
        w{j1} = w{j1} + x*full(SpL{iL1, iL2})*psiLR{j2}...
                            *full(transpose(SpR{iR2, iR1}'));
        end

        if (qL2 - qL1 == 1) && (qR1 - qR2 == 1)
        w{j1} = w{j1} + x*full(SpL{iL2, iL1}')*psiLR{j2}...
                            *full(transpose(SpR{iR1, iR2}));   
        end

    end
end

% Represent the output vector psiOut as a big column of size mSBs
psiOut = zeros(mSBs, 1);
p = 0;
for j = 1:k
    psiOut(p + 1: p + mLc(j)*mRc(j)) = ...
        reshape( permute(w{j}, [2, 1]), [mLc(j)*mRc(j), 1] );                              
    p = p + mLc(j)*mRc(j);
end

end





