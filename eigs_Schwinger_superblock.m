% v1. 09/29/2024. (Dempsey et all convention)
% This function finds the ground eigenstate (in LR form) and energy of 
% the superblock for the Schwinger model in a given total spin projection s
% using Dempsey et all convention
% Input: 
%        HL, QL, SL, SpL: left blocks
%        mL: array of non-zero dimensions of the left blocks
%        HR, QR, SR, SpR: right blocks 
%        mR: array of non-zero dimensions of the right blocks
%        N = [NL, NR]: number of DMRG sites of the left and right blocks
%        s: spin projection of the superblock
%        x, y, theta: Schwinger model parameters  
%        lvs: number of energy levels
%        kryl: dimension of the Krylov space
%        eig_tol: accuracy of the eigenvalues
%        maxit: maximum number of iterations of Lanczos
% Output:
%       psiLR: eigenfunctions in LR format
%       e0: diag matrix of eigenenergies 
%       idx_qc: list of indices of quantum combinations 
function [ psiLR, e0, idx_qc ] = eigs_Schwinger_superblock(HL, QL, ...
                                        SpL, mL, HR, SR, SpR, mR, ...
                                             N, s, x, y, theta, ...
                                               lvs, kryl, eig_tol, maxit)

% All combinations of L and R indices iL and iR which correspond to the qL
% and qR quantum numbers, which give total spin quantum sector q (qL + qR = q + 1)
idx_qc = [];
for iL = 1:size(mL, 2)   % grows in ascending order
    for iR = size(mR, 2):-1:1
        if mL(2, iL) + mR(2, iR) == sum(2*N)/2 + s + 2
            idx_qc = [idx_qc ; [iL, iR, mL(2, iL), mR(2, iR), mL(1, iL), mR(1, iR)]];
        end
    end
end
k = size(idx_qc, 1); % number of combinations of the L and R quantum numbers
% For each LR combination of indices find dimensions of L and R blocks  
mLc = mL(1, idx_qc(:, 1));
mRc = mR(1, idx_qc(:, 2));
mSBs = sum(mLc.*mRc);    % total dimension of the SuperBlock with spin s


OPTS.disp = 0;   % Diagnostic information display level
% Check that Krylov dimension is not larger than the superblock size mSBs
if kryl < mSBs
  % OPTS.p = kryl; % (Dimension of the Krylov spase) number of basis vectors
end   

levels = lvs;     % Basic number of energy levels
% Check that number of eigenvalues is less than the superblock size  mSBs
if lvs > mSBs
    levels = mSBs;
end

%OPTS.tol = eig_tol;   % Eigenvalue accuracy 
%OPTS.maxit = 1;      % Maximum number of iterations


% Define a function handle which implements multiplication by the Hamiltonian 
Afun = @(psiIn)apply_Schwinger_superblock(psiIn, HL, QL, SpL, mLc, ...
                             HR, SR, SpR, mRc, idx_qc, N, x, y, theta);
% compute eigenstates and eigenenergies
[psi0, e0] = eigs(Afun, mSBs, levels, 'SA', OPTS);   


% Represent the wave function(s) psi0 in LR form as a cells for each 
% energy level and for each combination of the spin projections   
psiLR = cell(levels, k);
for l = 1:levels
    p = 0;
    for j = 1:k
        psiLR{l, j} = permute(reshape(psi0(p + 1: p + mLc(j)*mRc(j), l), ...
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
% using Dempsey et all convention
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



