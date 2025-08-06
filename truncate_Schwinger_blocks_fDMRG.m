% v1. 09/29/2024. 
% This function implements DMRG truncation for the Schwinger model 
% L and R  blocks during the finite DMRG procedure
% Input:
%       HL, QL, SL, SpL: left blocks
%       mL:  array of dimensions of the left blocks
%       HR, QR, SR, SpR: right blocks 
%       mR: array of dimensions of the right blocks
%       psiLR: cell of eigenfunctions in L-R form
%       N = [NL, NR]: number of the DMRG sites of the left and right blocks
%       s: spin projection of the superblock
%       m: bond dimensions
%       dynm: 0 or 1, of dynm == 1 the bond dimension m is chosen such that 
%              truncation error is below the truncation tolerance
%       trun_tol: truncation tolerance, needed only for for dynm = 1
% Output:
%       HLb, QLb, SLb, SpLb: truncated left blocks
%       UL: cell of left truncation matrices
%       mLb:  array of dimensions of the truncated left blocks  
%       HRb, QRb, SRb, SpRb : truncated right blocks
%       UR: cell of right truncation matrices
%       mRb:  array of dimensions of the truncated right blocks 
%       m_current:  current bond dimension
function [HLb, QLb, SLb, SpLb, UL, mLb, ...
          HRb, QRb, SRb, SpRb, UR, mRb, m_current] = ... 
               truncate_Schwinger_blocks_fDMRG(HL, QL, SL, SpL, mL,...
                                               HR, QR, SR, SpR, mR, ...   
                                   psiLR, idx_qc, N, s, m, dynm, trun_tol)
                                  
lvs = size(psiLR, 1);  % number of energy levels                                        

k = size(idx_qc, 1); % number of combinations of the L and R quantum numbers
% For each LR combination of indices find dimensions of L and R blocks  
mLc = mL(1, idx_qc(:, 1));
mRc = mR(1, idx_qc(:, 2));
mSBs = sum(mLc.*mRc);    % total dimension of the SuperBlock with spin s
sum_mLc = sum(mLc);   % dimension of all the left blocks of combinations 
sum_mRc = sum(mRc);   % dimension of all the right blocks of combinations

% Find density matrices for each level and each combination of qnt numbers   
rhoL = cell(lvs, k);
rhoR = cell(lvs, k);
for l = 1:lvs
    for j = 1:k
        rhoL{l, j} = psiLR{l, j}*psiLR{l, j}';
        rhoR{l, j} = transpose(psiLR{l, j}'*psiLR{l, j});
    end
end

% Find average density matrices for each qnt numbers combination
rhoLav = cell(k, 1);
rhoRav = cell(k, 1);
for j = 1:k
    rhoLav{j} = (1/lvs)*rhoL{1, j};
    rhoRav{j} = (1/lvs)*rhoR{1, j};
    for l = 2:lvs
        rhoLav{j} = rhoLav{j} + (1/lvs)*rhoL{l, j};
        rhoRav{j} = rhoRav{j} + (1/lvs)*rhoR{l, j};
    end
end

% Eigenstates and probability matrices
UL = cell(k, 1);
pL = cell(k, 1);
UR = cell(k, 1);
pR = cell(k, 1);
for j = 1:k
    [UL{j}, pL{j}] = eig(rhoLav{j});
    [UR{j}, pR{j}] = eig(rhoRav{j});
end


% Combine all the probabilities pL and pR and construct their index map  
all_pL = zeros(1, sum_mLc);
all_pR = zeros(1, sum_mRc);
pL_idx = zeros(sum_mLc, 2);
pR_idx = zeros(sum_mRc, 2);
% Count indices for L and R
cdxL = 1;
cdxR = 1;
for j = 1:k
    % Combine all probabilities from combinations j into a single array
    all_pL(1, cdxL : cdxL + mLc(j) - 1) = diag(pL{j});
    all_pR(1, cdxR : cdxR + mRc(j) - 1) = diag(pR{j});
    % Create an index map between combinations arrays and total array  
    pL_idx(cdxL : cdxL + mLc(j) - 1, :) = [j*ones(mLc(j), 1), (1:mLc(j))'];
    pR_idx(cdxR : cdxR + mRc(j) - 1, :) = [j*ones(mRc(j), 1), (1:mRc(j))'];
    cdxL = cdxL + mLc(j);
    cdxR = cdxR + mRc(j);    
end
% Sort all the probabilities pL and pR and get the sorting index 
[all_pL_sorted, idx_sort_pL] = sort(all_pL, 'descend');
[all_pR_sorted, idx_sort_pR] = sort(all_pR, 'descend');


% Choose proper bond dimension m 
if dynm == 1  
    min_mL = find(cumsum(all_pL_sorted) > 1 - trun_tol, 1);
    if isempty(min_mL)
        min_mL = size(all_pL_sorted, 2);
    end
    min_mR = find(cumsum(all_pR_sorted) > 1 - trun_tol, 1);
    if isempty(min_mR)
        min_mR = size(all_pR_sorted, 2);
    end
    % min_mL and min_mR should be equal to each other, but some small
    % numerical fluctuation is possible
    m = min(min_mL, min_mR);
elseif dynm == 0
   if sum_mLc < m || sum_mRc < m
        m = min(sum_mLc, sum_mRc);
   end
end

% Compute truncation error 
tr_error_pL = 1 - sum(all_pL_sorted(1:m));
tr_error_pR = 1 - sum(all_pR_sorted(1:m));
% Compute positions of largest probabilities inside each combination 
idx_max_pL =  sortrows(pL_idx(idx_sort_pL(1:m), :), 1);
idx_max_pR =  sortrows(pR_idx(idx_sort_pR(1:m), :), 1);   
% Find the counts of indices of largest probabilities in each combination j 
counts_max_pL = histcounts(idx_max_pL(:, 1), 1:(k+1));
counts_max_pR = histcounts(idx_max_pR(:, 1), 1:(k+1));
% Create cells of highest probability indices in each combination j
idx_max_pLc = cell(k, 1);
idx_max_pRc = cell(k, 1);
idx_max_pLc = mat2cell(idx_max_pL(:, 2), counts_max_pL);
idx_max_pRc = mat2cell(idx_max_pR(:, 2), counts_max_pR);


m_current = m;
disp(['Bond dimension = ', num2str(m_current)])  
disp(['truncation error left = ', num2str(tr_error_pL)])
disp(['truncation error right = ', num2str(tr_error_pR)])


% Truncate matrices UL and UR for each combination of quantum numbers
% and update new truncated dimensions of the combinations
for j = 1:k
    % Left blocks
    UL{j} = UL{j}(:, idx_max_pLc{j}');
    % New truncated dimensions of the left blocks quantum sectors
    mLc(j) = size(idx_max_pLc{j}, 1);
    % Right blocks
    UR{j} = UR{j}(:, idx_max_pRc{j}');
    % New truncated dimensions of the right blocks quantum sectors
    mRc(j) = size(idx_max_pRc{j}, 1);
end

% Create arrays of truncated dimensions mLb and mRb
% and include also indices iL, iR of mL and mR as 3rd row
mLb = [mLc; idx_qc(:, 3)'; idx_qc(:, 1)'];
[~, sort_qRc] = sort(idx_qc(:,4));
mRb = [mRc(sort_qRc); idx_qc(sort_qRc, 4)'; idx_qc(sort_qRc, 2)'];

% Leave only non-zero dimensions in the truncated dimensions arrays
mLb = mLb(:, mLb(1, :) ~= 0);
mRb = mRb(:, mRb(1, :) ~= 0);

% Create new truncated left blocks 
HLb = cell(size(mLb, 2));
QLb = cell(size(mLb, 2));
SLb = cell(size(mLb, 2));
SpLb = cell(size(mLb, 2));
HLb  = HL(mLb(3, :), mLb(3, :));
QLb  = QL(mLb(3, :), mLb(3, :));
SLb = SL(mLb(3, :), mLb(3, :));
SpLb = SpL(mLb(3, :), mLb(3, :));

% Create new truncated right blocks 
HRb = cell(size(mRb, 2));
QRb = cell(size(mRb, 2));
SRb = cell(size(mRb, 2));
SpRb = cell(size(mRb, 2));
HRb  = HR(mRb(3, :), mRb(3, :));
QRb  = QR(mRb(3, :), mRb(3, :));
SRb = SR(mRb(3, :), mRb(3, :));
SpRb = SpR(mRb(3, :), mRb(3, :));


% Truncate left Hamiltonian blocks
% for all truncated combinations of quantum numbers
for i1 = 1:size(mLb, 2)
    for i2 = 1:size(mLb, 2)
        % Get charges for these indices and combination indices jL1 and jL2
        jL1 = find(idx_qc(:, 1) == mLb(3, i1));
        jL2 = find(idx_qc(:, 1) == mLb(3, i2));
        qL1 = mLb(2, i1);
        qL2 = mLb(2, i2);
        if qL1 == qL2
            HLb{i1, i2} = UL{jL1}'*full(HLb{i1, i2})*UL{jL2};
            QLb{i1, i2} = UL{jL1}'*full(QLb{i1, i2})*UL{jL2};
            SLb{i1, i2} = UL{jL1}'*full(SLb{i1, i2})*UL{jL2};
        else
            HLb{i1, i2} = sparse(mLb(1, i1), mLb(1, i2));
            QLb{i1, i2} = sparse(mLb(1, i1), mLb(1, i2));
            SLb{i1, i2} = sparse(mLb(1, i1), mLb(1, i2));
        end
        if qL1 == qL2 + 1 
            SpLb{i1, i2} = UL{jL1}'*full(SpLb{i1, i2})*UL{jL2};
        else
            SpLb{i1, i2} = sparse(mLb(1, i1), mLb(1, i2));
        end
    end
end

% Truncate right Hamiltonian blocks  
% for all truncated combinations of quantum numbers
for i1 = 1:size(mRb, 2)
    for i2 = 1:size(mRb, 2)
        % Get charges for these indices and combination indices jR1 and jR2
        jR1 = find(idx_qc(:, 2) == mRb(3, i1));
        jR2 = find(idx_qc(:, 2) == mRb(3, i2));
        qR1 = mRb(2, i1);
        qR2 = mRb(2, i2);
        if qR1 == qR2
            HRb{i1, i2} = UR{jR1}'*full(HRb{i1, i2})*UR{jR2};
            QRb{i1, i2} = UR{jR1}'*full(QRb{i1, i2})*UR{jR2};
            SRb{i1, i2} = UR{jR1}'*full(SRb{i1, i2})*UR{jR2};
        else
            HRb{i1, i2} = sparse(mRb(1, i1), mRb(1, i2));
            QRb{i1, i2} = sparse(mRb(1, i1), mRb(1, i2));
            SRb{i1, i2} = sparse(mRb(1, i1), mRb(1, i2));
        end
        if qR1 == qR2 + 1 
            SpRb{i1, i2} = UR{jR1}'*full(SpRb{i1, i2})*UR{jR2};
        else
            SpRb{i1, i2} = sparse(mRb(1, i1), mRb(1, i2));
        end
    end
end


% Create new truncated list of dimensions mLbz and mRbz
mLb = mLb(1:2, :);
mRb = mRb(1:2, :);

end


