% v1. 09/29/2024.
% This function implements DMRG truncation for the Schwinger model 
% L and R  blocks during the infinite DMRG procedure
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
%       HRb, QRb, SRb, SpRb: truncated right blocks
%       UR: cell of right truncation matrices
%       mRb:  array of dimensions of the truncated right blocks 
%       m_current: current bond dimension
function [HLb, QLb, SLb, SpLb, UL, mLb, ...
          HRb, QRb, SRb, SpRb, UR, mRb, m_current] = ... 
                    truncate_Schwinger_blocks_iDMRG(HL, QL, SL, SpL, mL,...
                                       HR, QR, SR, SpR, mR, ...   
                                       psiLR, N, s, m, dynm, trun_tol)                                  
lvs = size(psiLR, 1);  % number of energy levels                                        

% All combinations of L and R indices iL and iR which correspond to the qL
% and qR quantum numbers, which give total spin s (qL + qR = q + 1)
idx_qc = [];
for iL = 1:size(mL,2)
    for iR = 1:size(mR, 2)
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
    min_mR = find(cumsum(all_pR_sorted) > 1 - trun_tol, 1);
    m = max(min_mL, min_mR);
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
% and create new truncated dimensions of the quantum sectors
for j = 1:k
    % Left blocks
    UL{j} = UL{j}(:, idx_max_pLc{j}');
    % New truncated dimensions of the left blocks quantum sectors
    mLc(j) = size(idx_max_pLc{j}, 1);
    % Right blocks
    UR{j} = UR{j}(:, idx_max_pRc{j});
    % New truncated dimensions of the right blocks quantum sectors
    mRc(j) = size(idx_max_pRc{j}, 1);
end


% Leave in mLb and mRb only quantum numbers which are in mLc and mRc and
% are in the non-particpating arrays qL_np and qR_np.
% Arrays of non-participating quantum numbers are needed in future
% steps, but can not partitipate in the quantum combinations at this step
% Aslo include indices iL, iR of mL and mR as 3rd row in mLb and mRb
if s >= 0
    qL_np = [1:s];
    qR_np = [1:s];
    idx_qL_np = ismember(mL(2, :), qL_np);
    idx_qR_np = ismember(mR(2, :), qR_np);
    mLb = [ [ mL(:, idx_qL_np ); find(idx_qL_np)], ...
            [mLc; idx_qc(:, 3)'; idx_qc(:, 1)']];
    mRb = [ [ mR(:, idx_qR_np ); find(idx_qL_np)], ...
            [mRc; idx_qc(:, 4)'; idx_qc(:, 2)']];
    % sort mRb array such that the quantum numbers are in ascending order
    [~, qR_order] = sort(mRb(2, : ));
    mRb = mRb(:, qR_order);    
elseif s < 0 
    qL_np = [2*N(1) + 2 + s : 2*N(1) + 1]; 
    qR_np = [2*N(2) + 2 + s : 2*N(2) + 1]; 
    idx_qL_np = ismember(mL(2, :), qL_np);
    idx_qR_np = ismember(mR(2, :), qR_np);
    mLb = [[mLc; idx_qc(:, 3)' ; idx_qc(:, 1)'], ...
            [mL(:, idx_qL_np); find(idx_qL_np)]];
    mRb = [[mRc; idx_qc(:, 4)' ; idx_qc(:, 2)'], ...
            [mR(:, idx_qR_np); find(idx_qR_np)]];
    % sort mRb array such that the quantum numbers are in ascending order
    [~, qR_order] = sort(mRb(2, : ));
    mRb = mRb(:, qR_order);
end

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


% Truncate left and right Hamiltonian blocks
% for given combinations of quantum numbers
for j = 1:k
    % Check if size of the left blocks for combination j is not zero
    if size(idx_max_pLc{j}, 1) ~= 0
        % Get an index of the new truncated left blocks for combination j
        iLj = find(mLb(3, :) == idx_qc(j, 1)); 
        % Get a q.n. of the new truncated left blocks for combination j
        qLj = idx_qc(j, 3);
        % Go over all indices of the new truncated left blocks
        for iL = 1:size(mLb, 2)
            qL = mLb(2, iL);
            if iL == iLj
                HLb{iL, iLj} = UL{j}'*full(HLb{iL, iLj})*UL{j};
                QLb{iL, iLj} = UL{j}'*full(QLb{iL, iLj})*UL{j};
                SLb{iL, iLj} = UL{j}'*full(SLb{iL, iLj})*UL{j};
            else
                HLb{iL, iLj} = sparse(mLb(1, iL), mLb(1, iLj));
                HLb{iLj, iL} = sparse(mLb(1, iLj), mLb(1, iL));
                QLb{iL, iLj} = sparse(mLb(1, iL), mLb(1, iLj));
                QLb{iLj, iL} = sparse(mLb(1, iLj), mLb(1, iL));
                SLb{iL, iLj} = sparse(mLb(1, iL), mLb(1, iLj));
                SLb{iLj, iL} = sparse(mLb(1, iLj), mLb(1, iL));
            end
            
            if qL == qLj + 1 
                SpLb{iL, iLj} = full(SpLb{iL, iLj})*UL{j};
            else
                SpLb{iL, iLj} = sparse(mLb(1, iL), mLb(1, iLj));
            end
            if qL == qLj - 1             
                SpLb{iLj, iL} = UL{j}'*full(SpLb{iLj, iL});
            else
                SpLb{iLj, iL}  = sparse(mLb(1, iLj), mLb(1, iL));
            end        
        end
    end
    % Check if size of the right blocks for combination j is not zero
    if size(idx_max_pRc{j}, 1) ~= 0
        % Get an index of the new truncated right blocks for combination j
        iRj = find(mRb(3, :) == idx_qc(j, 2)); 
        % Get a q.n. of the new truncated left blocks for combination j
        qRj = idx_qc(j, 4);
        % Go over all indices of the new truncated right blocks
        for iR = 1:size(mRb, 2)
            qR = mRb(2, iR);
            if iR == iRj
                HRb{iR, iRj} = UR{j}'*full(HRb{iR, iRj})*UR{j}; 
                QRb{iR, iRj} = UR{j}'*full(QRb{iR, iRj})*UR{j};  
                SRb{iR, iRj} = UR{j}'*full(SRb{iR, iRj})*UR{j}; 
            else
                HRb{iR, iRj} = sparse(mRb(1, iR), mRb(1, iRj));
                HRb{iRj, iR} = sparse(mRb(1, iRj), mRb(1, iR));
                QRb{iR, iRj} = sparse(mRb(1, iR), mRb(1, iRj));
                QRb{iRj, iR} = sparse(mRb(1, iRj), mRb(1, iR));
                SRb{iR, iRj} = sparse(mRb(1, iR), mRb(1, iRj));
                SRb{iRj, iR} = sparse(mRb(1, iRj), mRb(1, iR));
            end    
    
            if qR == qRj + 1 
                SpRb{iR, iRj} = full(SpRb{iR, iRj})*UR{j};
            else
                SpRb{iR, iRj} = sparse(mRb(1, iR), mRb(1, iRj));
            end
    
            if qR == qRj - 1             
                SpRb{iRj, iR} = UR{j}'*full(SpRb{iRj, iR});
            else
                SpRb{iRj, iR} = sparse(mRb(1, iRj), mRb(1, iR));
            end
                                     
        end
    end
end

% Create new truncated list of dimensions mLbz and mRbz
mLb = mLb(1:2, :);
mRb = mRb(1:2, :);

end





