% Dempsey et all convention
% v1. 09/28/2024. 
% v2. 09/06/2024. (Bug with counts_of_q is fixed)
% This function builds left blocks in spin-spin projections 
% for the Schwinger model (Dempsey convention)
% Input:
%     HLb, QLb, SLb, SpLb: left block operators
%     mLb: quantum sectors dimensions of the initial left blocks 
%     n: Number of DMRG sites in the new built left blocks
%     N: Final number of DMRG sites in the left block at the end of iDMRG
%     s: Total projection of spin in the whole chain
%     x, mu, y, theta: Schwinger model parameters
% Output:
%     HL, QL, SL, SpL: new built left blocks
%     mL: quantum sectors dimensions of the new left blocks
%     TL: ledger of quantum numbers combinations and dimensions
function [ HL, QL, SL, SpL, mL, TL] = build_Schwinger_blocks_left(HLb, QLb,...
                                   SLb, SpLb, mLb, n, N, s, x, mu, y, theta)
% Single DMRG site quantum sector dimensions 
mS = [1, 2, 1; 1, 2, 3];
% Sum of all quantum sectors dimensions in the DMRG site: 
D = sum(mS(1, :));  
% DMRG site blocks' matrices (Dempsey convention)
H2_mat = sparse([theta^2/(4*pi^2),0,0,0; 0, theta^2/(4*pi^2),0,0; ...
                  0,0,1 + theta/pi + theta^2/(4*pi^2) ,0; ...
                  0,0,0, 1 + theta/pi + theta^2/(4*pi^2)]);
Q2_mat = sparse([-1,0,0,0; 0,0,0,0; 0,0,0,0; 0,0,0,1]);
S2_mat = sparse([-1,0,0,0; 0,-1,0,0; 0,0,1,0; 0,0,0,1]);
M2_mat = sparse([2,0,0,0; 0,0,0,0; 0,0,4,0; 0,0,0,2]);
HXY2_mat = sparse([0,0,0,0; 0,0,1,0; 0,1,0,0; 0,0,0,0]);
SpL2_mat = sparse([0,0,0,0; 1,0,0,0; 0,0,0,0; 0,0,1,0]);
SpR2_mat = sparse([0,0,0,0; 0,0,0,0; 1,0,0,0; 0,1,0,0]);
I4_mat = speye(D);

                                                                                   
% Sum of all non-zero sectors dimensions in the left blocks, i.e. 
% the bond dimension m:
m = sum(mLb(1, :));
                                                         
% Combine left blocks cells to matrices 
% Left blocks matrices
HLb_mat = sparse(cell2mat(HLb));
QLb_mat = sparse(cell2mat(QLb));
SLb_mat = sparse(cell2mat(SLb));
SpLb_mat = sparse(cell2mat(SpLb));
IdLb_mat = speye(m);
% Define auxiliarly operator  (Dempsey convention)
QLb_mat_Sq = QLb_mat*QLb_mat + (1/2 + theta/pi)*QLb_mat;


% Compute kronecker products of the blocks' matrices
HL_mat = kron(HLb_mat, I4_mat) + x*(kron(SpLb_mat, SpR2_mat') ...
              + kron(SpLb_mat', SpR2_mat)+ kron(IdLb_mat, HXY2_mat)) ...
         + (mu/2)*kron(IdLb_mat, M2_mat) + y*(kron(IdLb_mat, H2_mat) ... 
              + kron(QLb_mat, S2_mat) + 2*kron(QLb_mat_Sq, I4_mat) ...
              + (theta^2/(2*pi)^2)*kron(IdLb_mat, I4_mat));
QL_mat = kron(QLb_mat, I4_mat) + kron(IdLb_mat, Q2_mat);

SL_mat = kron(SLb_mat, I4_mat)+ kron(IdLb_mat, S2_mat) ...
                                + 4*kron(QLb_mat, I4_mat);
SpL_mat = kron(IdLb_mat, SpL2_mat);     

                                                            
% Create an array mL_all of all dimensions (including zeros)of blocks in HL
mLb_all = zeros(1, 2*n);
mLb_all(mLb(2, :)) = mLb(1, :);
mL_all = zeros(1, 2*n + 1);
mL_all = conv(mLb_all, mS(1, :));
% Remove zero elements from mL_all and obtain mL
qL_nz_mL = find(mL_all ~= 0); % get charges of non-zero elements in mL_all
mL = [mL_all(qL_nz_mL); qL_nz_mL];
length_mL = size(mL, 2);  % number of elements in mL


% Create an array  mLbmS of dimensions of blocks in part_kron(HLb, HS).
mLbmS = kron(mLb(1, :), mS(1, :));

% Blocks reordering in accordance to quantum sectors
% Generate all pairs (qL, qS), where qL = mLb(2, :) and qS = mS(2, :)=1,2,3
[QL, QS] = meshgrid(mLb(2, :), mS(2, :));
qLqS = [QL(:), QS(:)];
% Compute the sum (qL + qS) for each pair
sum_qLqS = sum(qLqS, 2);
% Get permutation of blocks by sorting the pairs based on the sum (qL + qS)
[sum_qLqS, perm] = sort(sum_qLqS);
% Get ordered pairs of quantum numbers qL and qS
qLqS = qLqS(perm, :);
% Construct permutation of indices 
idxP = zeros(sum(mLbmS), 1);
cdx = 1;    % current index
for i = 1:3*size(mLb, 2)
    idxP(cdx: cdx + mLbmS(perm(i)) - 1) = ...
        sum(mLbmS(1: perm(i) - 1)) + 1 : sum(mLbmS(1: perm(i))); 
    cdx = cdx + mLbmS(perm(i));
end

%%%%%%%%%% Permute matrices and compose cells for the left blocks %%%%%%%%%
KLS = K_mat(mS(1, :), mLb(1, :));
GLS = G_mat(mS(1, :), m);
PLS = speye(m*D);
PLS = PLS(idxP, :);
% Total permutation matrix
per_mat = PLS*KLS*GLS;
% Permutation of rows and columns
per = per_mat*(1:m*D)';

% Apply permutation to the left blocks matrices
HL_mat = HL_mat(per, per);
QL_mat = QL_mat(per, per);
SL_mat = SL_mat(per, per);
SpL_mat = SpL_mat(per, per);

% Pre-allocate new cells for the final left blocks
HL0 = cell(length_mL, length_mL);
QL0 = cell(length_mL, length_mL);
SL0 = cell(length_mL, length_mL);
SpL0 = cell(length_mL, length_mL);
% Represent matrices as cells for each pair of quantum numbers
HL0 = mat2cell(HL_mat, mL(1, :), mL(1, :));
QL0 = mat2cell(QL_mat, mL(1, :), mL(1, :));
SL0 = mat2cell(SL_mat, mL(1, :), mL(1, :));
SpL0 = mat2cell(SpL_mat, mL(1, :), mL(1, :));

%%%%%%%%%% Keep only needed quantum sectors in the left blocks %%%%%%%%%%%%
% Get needed spin domains 
[ sdn, ~ ] = spin_domain(2*n, 2*N, s);
% Compute new list of non-zero dimensions and their indices of blocks HL
% which are in the needed domain of quantum numbers sdn
idx_sdn = ismember(mL(2, :), sdn);
% New list of quantum sectors dimensions and their quantum numbers 
% within needed domain
mL = mL(:, idx_sdn);
length_mL = size(mL, 2);

% Pre-allocate cells for the finan left blocks within needed spins domain:
HL = cell(length_mL, length_mL);
QL = cell(length_mL, length_mL);
SL = cell(length_mL, length_mL);
SpL = cell(length_mL, length_mL);
% Construct new left blocks within needed spins domain
HL = HL0(idx_sdn, idx_sdn);
QL = QL0(idx_sdn, idx_sdn);
SL = SL0(idx_sdn, idx_sdn);
SpL = SpL0(idx_sdn, idx_sdn);


% Generate cells of ledgers TL 
% Get indices of quantum numbers qL from qLqS in the array mLb
[~, ind_qL] = ismember( qLqS(:, 1)', mLb(2, :));
% Get indices of quantum numbers qS  from qLqS in the array mS
[~, ind_qS] = ismember( qLqS(:, 2)', mS(2, :));
% Get dimensions mLc and mSc which are parts of combinations for ledger TL
mLbc = mLb(1, ind_qL);
mSc = mS(1, ind_qS);
% Create an array of (qL, qS, mLb, mS)
qL_qS_mLb_mS = [qLqS, mLbc', mSc'];
% Find logical index of (qL, qS) such that qL + qS - 1 is in sdn
[idx_q_in_sdn, ~] = ismember(sum_qLqS, sdn + 1);
% Trim the array qL_qS_mLb_mS excluding rows which have charges 
% sum of which is not in sdn
qL_qS_mLb_mS = qL_qS_mLb_mS(idx_q_in_sdn, :);
% Get number of combinations of qL and qS for each quantum number q in mL
[unique_vals, ~, idx_u_vals] = unique(sum_qLqS(idx_q_in_sdn), 'stable');
counts_of_q = accumarray(idx_u_vals, 1);
% Generate cells of ledgers TL
TL = cell(size(mL, 2), 1);
TL = mat2cell(qL_qS_mLb_mS, counts_of_q, 4);

end














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v1. 08/21/2024
% This script creates Kmn matrix:
% Input: mb = [m1, m2, ..., mr]  row dimensions of blocks 
%        nb = [n1, n2, ..., ns]  column dimensions of blocks 
function [ Kmn ] = K_mat(mb, nb)
m = sum(mb);
r = size(mb, 2); % length of row
n = sum(nb);
s = size(nb, 2); % length of row
% Create matrix Kmn 
idxK = zeros(m*n, 1);
for j = 1:s
    for i = 1:r
        idxK( m*sum(nb(1:j-1)) + sum(mb(1:i-1))*nb(j) + 1 : ...
              m*sum(nb(1:j-1)) + sum(mb(1:i))*nb(j)) =  ...
              sum(mb(1:i-1))*n + mb(i)*sum(nb(1:j-1)) + 1 : ...
              sum(mb(1:i-1))*n + mb(i)*sum(nb(1:j));
    end
end
Kmn = speye(m*n);
Kmn = Kmn(idxK, :);
end


% v1. 08/21/2024
% This script creates Gmn matrix:
% Input: mb = [m1, m2, ..., mr]  row dimensions of blocks 
%        n  = n1 + n2+ ... + ns  number of columns 
function [ Gmn ] = G_mat(mb, n)
m = sum(mb);
r = size(mb, 2); % length of row
% Create matrix Gmn 
idxG = zeros(n*m, 1);
IdG = reshape(1:m*n, m, n);
for i = 1:r
    idxG(sum(mb(1:i - 1))*n + 1: sum(mb(1:i))*n)  ...
            = reshape(IdG(sum(mb(1:i - 1)) + 1:sum(mb(1:i)), :), [], 1);
end
% Apply this permutation to the rows (i.e. to each column) of identity matrix
Gmn = speye(m*n);
Gmn = Gmn(idxG, :);
end



% v1. 08/28/2024
% This script creates domains of needed spin projections sdn,
% represented in the form of integer indices
% and also redundant spin projections sdr
% Input:
%       n: current number of spin-1/2 sites in the block
%       N: left or right total number of spin-1/2 sites in the chain
%       s: total spin of the whole chain of NL + NR sites
function [ sdn, sdr ] = spin_domain(n, N, s)
    if s >= 0
        if (n >= 2) && (n <= N - s)
            sdn = [-n/2: n/2] + n/2 + 1;
        elseif (n > N - s) && (n <= N + s)
            sdn = [n/2 - N + s: n/2] + n/2 + 1;
        elseif (n > N + s) && (n <= 2*N - 2)
            sdn = [n/2 - N + s : N + s - n/2] + n/2 + 1;
        end
    elseif s < 0
        if (n >= 2) && (n <= N - abs(s))
            sdn = [-n/2: n/2] + n/2 + 1;
        elseif (n > N - abs(s)) && (n <= N + abs(s))
            sdn = [-n/2: N - abs(s) - n/2] + n/2 + 1;
        elseif (n > N + abs(s)) && (n <= 2*N - 2)
            sdn = [n/2 - N - abs(s): N - abs(s) - n/2] + n/2 + 1;
        end    
    end
    sdr = setdiff([1: n + 1], sdn);
end
