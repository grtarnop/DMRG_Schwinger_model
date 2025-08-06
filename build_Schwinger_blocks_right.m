% Dempsey et all convention
% v1. 09/28/2024. 
% v2. 09/06/2024. (Bug with counts_of_q is fixed)
% This function builds right blocks in spin-spin projections 
% for the Schwinger model (Dempsey convention)
% Input:
%     HRb, QRb, SRb, SpRb: right block operators
%     mRb: quantum sectors dimensions of the initial right blocks 
%     n: Number of sites in the new built right blocks
%     N: Final number of sites in the left block at the end of iDMRG
%     s: Total projection of spin in the whole chain
%     x, mu, y, theta: Schwinger model parameters
% Output:
%     HR, QR, SR, SpR: new built right blocks
%     mR: quantum sectors dimensions of the new right blocks
%     TR: ledger of quantum numbers combinations and dimensions
function [ HR, QR, SR, SpR, mR, TR] = build_Schwinger_blocks_right(HRb, QRb,...
                                 SRb, SpRb, mRb, n, N, s, x, mu, y, theta)
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
% Define auxiliarly operator (Dempsey convention)
Q2_mat_Sq = Q2_mat*Q2_mat + (1/2 + theta/pi)*Q2_mat;
                                                                                      
% Sum of all non-zero sectors dimensions in the left blocks, i.e. 
% the bond dimension m:
m = sum(mRb(1, :));
                                                         
% Combine left blocks cells to matrices 
% Left blocks matrices
HRb_mat = sparse(cell2mat(HRb));
QRb_mat = sparse(cell2mat(QRb));
SRb_mat = sparse(cell2mat(SRb));
SpRb_mat = sparse(cell2mat(SpRb));
IdRb_mat = speye(m);

% Compute kronecker products of the blocks' matrices
HR_mat = kron(I4_mat, HRb_mat) + x*(kron(SpL2_mat', SpRb_mat) ...
            + kron(SpL2_mat, SpRb_mat') + kron(HXY2_mat, IdRb_mat)) ...
            + (mu/2)*kron(M2_mat, IdRb_mat) ...
            + y*(kron(H2_mat, IdRb_mat) + kron(Q2_mat, SRb_mat) ...
            + (2*n - 2)*kron(Q2_mat_Sq, IdRb_mat) ...
            + (theta^2/(2*pi)^2)*kron(I4_mat, IdRb_mat));            
% build new QR             
QR_mat = kron(I4_mat, QRb_mat) + kron(Q2_mat, IdRb_mat); 
             
% build new SR
SR_mat = kron(I4_mat, SRb_mat) + kron(S2_mat, IdRb_mat) ...
                 + 2*(2*n - 2)*kron(Q2_mat, IdRb_mat);
% build new SpR and IdR
SpR_mat = kron(SpR2_mat, IdRb_mat);                          
                                                            
% Create an array mR_all of all dimensions (including zeros)of blocks in HR
mRb_all = zeros(1, 2*n);
mRb_all(mRb(2, :)) = mRb(1, :);
mR_all = zeros(1, 2*n + 1);
mR_all = conv(mS(1, :), mRb_all);
% Remove zero elements from mR_all and obtain mR
qR_nz_mR = find(mR_all ~= 0); % get charges of non-zero elements in mR_all
mR = [mR_all(qR_nz_mR); qR_nz_mR];
length_mR = size(mR, 2);  % number of elements in mR


% Create an array mSmRb of dimensions of blocks in part_kron(HS, HRb).
mSmRb = kron(mS(1, :), mRb(1, :));

% Blocks reordering in accordance to quantum sectors
% Generate all pairs (qS, qR), where qS = mS(2, :)=1,2,3 and qR = mRb(2, :) 
[QS, QR] = meshgrid(mS(2, :), mRb(2, :));
qSqR = [QS(:), QR(:)];
% Compute the sum (qS + qR) for each pair
sum_qSqR = sum(qSqR, 2);
% Get permutation of blocks by sorting the pairs based on the sum (qS + qR)
[sum_qSqR, perm] = sort(sum_qSqR);
% Get ordered pairs of quantum numbers qS and qR
qSqR = qSqR(perm, :);
% Construct permutation of indices 
idxP = zeros(sum(mSmRb), 1);
cdx = 1;    % current index
for i = 1:3*size(mRb, 2)
    idxP(cdx: cdx + mSmRb(perm(i)) - 1) = ...
        sum(mSmRb(1: perm(i) - 1)) + 1 : sum(mSmRb(1: perm(i))); 
    cdx = cdx + mSmRb(perm(i));
end


%%%%%%%%%% Permute matrices and compose cells for the right blocks %%%%%%%%
KLS = K_mat(mRb(1, :), mS(1, :));
GLS = G_mat(mRb(1, :), D);
PLS = speye(D*m);
PLS = PLS(idxP, :);
% Total permutation matrix
per_mat = PLS*KLS*GLS;
% Permutation of rows and columns
per = per_mat*(1:D*m)';

% Apply permutation to the right blocks matrices
HR_mat = HR_mat(per, per);
QR_mat = QR_mat(per, per);
SR_mat = SR_mat(per, per);
SpR_mat = SpR_mat(per, per);

% Pre-allocate new cells for the final right blocks
HR0 = cell(length_mR, length_mR);
QR0 = cell(length_mR, length_mR);
SR0 = cell(length_mR, length_mR);
SpR0 = cell(length_mR, length_mR);
% Represent matrices as cells for each pair of quantum numbers
HR0 = mat2cell(HR_mat, mR(1, :), mR(1, :));
QR0 = mat2cell(QR_mat, mR(1, :), mR(1, :));
SR0 = mat2cell(SR_mat, mR(1, :), mR(1, :));
SpR0 = mat2cell(SpR_mat, mR(1, :), mR(1, :));

%%%%%%%%%% Keep only needed quantum sectors in the right blocks %%%%%%%%%%%
% Get needed spin domains 
[ sdn, ~ ] = spin_domain(2*n, 2*N, s);
% Compute new list of non-zero dimensions and their indices of blocks HR
% which are in the needed domain of quantum numbers sdn
idx_sdn = ismember(mR(2, :), sdn);
% New list of quantum sectors dimensions and their quantum numbers 
% within needed domain
mR = mR(:, idx_sdn);
length_mR = size(mR, 2);

% Pre-allocate cells for the final right blocks within needed spins domain:
HR = cell(length_mR, length_mR);
QR = cell(length_mR, length_mR);
SR = cell(length_mR, length_mR);
SpR = cell(length_mR, length_mR);
% Construct new right blocks within needed spins domain
HR = HR0(idx_sdn, idx_sdn);
QR = QR0(idx_sdn, idx_sdn);
SR = SR0(idx_sdn, idx_sdn);
SpR = SpR0(idx_sdn, idx_sdn);


% Generate cells of ledgers TR, ordered as charges in mR
% Get indices of quantum numbers qS from qSqR in the array mS
[~, ind_qS] = ismember( qSqR(:, 1)', mS(2, :));
% Get indices of quantum numbers qR from qSqR in the array mRb
[~, ind_qR] = ismember( qSqR(:, 2)', mRb(2, :));
% Get dimensions mRbc and mSc which are parts of combinations for ledger TR
mRbc = mRb(1, ind_qR);
mSc = mS(1, ind_qS);
% Create an array of (qS, qR, mS, mRb)
qS_qR_mS_mRb = [qSqR, mSc', mRbc'];
% Find logical index of (qS, qR) such that qS + qR - 1 is in sdn
[idx_q_in_sdn, ~] = ismember(sum_qSqR, sdn + 1);
% Trim the array qS_qR_mS_mRb excluding rows which have charges 
% sum of which is not in sdn
qS_qR_mS_mRb = qS_qR_mS_mRb(idx_q_in_sdn, :);
% Get number of combinations of qS and qR for each quantum number q in mR
[unique_vals, ~, idx_u_vals] = unique(sum_qSqR(idx_q_in_sdn), 'stable');
counts_of_q = accumarray(idx_u_vals, 1);
% Generate cells of ledgers TR
TR = cell(size(mR, 2), 1);
TR = mat2cell(qS_qR_mS_mRb, counts_of_q, 4);

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


