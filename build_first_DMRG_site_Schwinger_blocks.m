% Dempsey et all convention
% v1. 09/28/2024 
% This function initializes first DMRG site blocks for the Schwinger model
% (Dempsey convention)
function [HLb, QLb, SLb, SpLb, mLb, HRb, QRb, SRb, SpRb, mRb] = ...
                    build_first_DMRG_site_Schwinger_blocks(x, mu, y, theta)


% Single spin-1/2 site operators 
I2 = sparse([1 0; 0 1]);       % Identity matrix
sp = sparse([0 1; 0 0]);       % sigma plus matrix
sz = sparse([1 0; 0 -1]);      % Sz Pauli matrix

% Basic double-site operators 
% Hamiltonian operator HC for 2 sites  (Dempsey convention)       
H2 = (1/2)*kron((1 + theta/pi)*(I2 + sz) + theta^2/(2*pi^2)*I2, I2);      
Q2 = (1/2)*(kron(sz, I2) + kron(I2, sz));   % charge operator Q
S2 = kron(sz, I2);                          % operator S             
M2 = kron(I2 + sz, I2) + kron(I2, I2 - sz); % mass operator for 2 sites
HXY2 = kron(sp, sp') + kron(sp', sp);       % XY Hamiltonian for 2 sites
SpL2 = kron(I2, sp);
SpR2 = kron(sp, I2);


% Values for the initial double-site operators
% Left blocks (Dempsey et all convention)
HL1 = x*(kron(sp, sp') + kron(sp', sp)) ...
        + (mu/2)*(kron(I2 + sz, I2) + kron(I2, I2 - sz)) + y*H2;
QL1 = Q2;
SL1 = S2;
SpL1 = kron(I2, sp);
% Right blocks (Dempsey et all convention)
HR1 = x*(kron(sp, sp') + kron(sp', sp)) ... 
        + (mu/2)*(kron(I2 + sz, I2) + kron(I2, I2 - sz)) + y*H2;
QR1 = Q2;
SR1 = S2;
SpR1 = kron(sp, I2);


%%%% Construct all Hamiltonians in spin-spin sectors for 1 DMRG site %%%%%%
% Set of vectors with a given Sz value  
v = cell(3, 1);
for q = 1:3
    v{q} = subspace_fixed_Sz(2, -1 + q - 1);
end
% Left blocks 
HLb = cell(3, 3); QLb = cell(3, 3); SLb = cell(3, 3); SpLb = cell(3, 3);
% Right blocks
HRb = cell(3, 3); QRb = cell(3, 3); SRb = cell(3, 3); SpRb = cell(3, 3);  

for q1 = 1:3      % runs over all 3 different spin projections 
    for q2 = 1:3   % runs over all 3 different spin projections
        % Left blocks
        HLb{q1, q2} = v{q1}'*HL1*v{q2}; 
        QLb{q1, q2} = v{q1}'*QL1*v{q2};
        SLb{q1, q2} = v{q1}'*SL1*v{q2}; 
        SpLb{q1, q2} = v{q1}'*SpL1*v{q2};
        % Right blocks
        HRb{q1, q2} = v{q1}'*HR1*v{q2}; 
        QRb{q1, q2} = v{q1}'*QR1*v{q2};
        SRb{q1, q2} = v{q1}'*SR1*v{q2}; 
        SpRb{q1, q2} = v{q1}'*SpR1*v{q2};
    end
end


% For a single DMRG site create dimensions of blocks 
for q = 1:3
    mLb(1, q) = size(HLb{q, 1}, 1);
    mRb(1, q) = size(HRb{q, 1}, 1);
end
mLb = [mLb; 1:3 ];
mRb = [mRb; 1:3 ];


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates a set of states with the spin projection s
% and chain size N
function [ vec ] = subspace_fixed_Sz(N, s)

% compute all possible configurations of 0 and 1
% with  q  ones and N - q zeros and make a set of vectors vec
q = N/2 + s; 
vec = cell(nchoosek(N,q),1);           % initialize a cell for vectors
count = 1;
   for k = 0:2^N-1
       state = (dec2bin(k,N)-'0');  % represent number k as a binary of length N
       if sum(state) == q           % check if state has charge q (i.e. q ones )
            vec{count} = 1;
                for i = 1:N
                    if state(i) == 1 
                        vec{count} = sparse(kron(vec{count},[1;0]));
                    else
                        vec{count} = sparse(kron(vec{count},[0;1]));  
                    end
                end    
       count = count + 1;
       end
   end
% combine all the vectors into a big matrix 2^N x C(N, N/2 + s)
% so each column has length 2^N and there are  C(N, N/2 + s) columns
vec = sparse([vec{:}]);

end



