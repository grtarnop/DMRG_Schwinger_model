% v1. 12/11/2023
% v2. 05/31/2024  (added multiple levels)
% This function implements DMRG truncation for the Schwinger model left and right 
% blocks

function [HL_new, QL_new, SL_new, SpL_new, IdL_new, UL_new, ...
          HR_new, QR_new, SR_new, SpR_new, IdR_new, UR_new] = ... 
              truncate_Schwinger_blocks_ss(HLss, QLss, SLss, SpLss, IdLss, ...
                                           HRss, QRss, SRss, SpRss, IdRss, ...   
                                           psiLR, N, s, Ntot, m, dynm, trun_tol)
                                         
% spin projections domain for Left Blocks(in the form of integer indices)
if s >= 0
    % spin projections domain for Left Blocks(in the form of integer indices)
    if (N(1) >= 2) && (N(1) <= Ntot - s)
        sdnL = [-N(1)/2: N(1)/2] + N(1)/2 + 1;
    elseif (N(1) > Ntot - s) && (N(1) <= Ntot + s)
        sdnL = [N(1)/2 - Ntot + s: N(1)/2] + N(1)/2 + 1;
    elseif (N(1) > Ntot + s) && (N(1) <= 2*Ntot - 2)
        sdnL = [N(1)/2 - Ntot + s : Ntot + s - N(1)/2] + N(1)/2 + 1;
    end
     % spin projections domain for Right Blocks(in the form of integer indices)
    if (N(2) >= 2) && (N(2) <= Ntot - s)
        sdnR = [-N(2)/2: N(2)/2] + N(2)/2 + 1;
    elseif (N(2) > Ntot - s) && (N(2) <= Ntot + s)
        sdnR = [N(2)/2 - Ntot + s: N(2)/2] + N(2)/2 + 1;
    elseif (N(2) > Ntot + s) && (N(2) <= 2*Ntot - 2)
        sdnR = [N(2)/2 - Ntot + s : Ntot + s - N(2)/2] + N(2)/2 + 1;
    end
elseif s < 0
   % spin projections domain for Left Blocks(in the form of integer indices)
    if (N(1) >= 2) && (N(1) <= Ntot - abs(s))
        sdnL = [-N(1)/2: N(1)/2] + N(1)/2 + 1;
    elseif (N(1) > Ntot - abs(s)) && (N(1) <= Ntot + abs(s))
        sdnL = [-N(1)/2: Ntot - abs(s) - N(1)/2] + N(1)/2 + 1;
    elseif (N(1) > Ntot + abs(s)) && (N(1) <= 2*Ntot - 2)
        sdnL = [N(1)/2 - Ntot - abs(s): Ntot - abs(s) - N(1)/2] + N(1)/2 + 1;
    end    
    % spin projections domain for Right Blocks(in the form of integer indices)
    if (N(2) >= 2) && (N(2) <= Ntot - abs(s))
        sdnR = [-N(2)/2: N(2)/2] + N(2)/2 + 1;
    elseif (N(2) > Ntot - abs(s)) && (N(2) <= Ntot + abs(s))
        sdnR = [-N(2)/2: Ntot - abs(s) - N(2)/2] + N(2)/2 + 1;
    elseif (N(2) > Ntot + abs(s)) && (N(2) <= 2*Ntot - 2)
        sdnR = [N(2)/2 - Ntot - abs(s): Ntot - abs(s) - N(2)/2] + N(2)/2 + 1;
    end    
end
       
lvs = size(psiLR, 1);  % number of energy levels   
                                        
% create new left blocks by using loops over all spin sectors
for q1 = 1:(N(1) + 1)
    for q2 = 1:(N(1) + 1)
        if ismember(q1, sdnL) && ismember(q2, sdnL)
            HL_new{q1, q2} = HLss{q1, q2};
            QL_new{q1, q2} = QLss{q1, q2};
            SL_new{q1, q2} = SLss{q1, q2};
            SpL_new{q1, q2} = SpLss{q1, q2};
            IdL_new{q1, q2} = IdLss{q1, q2};
        end
    end
end            
% create new right blocks by using loops over all spin sectors
for q1 = 1:(N(2) + 1)
    for q2 = 1:(N(2) + 1)       
        if ismember(q1, sdnR) && ismember(q2, sdnR)
            HR_new{q1, q2} = HRss{q1, q2};
            QR_new{q1, q2} = QRss{q1, q2};
            SR_new{q1, q2} = SRss{q1, q2};
            SpR_new{q1, q2} = SpRss{q1, q2};
            IdR_new{q1, q2} = IdRss{q1, q2};
        end
    end
end                                
                                                                            
sct = spin_combinations(N, s);  % combinations of spin sectors
k = size(sct, 1);          % number of combinations of the LR spin sectors
% for each LR combination of the spin sectors find dimension of left and 
% right blocks  
for j = 1:k
    % dimensions of the left block for spin combination j
    dsl(j) = size(HLss{sct(j, 1) + N(1)/2 + 1, sct(j, 1) + N(1)/2 + 1}, 1);  
    % dimensions of the right block for spin combination j
    dsr(j) = size(HRss{sct(j, 2) + N(2)/2 + 1, sct(j, 2) + N(2)/2 + 1}, 1); 
end
ds = sum(dsl.*dsr);          % total dimension of the superblock sector s
dsl_tot = sum(dsl);          % dimension of all the left blocks
dsr_tot = sum(dsr);          % dimension of all the right blocks

% find density matrices for each level and each combination of spins 
rhoL = cell(lvs, k);
rhoR = cell(lvs, k);
for l = 1:lvs
    for j = 1:k
        rhoL{l, j} = psiLR{l, j}*psiLR{l, j}';
        rhoR{l, j} = transpose(psiLR{l, j}'*psiLR{l, j});
    end
end

% find average density matrices for each spin combination
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


% eigenstates and probabilities matrices
UL = cell(k, 1);
pL = cell(k, 1);
UR = cell(k, 1);
pR = cell(k, 1);
for j = 1:k
    [UL{j}, pL{j}] = eig(rhoLav{j});
    [UR{j}, pR{j}] = eig(rhoRav{j});
end

% choose proper bond dimension m 
if dynm == 1 
    for j = 1:dsl_tot
        m = j;
        [posL, tr_errorL] = find_highest_probabilities(pL, m, dsl);
        [posR, tr_errorR] = find_highest_probabilities(pR, m, dsr);
            if tr_errorL < trun_tol &&  tr_errorR < trun_tol
                break
            end
    end
elseif dynm == 0
   if dsl_tot < m
        m = dsl_tot;
   end
   [posL, tr_errorL] = find_highest_probabilities(pL, m, dsl);
   [posR, tr_errorR] = find_highest_probabilities(pR, m, dsr);
end
disp(['Bond dimension = ', num2str(m)])  
disp(['truncation error left = ', num2str(tr_errorL)])
disp(['truncation error right = ', num2str(tr_errorR)])


% truncate matrices UL and UR in each combination of the spin sectors
for j = 1:k
    UL{j} = UL{j}(:, posL{j});
    UR{j} = UR{j}(:, posR{j});
end


% truncate Hamiltonian blocks in the given combinations of spin sectors
% right side truncation
for j2 = 1:k
    for q1 = 1:(N(1) + 1)
        if ismember(q1, sdnL)
            ql2 = sct(j2, 1) + N(1)/2 + 1; 
            HL_new{q1, ql2} = HL_new{q1, ql2}*UL{j2};
            QL_new{q1, ql2} = QL_new{q1, ql2}*UL{j2};
            SL_new{q1, ql2} = SL_new{q1, ql2}*UL{j2};
            SpL_new{q1, ql2} = SpL_new{q1, ql2}*UL{j2};
            IdL_new{q1, ql2} = IdL_new{q1, ql2}*UL{j2};
        end
    end
    for q1 = 1:(N(2) + 1)
        if ismember(q1, sdnR)
            qr2 = sct(j2, 2) + N(2)/2 + 1; 
            HR_new{q1, qr2} = HR_new{q1, qr2}*UR{j2}; 
            QR_new{q1, qr2} = QR_new{q1, qr2}*UR{j2};
            SR_new{q1, qr2} = SR_new{q1, qr2}*UR{j2};
            SpR_new{q1, qr2} = SpR_new{q1, qr2}*UR{j2};
            IdR_new{q1, qr2} = IdR_new{q1, qr2}*UR{j2};        
        end
    end
end
% left side truncation
for j1 = 1:k
    for q2 = 1:(N(1) + 1)
        if ismember(q2, sdnL)
            ql1 = sct(j1, 1) + N(1)/2 + 1; 
            HL_new{ql1, q2} = UL{j1}'*HL_new{ql1, q2};
            QL_new{ql1, q2} = UL{j1}'*QL_new{ql1, q2};
            SL_new{ql1, q2} = UL{j1}'*SL_new{ql1, q2};
            SpL_new{ql1, q2} = UL{j1}'*SpL_new{ql1, q2};
            IdL_new{ql1, q2} = UL{j1}'*IdL_new{ql1, q2};
        end
    end
    for q2 = 1:(N(2) + 1)
        if ismember(q2, sdnR)
            qr1 = sct(j1, 2) + N(2)/2 + 1; 
            HR_new{qr1, q2} = UR{j1}'*HR_new{qr1, q2}; 
            QR_new{qr1, q2} = UR{j1}'*QR_new{qr1, q2};
            SR_new{qr1, q2} = UR{j1}'*SR_new{qr1, q2};
            SpR_new{qr1, q2} = UR{j1}'*SpR_new{qr1, q2};
            IdR_new{qr1, q2} = UR{j1}'*IdR_new{qr1, q2};      
        end
    end
end

UL_new = UL;
UR_new = UR;

end


% This function finds m highest probabilities and positions in their spin
% sector combinations
% Input: cells (for each sector) of probabilities: p
%        number of states kept: m
%        list of dimensions of each combination of spins: ds
% Output: cells (for each sector) of positions of the highest probabilites
function [pos, truncation_error] = find_highest_probabilities(p, m, ds)


k = size(p, 1);  % number of combinations of the spin sectors
p_list = [];     % initialize a list of all probabilities 
% add probabilities from each combination of spin sectors to a single list 
for j = 1:k
    p_list = vertcat(p_list, diag(p{j}));  
end

% if the size of the total proablility list turned out to be smaller than m
% replace m by the total size of the list
if size(p_list, 1) < m
    m = size(p_list, 1);
end

% sort the list of all probabilities and get the sorting index 
[p_sort, ind] = sort(p_list, 'descend');
% compute truncation error 
truncation_error = 1 - sum(p_sort(1:m));
% construct a list of combinations and positions of the highest probabilities
sec_list = [];
for i = 1:m
    % find a combination and position within it 
    % for a given absolute position in the full list
    [js, ls] = find_sector_and_position(ind(i), ds);
    sec_list = vertcat(sec_list, [js, ls]);
end
% sort the list of sectors and positions 
sec_list = sortrows(sec_list);

% initialize k cells of positions of the highest probabilities
pos = cell(k, 1);
% add positions of the highest eigenvalues for each combination cell
for i = 1:m
    pos{sec_list(i,1)} = horzcat(pos{sec_list(i,1)}, sec_list(i,2));
end


end


% This function finds spin combination j (j = 1,...., k) and position l_j within
% this spin combination of an elemenent whose absolute index is i
% Input: i and list of sector dimensions ds = [ds(1), ds(2),..., ds(k)] 
% Output: spin combination j and position within this combination l
function [j, l] = find_sector_and_position(i, ds)

l = i;  % set initial position to i
j = 1;  % set initial spin combination to 1
while l > 0  % while position bigger than zero
    l = l - ds(j);  % subtract dimension ds(j) of the combination j
    j = j + 1;      % move to the next spin combination
end

j = j - 1;          % move one spin combination back
l = l + ds(j);      % and add previous dimension to the position
  
end


