%% Computation of the gradient of OIR and OIR for blocks of processes 

%%% INPUT
% Am_n: Time-varying VAR coefficients matrix Q x Qp xN
% Su_n: Time-varying residuals covariance matrix Q x Q x N
% Mv: vector specifying how many processes for each block
% q:  lag to truncate autocorrelation function

%%% OUTPUT
% dO: time domain delta OIR at each time step
% OIR: time domain OIR at each time step

function out = tv_lrp_oir(Am_n,Su_n,q,Mv)

% check inputs
if nargin < 3
    error('Not enough input arguments')
end

M = length(Mv); % number of the blocks
p=size(Am_n,2)/size(Am_n,1); % model order
allM=cell(M,1); % number of combination for all orders
dO=cell(M,1); % deltaOIR

for N=3:M
    comb=nchoosek(1:M,N); % all multiplets of size N
    allM{N}=comb;
    for n=1:size(comb,1) % cycle on all the combinations of that order
        ij=comb(n,:); % select the combination
        for k=1:size(comb,2) % vary target inside the multiplet
            for t=p+1:size(Su_n,3)
                out1=lrp_oir_deltaO(Am_n(:,:,t),Su_n(:,:,t),q,Mv,ij,ij(k));
                dO12(n,k,t)=out1.dO;
            end
        end
    end
    dO{N-1}=dO12;
end

OIR=cell(M,1); % OIR
OIRf=cell(M,1); % spectral OIR
OIR{3}=squeeze(dO{2}(:,1,:)); % other columns are the same

for N=4:M
    t1=1; % index of Xj (can be any number btw 1 and N)
    t2=setdiff(1:N,t1); % index of X-j
    for cnt=1:size(allM{N},1) % cycle on all combinations of order N
        tmp=allM{N}(cnt,:);
        iN1=find(sum(tmp(t2)==allM{N-1},2)==N-1); % position of X-j in the O-info one step back
        OIR{N}(cnt,:) = squeeze(OIR{N-1}(iN1,:)) + squeeze(dO{N-1}(cnt,t1,:))';
        
    end
end

%% OUTPUT
out.dO=dO;
out.OIR=OIR; 
out.allM=allM;

end

