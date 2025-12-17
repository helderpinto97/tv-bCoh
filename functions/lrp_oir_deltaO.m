%% OIR - computation of delta of the O-information rate when the block Xj is added to the group X-j

% Am: VAR coefficients - to be [A(1),..,A(p)]
% Su: Covariance matrix of residuals Q x Q
% q:  lag to truncate autocorrelation function
% Mv: vector specifying how many processes for each block
% ij: complete vector of indexes
% j:  index of the target to analyze within ij

%%% OUTPUT
% dO: delta OIR
% dOf: spectral delta OIR

function out = lrp_oir_deltaO(Am,Su,q,Mv,ij,j)

%%% computation of dO(X-j;Xj)
assert(ismember(j,ij)); % verify target belongs to group
ii=setdiff(ij,j); % driver indexes
N=length(ij); % order of the OIR to compute

out1=lrp_oir_mir(Am,Su,q,Mv,ii,j);

i_cs=nchoosek(ii,N-2); % number of combinations to compute the sum of the deltaO 
dO=0;
for cnt=1:N-1
    outtmp=lrp_oir_mir(Am,Su,q,Mv,i_cs(cnt,:),j);% MIR between Xj and X-ij
    % time domain measure
    dO=dO+outtmp.I12;
end
% time domain measures: last term 
dO=dO+(2-N)*out1.I12; % dO(X-j;Xj)

%% output
out.dO=dO;

end

