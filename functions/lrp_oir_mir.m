%% OIR - computation of MIR and of its time-domain decomposition and frequency-domain expansion for two blocks of processes

function out=lrp_oir_mir(Am,Su,q,Mv,i_1,i_2)

% indexes of the two blocks to analyze inside the Q time series
[i1,i2]=oir_subindexes(Mv,i_1,i_2);

% reduced model with the two blocks to analyze - [i1 i2]
ret = lrp_MIR(Am,Su,q,[i1],[i2]);
out.I12=ret.Ixy;

end