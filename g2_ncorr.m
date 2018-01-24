function g2 = g2_ncorr(N_A,N_B)
%Evaluate g2 correlation function given counted detection events
%   
%   N_A: n-dim array; number of observed detection events in A; last dim separates experimental shots; 
%   N_B: n-dim array; number of observed detection events in B; last dim separates experimental shots; 
%
%   g2: n-1 dim array; g2 correlation function
%

adim=ndims(N_A);    % exp-shot dimension; i.e. dim to average in

g2 = mean(N_A.*N_B,adim)./(mean(N_A,adim).*(mean(N_B,adim)));
% omitnan is unnecessary since N is robust against empty-count NaN problems

end