%% (Exercise) Build a single level fast multipole scheme for Laplace's
% equation in two dimensions, using only multipole expansions. 
%
% For fixed N, large enough, say 100k, how does the speed of your method
% change as you vary the number of boxes. Is the optimal value close
% to what you would've expected.
%
% Using the optimal number of boxes in each direction, find the cross-over
% point at which the single level fmm starts beating an N^2 version
% of the algorithm. What is the cross-over point?
%
% Repeat the same exercises but for the helmholtz equation in two regimes.
% In one case fix the wave number and increase N.
%
% In the second regime, suppose that k scales like \sqrt{N}. Explain
% your findings.
