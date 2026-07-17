%D6M_EXERCISE_OPEN_PROBLEMS Day 6 morning open problems.
%
% Choose ONE of the projects below.  These are research problems, not routine
% programming exercises.  A useful submission may contain a construction, a
% counterexample, numerical evidence, an implementation attempt, or a clear
% explanation of what prevented further progress.
%
% You may use AI tools, but you are responsible for checking every statement,
% experiment, and citation in your submission.

%% Project 1. Low-rank completions of positive sparse matrices
%
% Lawrence Saul proposed the following conjecture.  Start with a matrix A
% whose nonzero entries are positive.  Keep these entries fixed.  Replace each
% zero entry by a negative real number.  Can the completed matrix always be
% made low rank?
%
% Write the question precisely for an m-by-n matrix.  Let S be the locations
% of the positive entries.  Seek a matrix X such that
%
%   X(i,j) = A(i,j)  for (i,j) in S,
%   X(i,j) < 0       for (i,j) not in S.
%
% Investigate how small rank(X) can be.  You must say what "low rank" means in
% your experiments.  For example, try to find completions of rank 1, 2, or 3,
% or compare rank(X) with min(m,n).
%
% Possible directions:
%
%   * Work out small examples by hand, such as diagonal, banded, or block
%     patterns.
%   * Look for sparsity patterns that rule out rank 1 or rank 2.
%   * Write a numerical search for low-rank completions with the required
%     signs, then test families of sparse matrices.
%   * State a narrower conjecture suggested by your evidence, or give a
%     counterexample to a proposed version of the conjecture.
%
% Numerical rank depends on a tolerance.  State that tolerance and do not call
% a matrix exactly low rank merely because one singular value is small.

%% Project 2. RCIP for three-dimensional edges and vertices
%
% Recursive compressed inverse preconditioning (RCIP) resolves singular
% behavior near corners in two-dimensional boundary integral equations.  Think
% about extending this idea to a three-dimensional boundary with an edge or a
% vertex.
%
% Choose a simple model problem, for example the Laplace equation on a
% polyhedral domain, and describe a possible method.  Address the following
% questions in plain language:
%
%   * What local patch near an edge or vertex would be refined?
%   * What boundary integral equation and unknown density would be used?
%   * What local operator would be compressed or inverted recursively?
%   * How would you check whether the method improves accuracy or conditioning?
%
% An implementation is welcome but is not required.  A careful plan based on a
% small model geometry is more valuable than a large unsupported claim.

%% Project 3. MFS with learned source placement near 3D singularities
%
% The method of fundamental solutions (MFS) approximates a solution using
% fundamental solutions centered at sources outside the physical domain.  A
% possible research direction is to use machine learning to choose source
% locations or other MFS parameters for problems with edges or vertices in
% three dimensions, in the spirit of lightning-type methods.
%
% Pick a simple family of geometries and boundary data.  Compare any learned
% strategy with an ordinary MFS baseline using fixed source locations.  Report
% both the error and the conditioning of the linear system.  Keep separate
% test geometries or boundary data that were not used to choose the learned
% parameters.
%
% A useful result may be that the learned strategy does not help, provided the
% experiment is reproducible and the comparison is fair.

%% Submission
%
% Submit a short research note with:
%
%   1. A precise problem statement and the project you chose.
%   2. Your method, calculations, or numerical experiment.
%   3. Evidence: examples, plots, tables, or a counterexample.
%   4. Limitations, negative results, and the next question you would study.
%
% Include runnable code when you write code, and give complete references for
% any material outside the course notes.
