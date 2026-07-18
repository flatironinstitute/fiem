% Building generalized gaussian quadratures using the old methods
% Consider the family of functions J_{0}(x) e^{-xt} for t \in [1,4],
% x \in [0,36]. $t$ is the parameter, and x should be treated as the
% domain of the functions. Note that for stability reasons, you should
% not be using a tolerance less than 1e-7 for these tests. These codes
% are often run in extended precision to obtain good double precision
% quadrature nodes and weights.
%
% (Exercise 1) Using the procedure in d2a_trunc_lap.m, determine the numerical
% rank of these functions.
%
% (Exercise 2) Suppose that the numerical rank is m. Let \psi_{0}(x), .. \psi_{m-1}(x)
% denote the orthogonal functions which span your family of functions. 
%
% Let P_{0}(x), ... P_{m-1}(x) denote scaled versions of the Legendre polynomials
% on the interval [0,36]
%
% Let \phi_{j} = (1-\alpha) P_{j} + \alpha \psi_{j}. Starting at \alpha =0, and using
% x_{j}, w_{j} to be the Gauss-Legendre nodes and weights, slowly increase \alpha,
% and solve the non-linear optimization problem to determine a generalized gaussian
% rule for \phi_{j}. For each new \alpha, use the solution from the previous \alpha
% as an intial guess for the Newton problem. 
%
%
% (Exercise 3) Implement the same procedure for the family of functions $x^{\mu}$ 
% for \mu \in [0.5,40]. What do you observe?
%
% (Exercise 4) If you are brave, and have the time, implement the node-eliminiation
% procedure to obtain a generalized Gaussian rule that we discussed at the end
% of the class. 
