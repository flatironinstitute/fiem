%WORKSHEET4_different_bcs
%
% Requirements:
%   - This worksheet will use the chunkIE package to demonstrate the ideas.
%   Please download and install the package according to the directions
%   here: https://chunkie.readthedocs.io/en/latest/getchunkie.html
%
% Scope:
%   - This worksheet demonstrates the need for the right preconditioned
%     combined field representation for the exterior Neumann problem
%   - This worksheet has an exercise for solving boundary value problems
%     with different boundary conditions
%

%% (Exercise 1) Exterior Neumann problem using the Helmholtz single layer
%  potential
%
%  Set up the Helmholtz exterior Neumann boundary value problem using the
%  single layer potential on the unit disk
%
%  Let k be the wavenumber. Let A_{k} be the discretized matrix that you obtain
%  using chunkie. For k\in [3,4] compute the inverse inner product 
%  g(k) = 1/(u, A_{k}^{-1} v) where u and v are smooth functions. 
%
%  Does the function g(k) have any zeros? You can use chebfun to find 
%  the zeros of g(k). If so, find them.
%  If you found a zero, try solving the integral equation at that wavenumber.
%
%  Do you get the desired accuracy? Just to make sure there isn't a bug in your code
%  try running the same test when g(k) \neq 0. 
%
%  Can you explain what goes wrong here? 
%  
%


%% (Exercise 2) The right preconditioned integral representation  
%
% Now try solving the exterior Neumann boundary value problem, using the
% following integral equation
%
% u = \beta(S_{k} + i \alpha D_{k} S_{ik}) \sigma
%
% with \beta = -1.0/(0.5 + 1i*0.25*alpha)
% 
% On imposing the boundary conditions, we get the following integral 
% equation
%
%  du/dn = I + \beta S_{k}'[\sigma] + 
%          i\beta \alpha(D'_{k} - D'_{ik})(S_{ik}[\sigma]) + 
%          i\beta \alpha(S_{ik}')^2 [\sigma];
% 
% Set it up as a system of integral equations for \sigma, \tau, and \mu where
% -S_{ik}[\sigma] + \tau = 0, and - S_{ik}'[\sigma] + \mu=0
%
% Try solving the exterior Neumann problem with this representation when g(k) = 0?
%
% What do you observe?
%


%% (Exercise 3) Multiply connected problems
% 1. Set up a multiply connected geometry using an array of chunkie objects
% 2. Solve the exterior problem specifying Dirichlet boundary conditions on
%    some of the components and Neumann boundary conditions on the rest. 
% 3. Verify your solution using an analytic solution test. 
