% Test implementation of an equilibrium calculator
%
% Equilibrium conditions are P, T and Ni
% Only single-sublattice phases are considered
% All phases dissolve all elements

%----------------------------------
% Number of elements
nel=2;
% Total number of phases in system
nphtot=2;
% Temperature (equilibrium condition)
T=1000.0;
% Pressure (equilibrium condition)
P=1.0e5;
% Number of moles of each component (equilibrium conditions)
Ni=[0.1 0.9];



% 1. Make an initial guess of the set of variables and the variable values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 2. Convergence loop start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 3. Construct the system of linear equations from the present state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 4. Scale the equations if necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: This is not strictly necessary.


% 5. Solve for updates of the variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note: The default "slash" solver did not converge. Use, for example, qmr


% 6. Do a line search and select an optimum step length for the variable update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: This is not strictly necessary.
% Note: Check that variables do not obtain un-physical values
% Note: Also update the dependent phase composition variables

% 7. Update the variable values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 8. Check if the set of variables should be changed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 9. Check convergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


