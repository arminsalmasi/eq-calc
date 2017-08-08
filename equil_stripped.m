clear variables; close all; clc;
% Test implementation of an equilibrium calculator
%
% Equilibrium conditions are P, T and Ni
% Only single-sublattice phases are considered
% All phases dissolve all elements

%----------------------------------
% Number of elements
nel = 2;
% Total number of phases in system
nphtot = 2;
% Temperature (equilibrium condition)
T = 1000.0;
% Pressure (equilibrium condition)
P = 1.0e5;
% Number of moles of each component (equilibrium conditions)
Ni = [0.4,0.6];

% 1. Make an initial guess of the set of variables and the variable values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nph is number of phases currently assumed to be stable
nph = 2;
% status specifies whether a phase is assumed stable (=1), or not (=0)
status = [1,1];
% xp contain the phase compositions
xp = [0.5, 0.5 ; 0.3, 0.7]   % Initial guess , 2 phases are pure
% np is the number of moles of each phase
np = [0.4, 0.6]  % Initial guess , 2 phases are pure

%%
% For phase 'index' evaluate
% Gibbs energy 'G'
% first derivatives 'dG' in order dG/dx dG/dn
% second derivatives 'd2G' in order
% d2G/dx2 d2G/dxdn d2G/dn2
% Chemical potentials mu
% The independent composition variable is given by scalar 'x'
% The number of moles form ula units are given by scalar 'n'
%  [G,dG,d2G,mu]=binary_regular_solution(index,x,n,T)
for i = 1 : nph 
    [G(i),dG(i,:),d2G(i,:),mu(i,:)] = binary_regular_solution(i,xp(i,1),np(1,i),T); %phase 2, x12,n2
end



%%
% the constraints g(nel,1)
% the 1st derivatives of the constraints dg(nel,nel*nph)
% the 2nd derivatives of the constraints d2g(nel,(n^2+n)/2) n=nel*nph
%     d2g is for each row given in order 
%     d2g/dy{1}dy{1} d2g/dy{2}dy{1} d2g/dy{2}dy{2} d2g/dy{3}dy{1} etc
[g,dg,d2g] = constraints(nel,nph,nphtot,Ni,status,xp,np);
%%%%%%% assume lamb = [1 ;1]
%%%%%L = Gtot + sum(g(:))
% 2. Convergence loop start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conv = false;
delta = zeros(6,1);
c(1,1) = dG(1,1) + mu(1,1) * dg(1,1) + mu(1,2)*dg(2,1);               %2
c(2,1) = G(1,1)/(-dg(1,1)) + mu(1,1) * dg(1,2) + mu(1,2)*dg(2,2);  %1
c(3,1) = dG(2,1) + mu(2,1) * dg(1,3) + mu(2,2)*dg(2,3); %mu?          %4
c(4,1) = G(1,2)/(-dg(1,3)) + mu(1,1) * dg(1,4)+mu(1,2)*dg(2,4); %mu? %3
c(5,1)= g(1,1); 
c(6,1)= g(2,1); 
c = -1*c; 
cnt=1;
while (cnt<100)
  % 3. Construct the system of linear equations from the present state
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  A= [zeros(4,4),(dg)';dg,zeros(2,2)];
  A(1,1) = d2G(1,1);
  A(1,2) = dG(1,1)-mu(1,1)+mu(1,2); %mu?
  A(2,1) = A(1,2);
  A(3,3) = d2G(2,1);
  A(3,4) = dG(2,1)-mu(1,1)+mu(1,2);% mu?
  A(4,3) = A(3,4);
  % 4. Scale the equations if necessary
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Note: This is not strictly necessary.


  % 5. Solve for updates of the variables
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Note: The default "slash" solver did not converge. Use, for example, qmr
  delta=A/c';
  
     
  % 6. Do a line search and select an optimum step length for the variable update
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Note: This is not strictly necessary.
  % Note: Check that variables do not obtain un-physical values
  % Note: Also update the dependent phase composition variables
  %[x11,n1,x12,n2,lam1,lam2]
  
  
%% loop to minimize r and find teta
    xptmp(1,:) = xp(1,:) ;    %2 
    nptmp(1) = np(1);   %1
    xptmp(2,:) = xp(2,:) ;   %4
    nptmp(2) = np(2)  ;     %3
     k = 1 ;
     teta = 0.000001;
     curl = 1;
     while curl > 0.00001
        pk(:,1) = sign(delta(:,1));

        xptmp(1,1) = xptmp(1,1) + k * teta * pk(1,1);    %2 
        nptmp(1) = nptmp(1) + k * teta * pk(2,1)   ;   %1
        xptmp(2,1) = xptmp(2,1) + k * teta * pk(3,1) ;   %4
        nptmp(2) = nptmp(2) + k * teta * pk(4,1) ;     %3
        r(k)= rp(xp,np,nel,nphtot,Ni,status,nph,T);
        k=k+1;

         for i = 1 : nph 
            [G(i),dG(i,:),d2G(i,:),mu(i,:)] = binary_regular_solution(i,xptmp(i,1),nptmp(1,i),T); %phase 2, x12,n2
         end
         [g,dg,d2g] = constraints(nel,nph,nphtot,Ni,status,xptmp,nptmp);
         curl = sum(d2G(1,:)) + sum(d2G(2,:)) + sum(d2g(1,:)) +sum(d2g(2,:));
      end
    teta = k *teta;  
    %plot(r(:))

%% end of the loop to minimize r and find teta

  % 7. Update the variable values
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xp(1,1) = xp(1,1) + teta * delta(1,1);    %2 
    np(1) = np(1) +  teta * delta(2,1)   ;   %1
    xp(2,1) = xp(2,1) + teta* delta (3,1) ;   %4
    np(2) = np(2) + teta * delta(4,1) ;     %3
    xp(1,2) =  1- xp(1,1);
    xp(2 ,2) = 1- xp(2,1);  
  
  % 8. Check if the set of variables should be changed
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % 9. Check convergence
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cnt= cnt+1;
end
xp=xp
np=np



function [r] = rp(xp,np,nel,nphtot,Ni,status,nph,T)
      for i = 1 : nph 
        [G(i),dG(i,:),d2G(i,:),mu(i,:)] = binary_regular_solution(i,xp(i,1),np(1,i),T); %phase 2, x12,n2
      end
      [g,dg,d2g] = constraints(nel,nph,nphtot,Ni,status,xp,np);
      %mu(5,1)= mu(1,1) + teta * c(5,1); 
      %mu(6,1)= mu(1,2) + teta * c(6,1); 
      %r = sum( temp(1:4,1) .*temp(1:4,1) ) + sum(temp(5:6,1).*c(5:6,1));
      r = ...
        ( G(1,1) + mu(1,1) * dg(1,1) + mu(1,2)*dg(2,1) )^2 + ...
        ( G(1,1)/(-dg(1,1)) + mu(1,1) * dg(1,2) + mu(1,2)*dg(2,2) )^2 + ...
        ( dG(2,1) + mu(2,1) * dg(1,3) + mu(2,2)*dg(2,3) )^2 + ...
        ( G(1,2)/(-dg(1,3)) + mu(1,1) * dg(1,4)+mu(1,2)*dg(2,4) )^2 + ...
        ( g(1,1) )^2 + ...
        ( g(2,1) )^2 ;
end
%%
function alphak = linesearch(f,d,x,rho,c)
% function alphak = linesearch(f,d,x,rho,c)
% Backtracking line search
% See Algorithm 3.1 on page 37 of Nocedal and Wright
% Input Parameters :
% f: MATLAB file that returns function value
% d: The search direction
% x: previous iterate
% rho :- The backtrack step between (0,1) usually 1/2
% c: parameter between 0 and 1 , usually 10^{-4}
% Output :
% alphak: step length calculated by algorithm
% Kartik's MATLAB code (27/1/08)

alphak = 1;
[fk, gk] = feval(f,x);
xx = x;
x = x + alphak*d;
fk1 = feval(f,x);
while fk1 > fk + c*alphak*(gk'*d)
  alphak = alphak*rho;
  x = xx + alphak*d;
  fk1 = feval(f,x);
end
end