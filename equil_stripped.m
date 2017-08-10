clear variables; close all; clc;
% Test implementation of an equilibrium calculator
% Equilibrium conditions are P, T and Ni
% Only single-sublattice phases are considered 

%% Armin: y==x 
%%

% All phases dissolve all elements
%----------------------------------

%% initialization of variables
% Number of elements
nel = 2;
% Total number of phases in system
nphtot_hld(1) = 2;
% Temperature (equilibrium condition)
T = 1000.0;
% Pressure (equilibrium condition)
P = 1.0e5;
% Number of moles of each component (equilibrium conditions)
Ni = [0.5,0.5];
%% Armin: this condition should always hold
Ntot =sum(Ni); 
%%

% 1. Make an initial guess of the set of variables and the variable values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Armin: initialize holder variables (to keep values of each itteration)
% nph is number of phases currently assumed to be stable
nph_hld(1) = nphtot_hld(1); % Armin:statrting guess-phases 1&2 are statble
% status specifies whether a phase is assumed stable (=1), or not (=0)
%% Armin: Todo : 1 phase regions????
status_hld(1,:) = [1,1]; %statrting guess - both phases are statble 
xp_hld(1,:,:) = [0.49 0.51; 0.51 0.49]; % Armin: Initial-2 statble phases
% np is the number of moles of each phase
np_hld(1,:) = [0.49*Ntot 0.51*Ntot] ; % Armin: Initial-2 statble phases
% lmb = lagrangian multipliers, 1*2 vector(2 constraints)
lmb_hld(1,:) = [1e-3 1e-3]; %% % Armin: Initial-2 statble phases
                           
% 2. Convergence loop start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ct1 = 1; %iteration counter
while 1==1
  %% Armin: initialize step variables (form the holder variables)
  xp(:,:) = xp_hld(ct1,:,:);
  np(:) = np_hld(ct1,:);
  status  = status_hld(ct1,:);
  nph = nph_hld(ct1);
  nphtot = nph_hld(ct1);
  lmb(:) = lmb_hld(ct1,:);

  %% Armin: Calculate Gibbs free energy and derivitives (both phases)
  %% Armin: if driving force for one of the phases is zero then??
  %% Armin: d2G should be [nel*nph nel*nph] matrix) : DDxijxij, DDxijnj
  for idx = 1 : nph 
    [G(idx),dG(idx,:),d2G(idx,:),mu(idx,:)] = ...
           binary_regular_solution(idx, xp(idx, 1), np(idx), T);
  end
  %% Armin: constraints and derivitives regards to variables
  %% Armin: variables order[x11,n1,x12,n2,lm1,lm2]:different in 2015 paper
  [g,dg,d2g] = constraints(nel,nph,nphtot,Ni,status,xp,np);
  
  
  % 3. Construct the system of linear equations from the present state
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Armin: transform d2G to 4*4 system with:
  %% Armin: d2G(y1,y2)=0,d2G(y1,n2)=0,d2G(n1,y2)=0,d2G(n1,n2) =0 
  %% Armin: Todo : this should be automatic - look at the G function 
  d2G_A(:) = ...
           [d2G(1,1) d2G(1,2) d2G(1,3) 0 0 d2G(2,1) 0 0 d2G(2,2) d2G(2,3)];
  %% Armin: dimension 6*6
  nd = nel*nph;  
  A(:,:) = [zeros(nd,nd),(dg(:,:))';dg(:,:),zeros(nel,nel)];
  %% Armin: 4*4 G''+g'' part of the matrix
  ct2 = 1; %% Armin: counter variable
  for nd=1:nd %%probe dimension == 4 line of equations
    for ct3=1:nd %% Armin: lower triangle => 4*4
      A(nd,ct3) = d2g(1, ct2)*lmb(1) + d2g(2, ct2)*lmb(2) + d2G_A(ct2);
      A(ct3,nd) = A(nd,ct3);
      ct2 = ct2 + 1;
    end
  end
  %% Armin: L' values (6*1)
  c = -[(dG(1,:) + lmb(1)*dg(1,1:2) + lmb(2)*dg(2,1:2) )';...
        (dG(2,:) + lmb(1)*dg(1,3:4) + lmb(2)*dg(2,3:4) )';...
        g];
  %% Armin: Deviaitions (6*1)
  dlt = zeros(nel*nph+nel,1); 
  
  % 4. Scale the equations if necessary
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Note: This is not strictly necessary.
  %% Armin: what is this part?!
  
  % 5. Solve for updates of the variables
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Note: The default "slash" solver did not converge. Use, for example, qmr
  %% Armin: results with / solver or vpasolve are completely different!?
  %% Armin: qmr command: 
    % X = qmr(A,B,TOL,MAXIT) specifies the maximum number of iterations
         %. If MAXIT is [] then qmr uses the default, min(N,20).
    % [X,FLAG] = qmr(A,B,...) also returns a convergence FLAG:
    % 0 qmr converged to the desired tolerance TOL within MAXIT iterations.
    % 1 qmr iterated MAXIT times but did not converge.
    % 2 preconditioner M was ill-conditioned.
    % 3 qmr stagnated (two consecutive iterates were the same).
    % 4 one of the scalar quantities calculated during qmr became too
    % small or too large to continue computing.
    % [X,FLAG,RELRES] = qmr(A,B,...) also returns the relative residual
    % NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL.
    % [X,FLAG,RELRES,ITER] = qmr(A,B,...) also returns the iteration number
    % at which X was computed: 0 <= ITER <= MAXIT.
    %[X,FLAG,RELRES,ITER,RESVEC] = qmr(A,B,...) also returns a vector of 
        %the residual norms at each iteration, including NORM(B-A*X0).
  [dlt,FLAG,RELRES,ITER,RESVEC] = qmr(A,c,[],1e6);
  %% Armin: break if solution contains complex answer!
  if not(isreal(dlt))
    err = 'complex resolts for delta'
    break % out of the convergance loop(while)
  end
  %% Armin: to update everything in one run!
  vrs = [xp(1,1); np(1); xp(2,1); np(2); lmb(1); lmb(2)]; 

  % 6. Do a line search and select an optimum step length for the variable
  %  update
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Note: This is not strictly necessary.
  % Note: Check that variables do not obtain un-physical values
  % Note: Also update the dependent phase composition variables
  %% Armin: this is not a line search Only search and elimination!:
  %%  calculate r for different values of theta and remove the thetas which
  %%  which violate constraints and selects the one whihc minimize r 
  %% Armin: ToDo: what if the desired theta is outside of this range?
  tta=[-logspace(1,-10,10000), logspace(-10,1,10000)];
  r=sum((A*(vrs+dlt*tta)).^2);
  ct3=1; %% Armin: itteration counter on lenght of theta
  %% remove the thetas which violate constraints
  while ct3<=length(tta) 
    for nd =1 : 2 : nel*nph
      %% Armin:check condition on 0<x<1  
      cnd =...
          (vrs(nd)+tta(ct3)*dlt(nd)>1 || vrs(nd)+tta(ct3)*dlt(nd)<0);
      if cnd
        tta(ct3)=[];
        r(ct3)=[];
        ct3=ct3-1;
       break 
       %% to the for loop
      end
      %% Armin: check condition 0<np<Ntot
      cnd = ...
     (vrs(nd+1)+tta(ct3)*dlt(nd+1)>Ntot || vrs(nd+1)+tta(ct3)*dlt(nd+1)<0);
      if cnd
        tta(ct3)=[];
        r(ct3)=[];
        ct3=ct3-1;
        break 
        %% to the for loop
      end
    end
    ct3=ct3+1; 
  end
  %% Armin: find the theta which minimize r (in the remaining thetas)
  tta=tta(r==min(r)); 
  %% Armin: error handler
  if (isempty(tta)) 
    err = 'all thetas in the range violated xp and np constraints'  
    return 
    %% to the main
  end
  %% Armin: in case there are more than one mins (t1=t2=...)
  tta = tta(1);
    
  % 7. Update the variable values
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  vrs=vrs+tta*dlt;
  xp=[vrs(1) 1-vrs(1); vrs(3) 1-vrs(3)];
  np=[vrs(2) vrs(4)];
  lmb=[vrs(5); vrs(6)];

  %% Armin: save values of the current iteration (for tracking)
  for nd = 1: 2
    G_hld(ct1, nd)= G(nd);
    dG_hld(ct1,nd,:)= dG(nd,:);
    d2G_hld(ct1,nd,:)= d2G(nd,:);
    mu_hld(ct1, nd,:)= mu(nd,:);
  end
  g_hld(:,ct1) = g(:);
  dg_hld(:,:,ct1) = dg(:,:);
  d2g_hld(:,:,ct1) = d2g(:,:);
  d2GtoA_hld(ct1,:) = d2G_A(:);   
  delta_hld(:,ct1)= dlt;
  r_hld(ct1) = min(r);
  c_hld(:,ct1) = c;
  %% Armin: initialize variables for next itteration step
  lmb_hld(ct1+1,:) =lmb(:) ;
  status_hld(ct1+1,:) = status(:);
  nph_hld(ct1+1) = nph ;
  nphtot_hld(ct1+1) = nphtot  ;
  xp_hld(ct1+1,:,:) = xp(:,:) ;
  np_hld(ct1+1,:)= np(:);
  % 8. Check if the set of variables should be changed
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Armin: ?? is it checking the staus of the phases?!
  
  % 9. Check convergence
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Armin: ToDo: What is a good convergence criteria? 
  if ct1>1
    if (min(r) - r_hld(ct1-1))<0 || ct1>1000 %% 
      del_r = (min(r) - r_hld(ct1-1));
      break  
      %% to the main / end of itteration conditions meet
    end
  end
  %% Armin: update counter variable
  ct1 =ct1+1;            
  %% Armin: Todo: only works when two phase are stable! 
end


%% Armin: output to the screen

X = ['Composition c1 and c2: ',num2str(Ni(1)),' , ', num2str(Ni(2))];
disp(X)

X = ['x11 and x21 initial guess: ',num2str(xp_hld(1,1,1)),' , ', ...
    num2str(xp_hld(1,1,2))];
disp(X)

X = ['x11 and x21 result: ',num2str(xp_hld(end,1,1)),' , ', ...
    num2str(xp_hld(end,1,2))];
disp(X)

X = ['x12 and x22 initial guess: ',num2str(xp_hld(1,2,1)),' , ', ...
    num2str(xp_hld(1,2,2))];
disp(X)

X = ['x12 and x22 result: ',num2str(xp_hld(end,2,1)),' , ', ...
    num2str(xp_hld(end,2,2))];
disp(X)

X = ['n1 and n2 guess: ',num2str(np_hld(1,1)),' , ', num2str(np_hld(1,1))];
disp(X)

X = ['n1 and n2 result: ',num2str(np_hld(end,1)),' , ', ...
    num2str(np_hld(end,2))];
disp(X)

X = ['squared residue of the last itteration: ', num2str(r_hld(end))];
disp(X)

X = ['difference in squared residues of the last two itterations: ',...
    num2str(del_r)];
disp(X)

X = ['number of itterations: ',num2str(ct1)];
disp(X)
