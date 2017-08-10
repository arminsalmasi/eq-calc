clear variables; close all; clc;
% Test implementation of an equilibrium calculator
%
% Equilibrium conditions are P, T and Ni
% Only single-sublattice phases are considered
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
Ntot =sum(Ni); % this is a constraint 
% 1. Make an initial guess of the set of variables and the variable values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialize holder variables
% nph is number of phases currently assumed to be stable
nph_hld(1) = nphtot_hld(1); %statrting guess - both phases are statble
% status specifies whether a phase is assumed stable (=1), or not (=0)
status_hld(1,:) = [1,1]; %statrting guess - both phases are statble
% xp contain the phase compositions
xp_hld(1,:,:) = [0.49 0.51; 0.51 0.49]; % Initial - 2 statble phase
% np is the number of moles of each phase
np_hld(1,:) = [0.49*Ntot 0.51*Ntot] ; % Initial - 2 statble phase
% lmb = lagrangian multipliers here it is a 1*2 vector(2 constraints)
lmb_hld(1,:) = [1e-5 1e-5]; %% Initial guess - both phases are statble
                           
% 2. Convergence loop start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt = 1;
while 1==1
  %% initialize step variables form the holder variables
  xp(:,:) = xp_hld(cnt,:,:);
  np(:)   = np_hld(cnt,:);
  status  = status_hld(cnt,:);
  nph = nph_hld(cnt);
  nphtot = nph_hld(cnt);
  lmb(:) = lmb_hld(cnt,:);

  %% Calculate nGm of both pghases and derivitives of Gibbs free energy 
  for index = 1 : nph 
    [G(index),dG(index,:),d2G(index,:),mu(index,:)] = ...
           binary_regular_solution(index, xp(index, 1), np(index), T);
  end
  %%  Calculate constraints and their derivitives regards to all variables
  [g,dg,d2g] = constraints(nel,nph,nphtot,Ni,status,xp,np);
  
  % 3. Construct the system of linear equations from the present state
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % transform d2G to 4*4 system % d2G(y1,y2)=0,
  % d2G(y1,n2)=0,d2G(n1,y2)=0 , d2G(n1,n2) =0 
  % Todo : this must be automatic - look at the G function 
  d2G_A(:) = ...
           [d2G(1,1) d2G(1,2) d2G(1,3) 0 0 d2G(2,1) 0 0 d2G(2,2) d2G(2,3)];

  %% construct A
  nd = nel*nph;  %% dimension - 1 or cnt is it always 4*4?
  A(:,:)= [zeros(nd,nd),(dg(:,:))';dg(:,:),zeros(nel,nel)];
  cnt2=1;
  for nd=1:nd
    for cnt_theta=1:nd
      A(nd,cnt_theta)=d2g(1, cnt2)*lmb(1) + d2g(2, cnt2)*lmb(2) +...
                      d2G_A(cnt2);
      A(cnt_theta,nd)=A(nd,cnt_theta);
      cnt2=cnt2+1;
    end
  end

  c = -[(dG(1,:) + lmb(1)*dg(1,1:2) + lmb(2)*dg(2,1:2) )';...
        (dG(2,:) + lmb(1)*dg(1,3:4) + lmb(2)*dg(2,3:4) )';...
        g];
  delta = zeros(6,1);
  
  % 4. Scale the equations if necessary
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Note: This is not strictly necessary.

  % 5. Solve for updates of the variables
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Note: The default "slash" solver did not converge. Use, for example, qmr
% X = qmr(A,B,TOL,MAXIT) specifies the maximum number of iterations. If MAXIT is [] then qmr uses the default, min(N,20).
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
% [X,FLAG,RELRES,ITER,RESVEC] = qmr(A,B,...) also returns a vector of the
% residual norms at each iteration, including NORM(B-A*X0).

 %test=A/c';
 [delta,FLAG,RELRES,ITER,RESVEC] = qmr(A,c,[],1e6);
 if not(isreal(delta))
   err = 'complex resolts for delta'
   break % out of the convergance loop(while)
 end
 vars = [xp(1,1); np(1); xp(2,1); np(2); lmb(1); lmb(2)];


 % 6. Do a line search and select an optimum step length for the variable update
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Note: This is not strictly necessary.
 % Note: Check that variables do not obtain un-physical values
 % Note: Also update the dependent phase composition variables
 %% there are better ways of doing this! this is not a line search:
 %  this just calculate r for different values of theta and remove the ones
 %  which wiolate the constrains and then selects the one whihc results in
 %  r minimum
 theta=[-logspace(1,-10,10000), logspace(-10,1,10000)];
 r=sum((A*(vars+delta*theta)).^2);
 cnt_theta=1; %% itteration counter on lenght of theta
 %% remove the thetas which violate the conditions
 while cnt_theta<=length(theta) 
   for nd =1:2:nel*nph
     if vars(nd)+theta(cnt_theta)*delta(nd)>1 || vars(nd)+theta(cnt_theta)*delta(nd)<0 %check condition on xp <xp<1
       theta(cnt_theta)=[];
       r(cnt_theta)=[];
       cnt_theta=cnt_theta-1;
       break %% to the for loop
     end
     if vars(nd+1)+theta(cnt_theta)*delta(nd+1)>Ntot || vars(nd+1)+theta(cnt_theta)*delta(nd+1)<0 %check condition 0<np<Ntot=1 
       theta(cnt_theta)=[];
       r(cnt_theta)=[];
       cnt_theta=cnt_theta-1;
       break %% to the for loop
     end
   end
   cnt_theta=cnt_theta+1; 
 end
 %% find the theta which gives minimum of r in the remaining thetas
 theta=theta(r==min(r)); 
 if length(theta)<1 %% error handler
   err = 'all thetas in the range violated xp and np constraints'  
   return %% to command line
 end 
 theta = theta(1); %% in case more than one thetas meet the conditions
  % 7. Update the variable values
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  vars=vars+theta*delta;
  xp=[vars(1) 1-vars(1); vars(3) 1-vars(3)];
  np=[vars(2) vars(4)];
  lmb=[vars(5); vars(6)];

  %% save values of the current iteration (for tracking)
  for nd = 1: 2
    G_hld(cnt, nd)= G(nd);
    dG_hld(cnt,nd,:)= dG(nd,:);
    d2G_hld(cnt,nd,:)= d2G(nd,:);
    mu_hld(cnt, nd,:)= mu(nd,:);
  end
  g_hld(:,cnt) = g(:);
  dg_hld(:,:,cnt) = dg(:,:);
  d2g_hld(:,:,cnt) = d2g(:,:);
  d2GtoA_hld(cnt,:) = d2G_A(:);   
  delta_hld(:,cnt)= delta;
  r_hld(cnt) = min(r);
  c_hld(:,cnt) = c;
  %% initialize variables for next itteration step
  lmb_hld(cnt+1,:) =lmb(:) ;
  status_hld(cnt+1,:) = status(:);
  nph_hld(cnt+1) = nph ;
  nphtot_hld(cnt+1) = nphtot  ;
  xp_hld(cnt+1,:,:) = xp(:,:) ;
  np_hld(cnt+1,:)= np(:);
   
  % 9. Check convergence
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if min(r)<1000 || cnt>1000 %% 
     break  % while to the main %% end of itteration conditions meet
  end
   
  %% update counter variable
  cnt =cnt+1;            



  % 8. Check if the set of variables should be changed
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %% in this composition always two phase are stable! what if only one of 
  %  the two is stable?! 
end

X = ['Composition c1 and c2: ',num2str(Ni(1)),' , ', num2str(Ni(2))];
disp(X)


X = ['x11 and x21 initial guess: ',num2str(xp_hld(1,1,1)),' , ', num2str(xp_hld(1,1,2))];
disp(X)

X = ['x11 and x21 result: ',num2str(xp_hld(end,1,1)),' , ', num2str(xp_hld(end,1,2))];
disp(X)

X = ['x12 and x22 initial guess: ',num2str(xp_hld(1,2,1)),' , ', num2str(xp_hld(1,2,2))];
disp(X)

X = ['x12 and x22 result: ',num2str(xp_hld(end,2,1)),' , ', num2str(xp_hld(end,2,2))];
disp(X)


X = ['n1 and n2 guess: ',num2str(np_hld(1,1)),' , ', num2str(np_hld(1,1))];
disp(X)

X = ['n1 and n2 result: ',num2str(np_hld(end,1)),' , ', num2str(np_hld(end,2))];
disp(X)

X = ['squared residue of the last itteration: ', num2str(r_hld(end))];
disp(X)

X = ['number of itterations: ',num2str(cnt)];
disp(X)
