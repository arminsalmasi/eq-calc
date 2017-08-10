clear variables

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
Ni=1*[0.65 0.35];
Ntot=sum(Ni);   % total number of moles of components, and formula unit 
                % given with decimals summing up to 1


% 1. Make an initial guess of the set of variables and the variable values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nph=nphtot;
status=[1 1];
xp=[0.3 0.7; 0.5001 0.4999]; % first row phase a, sec row phase b, col el 1, 2 
np=[0.4999*Ntot 0.5001*Ntot];
y=xp(:,1);

lambda=[1e-5;1e-5];
x=[xp(1,1); np(1); xp(2,1); np(2); lambda(1); lambda(2)];


% 2. Convergence loop start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numiter=1;
np_iter(1,:)=np(:);
xp_iter(1,:)=[xp(1,:), xp(2,:)];
while 1
    
% 3. Construct the system of linear equations from the present state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Get constraints and derivatives
[g,dg,d2g]=constraints(nel,nph,nphtot,Ni,status,xp,np);
%   Get Gibbs energy for phases and derivatives 
for phindex=1:nph
    [G(phindex),dG(phindex,:),d2Gtemp(phindex,:),mu(phindex,:)]=...
        binary_regular_solution(phindex,y(phindex),np(phindex),T);
end


%   Reform d2G matrix
d2G=[d2Gtemp(1,1) d2Gtemp(1,2) d2Gtemp(1,3) 0 0 d2Gtemp(2,1) 0 0 d2Gtemp(2,2) d2Gtemp(2,3)];

%    Insert d2g into A
n=nel*nph;
A=zeros(n+nel,n+nel);
k=1;
for i=1:n
    i
    for j=1:i
        j
        A(i,j)=d2g(1, k)*lambda(1) + d2g(2, k)*lambda(2) + d2G(k);
        A(j,i)=A(i,j);
        k=k+1
    end
end
%   Insert dg into A
for i=1:n
    k=1;
    for j=n+1:n+nel
        A(i,j)=dg(k,i);
        A(j,i)=dg(k,i);
        k=k+1;
    end
end

%   Construct b vector
b=-[(dG(1,:) + lambda(1)*dg(1,1:2) + lambda(2)*dg(2,1:2) )';...
    (dG(2,:) + lambda(1)*dg(1,3:4) + lambda(2)*dg(2,3:4) )';...
    g];

% 4. Scale the equations if necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: This is not strictly necessary.


% 5. Solve for updates of the variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note: The default "slash" solver did not converge. Use, for example, qmr
[dx, FLAG(numiter,:),RELRES(numiter,:),ITER(numiter,:),RESVECtemp] = qmr(A,b,1e-6,50);
RESVEC{numiter,:}=RESVECtemp;
if imag(dx)
    dx
    numiter
    return
end

% 6. Do a line search and select an optimum step length for the variable update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: This is not strictly necessary.
% Note: Check that variables do not obtain un-physical values
% Note: Also update the dependent phase composition variables
theta=[-logspace(1,-20,1000), logspace(-20,1,1000)];
r=sum((A*(x+dx*theta)).^2);
j=1;
while j<=length(theta)
    for i=1:nel*nph
        if mod(i,nel)~=0
            if x(i)+theta(j)*dx(i)>1
                theta(j)=[];
                r(j)=[];
                j=j-1;
                break
            elseif x(i)+theta(j)*dx(i)<0
                theta(j)=[];
                r(j)=[];
                j=j-1;
                break
            end
        else
            if x(i)+theta(j)*dx(i)>Ntot
                theta(j)=[];
                r(j)=[];
                j=j-1;
                break
            elseif x(i)+theta(j)*dx(i)<0
                theta(j)=[];
                r(j)=[];
                j=j-1;
                break
            end
        end
    end
    j=j+1;
end

theta=theta(r==min(r));
if length(theta)<1
    numiter
    x
    dx
    A
    return
end
theta=theta(1);

% 7. Update the variable values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=x+theta*dx;

xp=[x(1) 1-x(1); x(3) 1-x(3)];
np=[x(2) x(4)];
y=xp(:,1);
lambda=[x(5); x(6)];

np_iter(numiter+1,:)=np(:);
xp_iter(numiter+1,:)=[xp(1,:), xp(2,:)];
dx_iter(numiter)=norm(dx);
r_iter(numiter)=min(r);
b_iter(numiter)=norm(b);
% 8. Check if the set of variables should be changed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for phindex=1:nph
%     [G(phindex),dG(phindex,:),d2Gtemp(phindex,:),mu(phindex,:)]=...
%         binary_regular_solution(phindex,y(phindex),np(phindex),T);
%     Gm(phindex)=G(phindex)/np(phindex);
%     dGm(phindex)=dG(phindex,1)/np(phindex);
% end
% G(1)+dGm(1)*(xp(1,2));

nph=nphtot;
status=[1 1];

% 9. Check convergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if min(r)<100 || numiter>=200
    break
end
numiter=numiter+1;
end

min_r=min(r)
numiter
np

i=1;
figure(i)
clf
plot(0:numiter,np_iter)
xlabel('Number of iterations')
legend('Phase a','Phase b')
ylim([0 1*Ntot])

i=i+1;
figure(i)
clf
plot(0:numiter,xp_iter)
xlabel('Number of iterations')
legend('X(1) in a','X(2) in a','X(1) in b','X(2) in b')
ylim([0 1])

i=i+1;
figure(i)
clf
plot(1:numiter,dx_iter)
title('Norm of dx')
xlabel('Number of iterations')

i=i+1;
figure(i)
clf
plot(1:numiter,r_iter)
title('Square residual')
xlabel('Number of iterations')

i=i+1;
figure(i)
clf
plot(1:numiter,b_iter)
title('Norm of b')
xlabel('Number of iterations')