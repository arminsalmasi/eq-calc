function [G,dG,d2G,mu]=binary_regular_solution(index,y,n,T)

% For phase 'index' evaluate
% Gibbs energy 'G'
% first derivatives 'dG' in order dG/dy dG/dn
% second derivatives 'd2G' in order
%
% d2G/dy2 d2G/dydn d2G/dn2
%
% Chemical potentials mu
%
% The independent composition variable is given by scalar 'y'
% The number of moles formula units are given by scalar 'n'

R=8.31451;
RT=R*T;

mu=zeros(2,1);

if (index==1)
%-------------------phase 1='a'-----------------
G0a1=0.0;
G0a2=10000.0;
La12=-10000.0;
ydep=1.0-y;

Gm=y*G0a1+ydep*G0a2+RT*(y*log(y)+ydep*log(ydep))+y*ydep*La12;

dGmdx=G0a1-G0a2+RT*(log(y)-log(ydep))+La12*(ydep-y);

d2Gmdx2=RT*((1.0/y)+(1.0/ydep))-2.0*La12;

G=n*Gm;

dG=zeros(2,1);
dG(1)=n*dGmdx;
dG(2)=Gm;


d2G=zeros(3,1);

d2G(1)=n*d2Gmdx2;
d2G(2)=dGmdx;
d2G(3)=0.0;

mu(1)=G0a1+RT*log(y)+0.5*ydep*La12;
mu(2)=G0a2+RT*log(ydep)+0.5*y*La12;

elseif (index==2)
%-------------------phase 2='b'-----------------
G0b1=10000.0;
G0b2=0.0;
Lb12=-20000.0;
ydep=1.0-y;

Gm=y*G0b1+ydep*G0b2+RT*(y*log(y)+ydep*log(ydep))+y*ydep*Lb12;

dGmdx=G0b1-G0b2+RT*(log(y)-log(ydep))+Lb12*(ydep-y);

d2Gmdx2=RT*((1.0/y)+(1.0/ydep))-2.0*Lb12;

G=n*Gm;

dG=zeros(2,1);
dG(1)=n*dGmdx;
dG(2)=Gm;


d2G=zeros(3,1);

%d2G/dx2
d2G(1)=n*d2Gmdx2;

%d2G/dxdn
d2G(2)=dGmdx;

%d2G/dn2
d2G(3)=0.0;

mu(1)=G0b1+RT*log(y)+0.5*ydep*Lb12;
mu(2)=G0b2+RT*log(ydep)+0.5*y*Lb12;

end

end