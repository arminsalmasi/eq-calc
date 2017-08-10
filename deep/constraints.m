function [g,dg,d2g]=constraints(nel,nph,nphtot,Ni,status,xp,np)

% nel is number of elements
% nph is number of phases currently assumed to be stable
% nphtot is the total number of phases in the system
% Ni is the number of moles of each component (equilibrium conditions)
% status specifies whether a phase is assumed stable (=1), or not (=0)
% xp contain the phase compositions
% np is the number of moles of each phase

% Evaluate:
% the constraints g(nel,1)
% the 1st derivatives of the constraints dg(nel,nel*nph)
% the 2nd derivatives of the constraints d2g(nel,(n^2+n)/2) n=nel*nph
%     d2g is for each row given in order 
%     d2g/dy{1}dy{1} d2g/dy{2}dy{1} d2g/dy{2}dy{2} d2g/dy{3}dy{1} etc

g=zeros(nel,1);
dg=zeros(nel,nel*nph);
d2g=zeros(nel,((nel*nph)^2+nel*nph)/2);

for i=1:nel
    g(i)=Ni(i);
    for j=1:nphtot
        if (status(j)==1)
            g(i)=g(i)-np(j)*xp(j,i);
        end
    end
end

for i=1:(nel-1)
    k=1;
    for j=1:nphtot
        if (status(j)==1)
            dg(i,(k-1)*nel+i)=-np(j);
            dg(i,(k-1)*nel+nel)=-xp(j,i);
            k=k+1;
        end
    end
end

i=nel;
k=1;
for j=1:nphtot
    if (status(j)==1)
        for m=1:nel-1
            dg(i,(k-1)*nel+m)=np(j);
        end
        dg(i,(k-1)*nel+nel)=-xp(j,nel);
        k=k+1;
    end
end
	   
m=1;
r=1;%phase counter for index i
s=1;%variable counter for index i
for i=1:nel*nph
    p=1;%phase counter for index j
    t=1;%variable counter for index j
    for j=1:i
        for k=1:nel
            % get derivative d^2g_k/dz_jdz_i where z is xp and np
            % d2g(k,m)
            if (k==nel)
                % derivatives of the dependent element
                if (p==r)
                    if (s==nel)
                        if (t<nel)
                            d2g(k,m)=1.0;
                        end
                    end
                    if (t==nel)
                        if (s<nel)
                            d2g(k,m)=1.0;
                        end
                    end
                end
            else
                % derivatives of non-dependent elements
                if (p==r)
                    if (s==nel)
                        if (k==t)
                            d2g(k,m)=-1.0;
                        end
                    end
                    if (t==nel)
                        if (k==s)
                            d2g(k,m)=-1.0;
                        end
                    end
                end
            end
        end
        m=m+1;
        t=t+1;
        if (t>nel)
            p=p+1;
            t=1;
        end
    end
    s=s+1;
    if (s>nel)
        r=r+1;
        s=1;
    end
end

end