function pj=paretotin(xb,sig,f1,r,n,tau,objcost)
pj=zeros(r,2);
for k=1:r
    if f1(k)==1
        f2=zeros(1,r);
        for i=[1:k-1,k+1:r]
            a=(xb(i,1)<=xb([1:k-1,k+1:r],1));
            b=(xb(i,2)<=xb([1:k-1,k+1:r],2));
            c=a|b;
            cn=sum(c);
            f2(i)=(cn==(r-1));
        end
    else
        f2=f1;
    end
    A=sort(xb(f2==1,1));
    B=sort(xb(f2==1,2),'descend');

    if f1(k)==1
        A2=[-inf;A;inf]; %%situation 3.1.1 in paper
        B2=[inf;B;-inf];
        xab(1)=max(A2(A2<xb(k,1)));yab(1)=max(B2(B2<xb(k,2)));
        xab(2)=min(A2(A2>xb(k,1)));yab(2)=min(B2(B2>xb(k,2)));
    else
        xab(1)=min(A(B<xb(k,2)));yab(1)=min(B(A<xb(k,1)));
        xab(2)=inf;yab(2)=inf;
    end
    pj(k,1:2)=tdin(xab,yab,xb,sig,k,n,tau);%%here%%the probability of no change
end
pj=(1-pj)./objcost;