n0=5;
matlabzero=10^(-6);
%Mu=[0.5 5.5;1.9 4.2;2.8 3.3;3 3;3.9 2.1;4.3 1.8;4.6 1.5;3.8 6.3;4.8 5.5;5.2 5;5.9 4.1;6.3 3.8;6.7 7.2;7 7;7.9 6.1;9 9];
Mu=[1.5 4.8;4 6;3 1;1 5;2.8 3;2 6];
objcost = [1,2];% cost 1 to evaluate ojt 1 and 2 to ojb 2
totalcostequal = [];
totalcostds = [];
totalcostpcs = [];
%Mu=[1 2;3 1;5 5];
[r,cl]=size(Mu);tsig=2*ones(r,cl);
for k=1:r 
    a=Mu(k,1)<Mu([1:k-1,k+1:r],1);
    b=Mu(k,2)<Mu([1:k-1,k+1:r],2);
    c=a|b;
    cn=sum(c);
    f0(k)=(cn==(r-1));%% if point k is on the pareto front, f0(k)=1,0 otherwise
end

jn=100;   %%how many different values of budget will be tested (a.k.a how many steps in each iteration)
ct=zeros(3,jn);%%the times of the right choice
T1=zeros(1,jn);%%how many times algorithm used with each budget
T2=zeros(1,jn);
budgets=1;%%the smallest budget
budgeti=1;%%budget step
mre=5000; %% the maximum repitation
record=zeros(r,cl,mre,jn);%% record pj(1:r,1:cl) at each step in each iteration and therefore we can know which alternative gets the new budgets in this iteration
record1=zeros(r,cl,mre,jn);
record2=zeros(r,mre,jn);
for re=1:mre
    for i=1:r
        for j=1:cl
            sps(i,j,1:jn)=normrnd(Mu(i,j),tsig(i,j),jn,1);
        end
    end
%%%%%%%%%%%%%%%the initial sample mean value and variance%%%%%%%%%%%%%%%%%%
    xb0=zeros(r,cl);sig0=zeros(r,cl);
    for i=1:n0
        X(1:r,2*i-1:2*i)= normrnd(Mu,tsig);
        xb0= (xb0+X(1:r,2*i-1:2*i));
    end
    xb0 = (xb0/n0);
    for i=1:n0
        sig0 = (sig0+(X(1:r,2*i-1:2*i)-xb0).^2);
    end
    sig0 =sqrt(sig0/(n0-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%euqal allocation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xb=xb0;n=n0*ones(r,2);budget=0;    
    for j=1:jn
        tbudget=budgets+(j-1)*budgeti;
        while(budget<tbudget)
            mn=mod(budget,(2*r))+1;
            indicator=mod(mn-1,2)+1;
            mn=fix((mn+1)/2);
            budget=budget+1;
            X=sps(mn,indicator,n(mn,indicator)-n0+1);
            xb(mn,indicator)= (n(mn,indicator)*xb(mn,indicator)+X)/(n(mn,indicator)+1);
            n(mn,indicator)=n(mn,indicator)+1;
            costequal(j) = objcost(indicator);
            if j==1
                 totalcostequal(j)= costequal(j);
            else
                totalcostequal(j)= costequal(j)+totalcostequal(j-1);
            end
        end
        for k=1:r
            a=xb(k,1)<xb([1:k-1,k+1:r],1);
            b=xb(k,2)<xb([1:k-1,k+1:r],2);
            c=a|b;
            cn=sum(c);
            f1(k)=(cn==(r-1));
        end
        ct(1,j)=ct(1,j)+(sum(f1==f0)==r);
    end
%%%%%%%%%%%%%%%%allocation by independent objective%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xb=xb0; sig=sig0;n=n0*ones(r,cl);budget=0;%%n is divided into nx and ny
    eqb=0;%%to count how many budgets have been allocated by equal allocation when IDO method failed
    for k=1:r
        a=xb(k,1)<xb([1:k-1,k+1:r],1);
        b=xb(k,2)<xb([1:k-1,k+1:r],2);
        c=a|b;
        cn=sum(c);
        f1(k)=(cn==(r-1));%%if point k is on the pareto front based on five observations, f1(k)=1,0 otherwise
    end
    for j=1:jn
        tbudget=budgets+(j-1)*budgeti;
        while budget<tbudget
            pj= paretotin(xb,sig,f1,r,n,1,objcost);
            record(:,:,re,j)=pj;%%record in each iteration which alternative gets the sampling and the pj accordingly
            if sum(pj)>matlabzero
                [m1,mn1]=max(pj(:,1));%%m is the largest number of pj and mn is the correspoding alternative's sequence number
                [m2,mn2]=max(pj(:,2));
                if m1>m2
                    indicator=1;
                    mn=mn1;
                else
                    indicator=2;
                    mn=mn2;
                end
                recordds(re,j)=mn+100;
                recorddso(re,j)= indicator+100;
            else
                pj=paretotin(xb,sig,f1,r,n,10,objcost);
                record1(:,:,re,j)=pj;
                if sum(pj)>matlabzero
                    [m1,mn1]=max(pj(:,1));%%m is the largest number of pj and mn is the correspoding alternative's sequence number
                    [m2,mn2]=max(pj(:,2));
                    if m1>m2
                        indicator=1;
                        mn=mn1;
                    else
                        indicator=2;
                        mn=mn2;
                    end
                    recordds(re,j)=mn+200;
                    recorddso(re,j)= indicator+200;
                else
                    mn=mod(eqb,(2*r))+1;
                    indicator=mod(mn-1,2)+1;
                    mn=fix((mn+1)/2);
                    eqb=eqb+1;
                    recordds(re,j)=mn+300;
                    recorddso(re,j)= indicator+300;
                end
            end
            X=sps(mn,indicator,n(mn,indicator)-n0+1);
            if indicator==1
                xb(mn,1)=(n(mn,1)*xb(mn,1)+X)/(n(mn,1)+1);
                sig(mn,1)=sqrt((n(mn,1)-1)/n(mn,1)*sig(mn,1).^2+1/(n(mn,1)+1)*(X-xb(mn,1)).^2);
                n(mn,1)=n(mn,1)+1;
            else  
                xb(mn,2)=(n(mn,2)*xb(mn,2)+X)/(n(mn,2)+1);
                sig(mn,2)=sqrt((n(mn,2)-1)/n(mn,2)*sig(mn,2).^2+1/(n(mn,2)+1)*(X-xb(mn,2)).^2);
                n(mn,2)=n(mn,2)+1;
            end
            budget=budget+1;
            costds(re,j) = objcost(indicator);
            if j==1
                 totalcostds(re,j)= costds(re,j);
            else
                totalcostds(re,j)= costds(re,j)+totalcostds(re,j-1);
            end
            for k=1:r 
                a=xb(k,1)<xb([1:k-1,k+1:r],1);
                b=xb(k,2)<xb([1:k-1,k+1:r],2);
                c=a|b;
                cn=sum(c);
                f1(k)=(cn==(r-1));
            end
        end
%         ct(2,j)=ct(2,j)+(sum(f1==f0)==r);
        mobadscs(re,j)= double(sum(f1==f0)==r); % correct selection in each iteration 
    end
    %%%%%%%%%%%%%%%%allocation by M-MOBA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xb=xb0; sig=sig0;n=n0*ones(r,1);budget=0;nonalgorithm=0;    
    for k=1:r
        a=xb(k,1)<xb([1:k-1,k+1:r],1);
        b=xb(k,2)<xb([1:k-1,k+1:r],2);
        c=a|b;
        cn=sum(c);
        f1(k)=(cn==(r-1));%%if point k is on the pareto front based on five observations, f1(k)=1,0 otherwise
    end
    for j=1:jn/2
        tbudget=budgets+(j-1)*budgeti;
        while(budget<tbudget)
            pj= (paretot(xb,sig,f1,r,n,1));
            record2(:,re,j)=pj;%%record in each iteration which alternative gets the sampling and the pj accordingly
            if sum(pj)~=0
                [m,mn]=max(pj);%%m is the largest number of pj and mn is the correspoding alternative's sequence number 
                %[m,mn]= min(pj); 
                T1(1,j)=T1(1,j)+1;%%how many times use algorithm
                recordmoba(re,j)=mn+100;
            else
                pj=paretot(xb,sig,f1,r,n,10);
                if sum(pj)~=0
                    [m,mn]=max(pj);
                    T1(1,j)=T1(1,j)+1;
                    recordmoba(re,j)=mn+200;
                else
                    T2(1,j)=T2(1,j)+1;
                    mn=mod(nonalgorithm,r)+1;
                    nonalgorithm=nonalgorithm+1;
                    recordmoba(re,j)=mn+300;
                end
            end
            
            X=sps(mn,:,n(mn)-n0+1);%
            sig(mn,1:2)=(sqrt((n(mn)-1)/n(mn)*sig(mn,1:2).^2+1/(n(mn)+1)*(X-xb(mn,1:2)).^2));
            xb(mn,1:2)=((n(mn)*xb(mn,1:2)+X)/(n(mn)+1));
            n(mn)=n(mn)+1;
            budget=budget+1;% 
            for k=1:r 
                a=xb(k,1)<xb([1:k-1,k+1:r],1);
                b=xb(k,2)<xb([1:k-1,k+1:r],2);
                c=a|b;
                cn=sum(c);
                f1(k)=(cn==(r-1));
            end 
        end
        ct(3,j)=ct(3,j)+(sum(f1==f0)==r);
    end

    fprintf('re=%d\n',re);
end

x=(1:jn)*budgeti;
costpcs=(1:jn/2)*budgeti*3;% MMOBA PCS evaluate two objs in each iteration
cp=ct/mre;
cs =zeros(1,jn*1.5);
for cost =1:(jn*1.5)
    a = unique(cost);
    frequent(cost) = histc(totalcostds(:),a); % how many times cost = a appears
    for row =1:mre
        if mobadscs(row,totalcostds(row,:)==cost) == 1
            cs(cost) = cs(cost)+1;
        end
    end
end
x = 1:jn*1.5;
pcsds = cs./frequent; % how many times is correctly selected in each cost
for i = 1:jn*1.5
    if isnan(pcsds(i))
        x(i) = NaN;
    end
end
pcsds(isnan(pcsds)) = []; % remove nan number
x(isnan(x)) = [];% remove the x accordingly

figure
plot(totalcostequal,cp(1,:),'r')
legend('Equal')
xlabel('cost')
ylabel('P(CS)')
xlim([0 750])

figure
plot(x,pcsds,'-.b')
legend('M-MOBA DS PCS')

figure
plot(costpcs,cp(3,1:jn/2),'--b')
legend('M-MOBA PCS')

figure
plot(totalcostequal,cp(1,:),'r', x,pcsds,'-.y',costpcs,cp(3,1:jn/2),'-b')
legend('Equal','M-MOBA DS PCS','M-MOBA PCS')
xlabel('cost')
ylabel('P(CS)')
xlim([0 750])

