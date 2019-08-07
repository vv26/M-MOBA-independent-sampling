function pj=tdin(xab,yab,xb,sig,k,n,tau)%%the probability of no change
kp_x=n(k,1)*(n(k,1)+tau)./(tau*sig(k,1).^2);
kp_y=n(k,2)*(n(k,2)+tau)./(tau*sig(k,2).^2);
mu_x=xb(k,1);
mu_y=xb(k,2);
nu_x=n(k,1)-1;
nu_y=n(k,2)-1;
txab=(xab-mu_x)*sqrt(kp_x);
tyab=(yab-mu_y)*sqrt(kp_y);
pj1=(tcdf(txab(2),nu_x)-tcdf(txab(1),nu_x));
pj2=(tcdf(tyab(2),nu_y)-tcdf(tyab(1),nu_y));
pj=[pj1,pj2];
 