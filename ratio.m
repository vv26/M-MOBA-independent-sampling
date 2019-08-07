s2 = 0;
s3= 0;
s1 =0;
s4 =0;
for re=1:mre
    for j=1:jn
        if max(record(:,1,re,j))<max(record(:,2,re,j))
            s2 = s2+1;
        else
            s1 = s1+1;
        end
        if max(record1(:,1,re,j))<max(record1(:,2,re,j))
            s4 = s4+1;
        else
            s3 = s3+1;
        end
        
    end
end
S1 =s1+s3;
S2 =s2+s4;
R=S1/S2;