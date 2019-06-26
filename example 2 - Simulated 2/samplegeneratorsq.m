function sample=samplegeneratorsq(n)
dd=randsample(2,1);
if n==1
        aa=randsample(2000,1,1)/10000+0.4;
        bb=randsample(1000,1,1)/1000;
        cc=randsample(1000,1,1)/10000;
     
      if dd==1
          sample=[aa bb cc];
      else
          sample=[bb aa cc];
      end      

else
        aa=randsample(2000,1,1)/10000+0.8;
        bb=randsample(2000,1,1)/10000;
        cc=randsample(1000,1,1)/10000;
        if n==2
          sample=[aa bb cc];
        else
          sample=[bb aa cc];
        end   
end
   
