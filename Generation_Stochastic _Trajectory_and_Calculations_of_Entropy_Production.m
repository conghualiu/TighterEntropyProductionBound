%Generation of Stochastic Trajectory and Calculation of Entropy Production
clear all;
clc;
qhcs=1;
ATE=zeros(qhcs,10);FD=zeros(qhcs,10);
tf1=1;tf2=1;zgjs=20000;tf=2;
bh=2;bc=1;dh=4;dc=2;
kh=exp(-1/3);kc=2-2*exp(-1/3);
ah=bh*kh;ac=bc*kc;
wc=0.2;d=0.2;

for abc=1:1:qhcs

%following is the initial state probability
fi00=0.6*(1-wc);fi10=0.4*wc;fi01=0.6*wc;fi11=0.4*(1-wc);

%folowing is the final state probability
%p00=0.36;p10=0.24;p01=0.24;p11=0.16;%d=0.2,wc=0
%p00=0.348;p10=0.232;p01=0.252;p11=0.168;%d=0.2,wc=0.1
p00=0.336;p10=0.22399;p01=0.263998;p11=0.176;%d=0.2,wc=0.2
%p00=0.324;p10=0.216;p01=0.276;p11=0.184;%d=0.2,wc=0.3
%p00=0.312;p10=0.208;p01=0.288;p11=0.192;%d=0.2,wc=0.4
%p00=0.3;p10=0.2;p01=0.3;p11=0.2;%d=0.2,wc=0.5
%p00=0.2880;p10=0.1920;p01=0.3120;p11=0.2080;%d=0.2,wc=0.6
%p00=0.2760;p10=0.1840;p01=0.3240;p11=0.2160;%d=0.2,wc=0.7
%p00=0.2640;p10=0.1760;p01=0.3360;p11=0.2240;%d=0.2,wc=0.8
%p00=0.2520;p10=0.1680;p01=0.3480;p11=0.2320;%d=0.2,wc=0.9
%p00=0.24;p10=0.16;p01=0.36;p11=0.24;%d=0.2,wc=1



%p00=0.336;p10=0.22399;p01=0.263998;p11=0.176;%d=0.2,wc=0.2
%d change, the initial probability in the feedback process is not changed





ate1=0;ate2=0;ate3=0;ate4=0;ate5=0;ate6=0;ate7=0;ate9=0;ate10=0;ate8=0;
fd1=0;fd2=0;fd3=0;fd4=0;fd5=0;fd6=0;fd7=0;fd9=0;fd10=0;fd8=0;
KK1=zeros(1,20000);KK2=zeros(1,20000);KK3=zeros(1,20000);KK4=zeros(1,20000);KK5=zeros(1,20000);KK6=zeros(1,20000);KK7=zeros(1,20000);KK8=zeros(1,20000);

syms t3 t2;
f=piecewise(t3<1,dh*(kh)^d,t3>=1,dh*kh)+piecewise(t3<1,dc*(kc)^d,t3>=1,dc*kc);
f1=piecewise(t3<1,dh*(kh)^d,t3>=1,dh*kh);
f2=piecewise(t3<1,dc*(kc)^d,t3>=1,dc*kc);
%note that the elements in the zhu dui jiao yuan does not change to negative
%w000=a0*(2+3/4*r*sin(2*pi*t/tm));w001=a0*(1-r/4*cos(2*pi*t/tm));w002=a0*(1-r*sin(2*pi*t/tm));
%w010=a0*(1+r*sin(2*pi*t/tm));w011=a0*(2+3/4*r*cos(2*pi*t/tm));w012=a0*(1+r/4*sin(2*pi*t/tm));
%w020=a0*(1-r/4*sin(2*pi*t/tm));w021=a0*(1+r*cos(2*pi*t/tm));w022=a0*(2-3/4*r*sin(2*pi*t/tm));
%w100=a0*(2-3/4*r*sin(2*pi*t/tm));w101=a0*(1+r/4*cos(2*pi*t/tm));w102=a0*(1+r*sin(2*pi*t/tm));
%w110=a0*(1-r*sin(2*pi*t/tm));w111=a0*(2-3/4*r*cos(2*pi*t/tm));w112=a0*(1-r/4*sin(2*pi*t/tm));
%w120=a0*(1+r/4*sin(2*pi*t/tm));w121=a0*(1-r*cos(2*pi*t/tm));w122=a0*(2+3*r/4*sin(2*pi*t/tm));
%random number a, which satisfies a=e^{-\int b(t)dt}

parfor lk=1:zgjs
    
if lk<zgjs*fi00+1
    
        t0=0;T=zeros(1,2000);A=zeros(1,2000);A(1,1)=0;B=zeros(1,2000);

R=rand(2,2000);
W0=[ah+ac,bh+bc;ah+ac,bh+bc];

W01=[f,dh+dc;f,dh+dc];

a=1;

  for m=1:1:2000
    r1=R(1,m);r2=R(2,m);
    t1=vpasolve(-vpaintegral(W0(a,a),t3,t0,t2)-log(r1)==0,t2);
     t0=t1;
    if t0>tf
        m1=m+1;
        break
    end
    
    T(1,m+1)=t0;
    if a==1
        if r2*W0(a,a)<ah
            b=1;
        else
            b=2;
        end
    else
        if r2*W0(a,a)<bh
            b=1;
        else
            b=2;
        end
    end
    
    if a==1
        a=2;
    else 
        a=1;
    end
    A(1,m+1)=a-1;B(1,m)=b;
  end
 
  A1=zeros(1,m1);T1=zeros(1,m1);A1(1,m1)=a-1;T1(1,m1)=tf;B1=zeros(1,m1-2);
  for m=1:m1-1
      A1(1,m)=A(1,m);T1(1,m)=T(1,m);
  end
  for m=1:1:m1-2
      B1(1,m)=B(1,m);
  end

fp1=exp(-vpaintegral(W0(a,a),t3,T1(1,m1-1),tf));
fp3=exp(-vpaintegral(W01(a,a),t3,T1(1,m1-1),tf));

for m=1:m1-2
    fp1=fp1*R(1,m);
  fp3=fp3*exp(-vpaintegral(W01(A1(1,m)+1,A1(1,m)+1),t3,T1(1,m),T1(1,m+1)));
  
end
%fp1:y=0 forward probability fp2:y=0 backward probability fp3:y=1 forward probability
%fp8:y=1 backward probability %fp4:y=0 dual forward probability fp5:y=0 dual backward probability 
%fp6:y=1 dual forward probability fp7:y=1 dual backward probability


fp2=fp1;
fp4=fp1;
fp5=fp1;
fp8=fp3;
fp6=fp3;
fp7=fp3;




  BE1=zeros(1,m1-2); BE2=zeros(1,m1-2);BE31=zeros(1,m1-2);BE32=zeros(1,m1-2);
  BE41=zeros(1,m1-2);BE42=zeros(1,m1-2);BE51=zeros(1,m1-2);BE52=zeros(1,m1-2);
  BE8=zeros(1,m1-2);
 %BE61=zeros(1,m1-2);BE62=zeros(1,m1-2);BE71=zeros(1,m1-2);BE72=zeros(1,m1-2);
  
  for k=2:m1-1
      
          %W2:y=0 transition rate  W21:y=1 transition rate
         %W3: y=0 transition rate  W4:y=1 transition rate
      

         %BE1:y=0 forward transition probability (not include staying) BE2: y=0 backward transition probability (not include staying)
         %BE51:y=1 forward transition probability (not include staying) BE52: y=1 backward transition probability (not include staying)
         
         %BE31: y=0 dual forward transition probability    BE32: y=0 dual backward transition probability
         %BE41:y=1 dual forward transition probability     BE42: y=0 dual backward transition probability
         
         
         u2=[ah+ac,bh+bc;ah+ac,bh+bc];
         u21=[subs(f,t3,T1(1,k)),dh+dc;subs(f,t3,T1(1,k)),dh+dc];
         
         
          W32=[-u2(1,1),u2(1,2);1,1];
         W33=[0;fi00+fi10];
         W34=W32\W33;

         W42=[-u21(1,1),u21(1,2);1,1];
         W43=[0;fi01+fi11];
         W44=W42\W43;

         
          if B(1,k-1)==1
              W2=[ah,bh;ah,bh];
              W21=[subs(f1,t3,T1(1,k)),dh;subs(f1,t3,T1(1,k)),dh];

              W3=[W2(1,1),W2(2,1)*W34(1,1)/W34(2,1);W2(1,2)*W34(2,1)/W34(1,1),W2(2,2)];
              W4=[W21(1,1),W21(2,1)*W44(1,1)/W44(2,1);W21(1,2)*W44(2,1)/W44(1,1),W21(2,2)];
          else
              W2=[ac,bc;ac,bc];
              W21=[subs(f2,t3,T1(1,k)),dc;subs(f2,t3,T1(1,k)),dc];
              W3=[W2(1,1),W2(2,1)*W34(1,1)/W34(2,1);W2(1,2)*W34(2,1)/W34(1,1),W2(2,2)];
              W4=[W21(1,1),W21(2,1)*W44(1,1)/W44(2,1);W21(1,2)*W44(2,1)/W44(1,1),W21(2,2)];
          end

          
          BE1(1,k-1)=W2(A1(1,k)+1,A1(1,k-1)+1);
          BE2(1,k-1)=W2(A1(1,k-1)+1,A1(1,k)+1);
          BE51(1,k-1)=W21(A1(1,k)+1,A1(1,k-1)+1);
          BE52(1,k-1)=W21(A1(1,k-1)+1,A1(1,k)+1);
         

          BE31(1,k-1)=W3(A1(1,k)+1,A1(1,k-1)+1);
          BE32(1,k-1)=W3(A1(1,k-1)+1,A1(1,k)+1);
          BE41(1,k-1)=W4(A1(1,k)+1,A1(1,k-1)+1);
          BE42(1,k-1)=W4(A1(1,k-1)+1,A1(1,k)+1);
         
  end
  
fp1=fi00*fp1;fp3=fi01*fp3;fp4=fi00*fp4;fp6=fi01*fp6;
if A(1,m1-1)==0
    fp2=p00*fp2;fp5=fp5*p00;
    fp7=p01*fp7;fp8=fp8*p01;
    
    MI=p00*(fi00+fi01)/((p00+p01)*fi00);
else
    fp2=fp2*p10;fp5=p10*fp5;
    fp7=p11*fp7;fp8=fp8*p11;
    
    MI=p10*(fi00+fi01)/((p10+p11)*fi00);
end
for m=1:m1-2
    fp1=fp1*BE1(1,m);
    fp2=fp2*BE2(1,m);
    fp3=fp3*BE51(1,m);
    fp8=fp8*BE52(1,m);
    

    fp4=fp4*BE31(1,m);
    fp5=fp5*BE32(1,m);
    fp6=fp6*BE41(1,m);
    fp7=fp7*BE42(1,m);
    
end
%ate1: toal           ate2:excess     ate3:housekeeping 
%ate4: marginal excess                ate5: marginal housekeeping  
%ate6: mixed excess                   ate7:mixed housekeeping  
te1=log(fp1/fp2);
te2=log(fp1/fp5);
te3=log(fp1/fp4);
te4=log((fp1+fp3)/(fp5+fp7));
te5=log((fp1+fp3)/(fp4+fp6));
te6=(log(fp1/fp5)-log((fp1+fp3)/(fp5+fp7)));
te7=(log(fp1/fp4)-log((fp1+fp3)/(fp4+fp6)));
te9=log((fp1+fp3)/(fp2+fp8));
te10=(log(fp1/fp2)-log((fp1+fp3)/(fp2+fp8)));


ate1=ate1+1/zgjs*log(fp1/fp2);
ate2=ate2+1/zgjs*log(fp1/fp5);
ate3=ate3+1/zgjs*log(fp1/fp4);
ate4=ate4+(1/zgjs)*log((fp1+fp3)/(fp5+fp7));
ate5=ate5+(1/zgjs)*log((fp1+fp3)/(fp4+fp6));
ate6=ate6+1/zgjs*(log(fp1/fp5)-log((fp1+fp3)/(fp5+fp7)));
ate7=ate7+1/zgjs*(log(fp1/fp4)-log((fp1+fp3)/(fp4+fp6)));
ate9=ate9+(1/zgjs)*log((fp1+fp3)/(fp2+fp8));
ate10=ate10+(1/zgjs)*(log(fp1/fp2)-log((fp1+fp3)/(fp2+fp8)));



fd1=fd1+1/zgjs*exp(-te1);
fd2=fd2+1/zgjs*exp(-te2);
fd3=fd3+1/zgjs*exp(-te3);
fd4=fd4+1/zgjs*exp(-te4);
fd5=fd5+1/zgjs*exp(-te5);
fd6=fd6+1/zgjs*exp(-te6);
fd7=fd7+1/zgjs*exp(-te7);
fd9=fd9+1/zgjs*exp(-te9);
fd10=fd10+1/zgjs*exp(-te10);



te8=log(MI*fp1/fp2);
ate8=ate8+1/zgjs*te8;
fd8=fd8+1/zgjs*exp(-te8);
%syms te2
%te2=te-te1;
%disp(double(te2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (zgjs*fi00<lk)&&(lk<zgjs*(fi00+fi10)+1)
        t0=0;T=zeros(1,2000);A=zeros(1,2000);A(1,1)=1;

R=rand(2,2000);
W0=[ah+ac,bh+bc;ah+ac,bh+bc];

W01=[f,dh+dc;f,dh+dc];
a=2;
for m=1:1:2000
    r1=R(1,m);r2=R(2,m);
    t1=vpasolve(-vpaintegral(W0(a,a),t3,t0,t2)-log(r1)==0,t2);
     t0=t1;
    if t0>tf
        m1=m+1;
        break
    end
    
    T(1,m+1)=t0;
    if a==1
        if r2*W0(a,a)<ah
            b=1;
        else
            b=2;
        end
    else
        if r2*W0(a,a)<bh
            b=1;
        else
            b=2;
        end
    end
    
    if a==1
        a=2;
    else 
        a=1;
    end
    A(1,m+1)=a-1;B(1,m)=b;
 end
 
  A1=zeros(1,m1);T1=zeros(1,m1);A1(1,m1)=a-1;T1(1,m1)=tf;B1=zeros(1,m1-2);
  for m=1:m1-1
      A1(1,m)=A(1,m);T1(1,m)=T(1,m);
  end
  for m=1:1:m1-2
      B1(1,m)=B(1,m);
  end

fp1=exp(-vpaintegral(W0(a,a),t3,T1(1,m1-1),tf));
fp3=exp(-vpaintegral(W01(a,a),t3,T1(1,m1-1),tf));

for m=1:m1-2
    fp1=fp1*R(1,m);
  fp3=fp3*exp(-vpaintegral(W01(A1(1,m)+1,A1(1,m)+1),t3,T1(1,m),T1(1,m+1)));
  
end
%fp1:y=0 forward probability fp2:y=0 backward probability fp3:y=1 forward probability
%fp8:y=1 backward probability %fp4:y=0 dual forward probability fp5:y=0 dual backward probability 
%fp6:y=1 dual forward probability fp7:y=1 dual backward probability
fp2=fp1;
fp4=fp1;
fp5=fp1;
fp8=fp3;
fp6=fp3;
fp7=fp3;




  BE1=zeros(1,m1-2); BE2=zeros(1,m1-2);BE31=zeros(1,m1-2);BE32=zeros(1,m1-2);
  BE41=zeros(1,m1-2);BE42=zeros(1,m1-2);BE51=zeros(1,m1-2);BE52=zeros(1,m1-2);
  BE8=zeros(1,m1-2);
 %BE61=zeros(1,m1-2);BE62=zeros(1,m1-2);BE71=zeros(1,m1-2);BE72=zeros(1,m1-2);
  
  for k=2:m1-1
      
          %W2:y=0 transition rate  W21:y=1 transition rate
         %W3: y=0 transition rate  W4:y=1 transition rate
      

         %BE1:y=0 forward transition probability (not include staying) BE2: y=0 backward transition probability (not include staying)
         %BE51:y=1 forward transition probability (not include staying) BE52: y=1 backward transition probability (not include staying)
         
         %BE31: y=0 dual forward transition probability    BE32: y=0 dual backward transition probability
         %BE41:y=1 dual forward transition probability     BE42: y=0 dual backward transition probability
         
         
         u2=[ah+ac,bh+bc;ah+ac,bh+bc];
         u21=[subs(f,t3,T1(1,k)),dh+dc;subs(f,t3,T1(1,k)),dh+dc];
         
         
          W32=[-u2(1,1),u2(1,2);1,1];
         W33=[0;fi00+fi10];
         W34=W32\W33;

         W42=[-u21(1,1),u21(1,2);1,1];
         W43=[0;fi01+fi11];
         W44=W42\W43;

         
          if B(1,k-1)==1
              W2=[ah,bh;ah,bh];
              W21=[subs(f1,t3,T1(1,k)),dh;subs(f1,t3,T1(1,k)),dh];

              W3=[W2(1,1),W2(2,1)*W34(1,1)/W34(2,1);W2(1,2)*W34(2,1)/W34(1,1),W2(2,2)];
              W4=[W21(1,1),W21(2,1)*W44(1,1)/W44(2,1);W21(1,2)*W44(2,1)/W44(1,1),W21(2,2)];
          else
              W2=[ac,bc;ac,bc];
              W21=[subs(f2,t3,T1(1,k)),dc;subs(f2,t3,T1(1,k)),dc];
             W3=[W2(1,1),W2(2,1)*W34(1,1)/W34(2,1);W2(1,2)*W34(2,1)/W34(1,1),W2(2,2)];
              W4=[W21(1,1),W21(2,1)*W44(1,1)/W44(2,1);W21(1,2)*W44(2,1)/W44(1,1),W21(2,2)];
          end


          
          BE1(1,k-1)=W2(A1(1,k)+1,A1(1,k-1)+1);
          BE2(1,k-1)=W2(A1(1,k-1)+1,A1(1,k)+1);
          BE51(1,k-1)=W21(A1(1,k)+1,A1(1,k-1)+1);
          BE52(1,k-1)=W21(A1(1,k-1)+1,A1(1,k)+1);
         

          BE31(1,k-1)=W3(A1(1,k)+1,A1(1,k-1)+1);
          BE32(1,k-1)=W3(A1(1,k-1)+1,A1(1,k)+1);
          BE41(1,k-1)=W4(A1(1,k)+1,A1(1,k-1)+1);
          BE42(1,k-1)=W4(A1(1,k-1)+1,A1(1,k)+1);
         
  end
fp1=fi10*fp1;fp3=fi11*fp3;fp4=fi10*fp4;fp6=fi11*fp6;
if A(1,m1-1)==0
    fp2=p00*fp2;fp5=fp5*p00;
    fp7=p01*fp7;fp8=fp8*p01;
    
    MI=p00*(fi10+fi11)/((p00+p01)*fi10);
else 
    fp2=fp2*p10;fp5=p10*fp5;
    fp7=p11*fp7;fp8=fp8*p11;
    
    MI=p10*(fi10+fi11)/((p10+p11)*fi10);

end
for m=1:m1-2
    fp1=fp1*BE1(1,m);
    fp2=fp2*BE2(1,m);
    fp3=fp3*BE51(1,m);
    fp8=fp8*BE52(1,m);
    
    fp4=fp4*BE31(1,m);
    fp5=fp5*BE32(1,m);
    fp6=fp6*BE41(1,m);
    fp7=fp7*BE42(1,m);
    
end

%ate1: toal           ate2:excess     ate3:housekeeping 
%ate4: marginal excess                ate5: marginal housekeeping  
%ate6: mixed excess                   ate7:mixed housekeeping  
te1=log(fp1/fp2);
te2=log(fp1/fp5);
te3=log(fp1/fp4);
te4=log((fp1+fp3)/(fp5+fp7));
te5=log((fp1+fp3)/(fp4+fp6));
te6=(log(fp1/fp5)-log((fp1+fp3)/(fp5+fp7)));
te7=(log(fp1/fp4)-log((fp1+fp3)/(fp4+fp6)));
te9=log((fp1+fp3)/(fp2+fp8));
te10=(log(fp1/fp2)-log((fp1+fp3)/(fp2+fp8)));


ate1=ate1+1/zgjs*log(fp1/fp2);
ate2=ate2+1/zgjs*log(fp1/fp5);
ate3=ate3+1/zgjs*log(fp1/fp4);
ate4=ate4+(1/zgjs)*log((fp1+fp3)/(fp5+fp7));
ate5=ate5+(1/zgjs)*log((fp1+fp3)/(fp4+fp6));
ate6=ate6+1/zgjs*(log(fp1/fp5)-log((fp1+fp3)/(fp5+fp7)));
ate7=ate7+1/zgjs*(log(fp1/fp4)-log((fp1+fp3)/(fp4+fp6)));
ate9=ate9+(1/zgjs)*log((fp1+fp3)/(fp2+fp8));
ate10=ate10+(1/zgjs)*(log(fp1/fp2)-log((fp1+fp3)/(fp2+fp8)));



fd1=fd1+1/zgjs*exp(-te1);
fd2=fd2+1/zgjs*exp(-te2);
fd3=fd3+1/zgjs*exp(-te3);
fd4=fd4+1/zgjs*exp(-te4);
fd5=fd5+1/zgjs*exp(-te5);
fd6=fd6+1/zgjs*exp(-te6);
fd7=fd7+1/zgjs*exp(-te7);
fd9=fd9+1/zgjs*exp(-te9);
fd10=fd10+1/zgjs*exp(-te10);

te8=log(MI*fp1/fp2);
ate8=ate8+1/zgjs*te8;
fd8=fd8+1/zgjs*exp(-te8);
%syms te2
%te2=te-te1;
%disp(double(te2));


%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
elseif ((fi00+fi10)*zgjs<lk)&&(lk<(fi00+fi10+fi01)*zgjs+1)
       
        t0=0;T=zeros(1,2000);A=zeros(1,2000);A(1,1)=0;

R=rand(2,2000);
W0=[f,dh+dc;f,dh+dc];
W01=[ah+ac,bh+bc;ah+ac,bh+bc];

a=1;
        
  for m=1:1:2000
    r1=R(1,m);r2=R(2,m);
    t1=vpasolve(-vpaintegral(W0(a,a),t3,t0,t2)-log(r1)==0,t2);
     t0=t1;
    if t0>tf
        m1=m+1;
        break
    end
    
    T(1,m+1)=t0;
    W1=[subs(f,t3,t0),dh+dc;subs(f,t3,t0),dh+dc];
    if a==1
        if r2*W1(a,a)<subs(f1,t3,t0)
            b=1;
        else
            b=2;
        end
    else
        if r2*W1(a,a)<dh
            b=1;
        else
            b=2;
        end
    end
    
    if a==1
        a=2;
    else 
        a=1;
    end
    A(1,m+1)=a-1;B(1,m)=b;
 end
 
  A1=zeros(1,m1);T1=zeros(1,m1);A1(1,m1)=a-1;T1(1,m1)=tf;B1=zeros(1,m1-2);
  for m=1:m1-1
      A1(1,m)=A(1,m);T1(1,m)=T(1,m);
  end
  for m=1:1:m1-2
      B1(1,m)=B(1,m);
  end

fp1=exp(-vpaintegral(W0(a,a),t3,T1(1,m1-1),tf));
fp3=exp(-vpaintegral(W01(a,a),t3,T1(1,m1-1),tf));

for m=1:m1-2
    fp1=fp1*R(1,m);
  fp3=fp3*exp(-vpaintegral(W01(A1(1,m)+1,A1(1,m)+1),t3,T1(1,m),T1(1,m+1)));
  
end
%fp1:y=1 forward probability fp2:y=1 backward probability fp3:y=0 forward probability
%fp8:y=0 backward probability %fp4:y=1 dual forward probability fp5:y=1 dual backward probability 
%fp6:y=0 dual forward probability fp7:y=0 dual backward probability
fp2=fp1;
fp4=fp1;
fp5=fp1;
fp8=fp3;
fp6=fp3;
fp7=fp3;




  BE1=zeros(1,m1-2); BE2=zeros(1,m1-2);BE31=zeros(1,m1-2);BE32=zeros(1,m1-2);
  BE41=zeros(1,m1-2);BE42=zeros(1,m1-2);BE51=zeros(1,m1-2);BE52=zeros(1,m1-2);
  BE8=zeros(1,m1-2);
 %BE61=zeros(1,m1-2);BE62=zeros(1,m1-2);BE71=zeros(1,m1-2);BE72=zeros(1,m1-2);
  
  for k=2:m1-1
      %W2:y=1 transition rate  W21:y=0 transition rate
         %W3: y=1 transition rate  W4:y=0 transition rate
         %BE1:y=1 forward transition probability (not include staying) BE2: y=1 backward transition probability (not include staying)
         %BE51:y=0 forward transition probability (not include staying) BE52: y=0 backward transition probability (not include staying)
         %BE31: y=1 dual forward transition probability    BE32: y=1 dual backward transition probability
         %BE41:y=0 dual forward transition probability     BE42: y=1 dual backward transition probability
      
         u2=[subs(f,t3,T1(1,k)),dh+dc;subs(f,t3,T1(1,k)),dh+dc];
         u21=[ah+ac,bh+bc;ah+ac,bh+bc];
     

          W32=[-u2(1,1),u2(1,2);1,1];
         W33=[0;fi01+fi11];
         W34=W32\W33;

         W42=[-u21(1,1),u21(1,2);1,1];
         W43=[0;fi00+fi10];
         W44=W42\W43;

        

         if B(1,k-1)==1
              W2=[subs(f1,t3,T1(1,k)),dh;subs(f1,t3,T1(1,k)),dh];
              W21=[ah,bh;ah,bh];
              W3=[W2(1,1),W2(2,1)*W34(1,1)/W34(2,1);W2(1,2)*W34(2,1)/W34(1,1),W2(2,2)];
              W4=[W21(1,1),W21(2,1)*W44(1,1)/W44(2,1);W21(1,2)*W44(2,1)/W44(1,1),W21(2,2)];
         else
             W2=[subs(f2,t3,T1(1,k)),dc;subs(f2,t3,T1(1,k)),dc];
             W21=[ac,bc;ac,bc];
              W3=[W2(1,1),W2(2,1)*W34(1,1)/W34(2,1);W2(1,2)*W34(2,1)/W34(1,1),W2(2,2)];
              W4=[W21(1,1),W21(2,1)*W44(1,1)/W44(2,1);W21(1,2)*W44(2,1)/W44(1,1),W21(2,2)];
          end

          
          BE1(1,k-1)=W2(A1(1,k)+1,A1(1,k-1)+1);
          BE2(1,k-1)=W2(A1(1,k-1)+1,A1(1,k)+1);
          BE51(1,k-1)=W21(A1(1,k)+1,A1(1,k-1)+1);
          BE52(1,k-1)=W21(A1(1,k-1)+1,A1(1,k)+1);
         

          BE31(1,k-1)=W3(A1(1,k)+1,A1(1,k-1)+1);
          BE32(1,k-1)=W3(A1(1,k-1)+1,A1(1,k)+1);
          BE41(1,k-1)=W4(A1(1,k)+1,A1(1,k-1)+1);
          BE42(1,k-1)=W4(A1(1,k-1)+1,A1(1,k)+1);
         
  end
 
fp1=fi01*fp1;fp3=fi00*fp3;fp4=fi01*fp4;fp6=fi00*fp6;
if A(1,m1-1)==0
    fp2=p01*fp2;fp5=fp5*p01;
    fp7=p00*fp7;fp8=fp8*p00;
    MI=p01*(fi00+fi01)/((p00+p01)*fi01);
else
    fp2=fp2*p11;fp5=p11*fp5;
    fp7=p10*fp7;fp8=fp8*p10;
    
    MI=p11*(fi00+fi01)/((p10+p11)*fi01);
end
for m=1:m1-2
    fp1=fp1*BE1(1,m);
    fp2=fp2*BE2(1,m);
    fp3=fp3*BE51(1,m);
    fp8=fp8*BE52(1,m);
    

    fp4=fp4*BE31(1,m);
    fp5=fp5*BE32(1,m);
    fp6=fp6*BE41(1,m);
    fp7=fp7*BE42(1,m);
    
end
%ate1: toal           ate2:excess     ate3:housekeeping 
%ate4: marginal excess                ate5: marginal housekeeping  
%ate6: mixed excess                   ate7:mixed housekeeping  
te1=log(fp1/fp2);
te2=log(fp1/fp5);
te3=log(fp1/fp4);
te4=log((fp1+fp3)/(fp5+fp7));
te5=log((fp1+fp3)/(fp4+fp6));
te6=(log(fp1/fp5)-log((fp1+fp3)/(fp5+fp7)));
te7=(log(fp1/fp4)-log((fp1+fp3)/(fp4+fp6)));
te9=log((fp1+fp3)/(fp2+fp8));
te10=(log(fp1/fp2)-log((fp1+fp3)/(fp2+fp8)));


ate1=ate1+1/zgjs*log(fp1/fp2);
ate2=ate2+1/zgjs*log(fp1/fp5);
ate3=ate3+1/zgjs*log(fp1/fp4);
ate4=ate4+(1/zgjs)*log((fp1+fp3)/(fp5+fp7));
ate5=ate5+(1/zgjs)*log((fp1+fp3)/(fp4+fp6));
ate6=ate6+1/zgjs*(log(fp1/fp5)-log((fp1+fp3)/(fp5+fp7)));
ate7=ate7+1/zgjs*(log(fp1/fp4)-log((fp1+fp3)/(fp4+fp6)));
ate9=ate9+(1/zgjs)*log((fp1+fp3)/(fp2+fp8));
ate10=ate10+(1/zgjs)*(log(fp1/fp2)-log((fp1+fp3)/(fp2+fp8)));



fd1=fd1+1/zgjs*exp(-te1);
fd2=fd2+1/zgjs*exp(-te2);
fd3=fd3+1/zgjs*exp(-te3);
fd4=fd4+1/zgjs*exp(-te4);
fd5=fd5+1/zgjs*exp(-te5);
fd6=fd6+1/zgjs*exp(-te6);
fd7=fd7+1/zgjs*exp(-te7);
fd9=fd9+1/zgjs*exp(-te9);
fd10=fd10+1/zgjs*exp(-te10);

te8=log(MI*fp1/fp2);
ate8=ate8+1/zgjs*te8;
fd8=fd8+1/zgjs*exp(-te8);
%syms te2
%te2=te-te1;
%disp(double(te2));
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
else
        t0=0;T=zeros(1,2000);A=zeros(1,2000);A(1,1)=1;

R=rand(2,2000);
W0=[f,dh+dc;f,dh+dc];
W01=[ah+ac,bh+bc;ah+ac,bh+bc];
        a=2;
for m=1:1:2000
    r1=R(1,m);r2=R(2,m);
    t1=vpasolve(-vpaintegral(W0(a,a),t3,t0,t2)-log(r1)==0,t2);
     t0=t1;
    if t0>tf
        m1=m+1;
        break
    end
    
    T(1,m+1)=t0;
    W1=[subs(f,t3,t0),dh+dc;subs(f,t3,t0),dh+dc];
    if a==1
        if r2*W1(a,a)<subs(f1,t3,t0)
            b=1;
        else
            b=2;
        end
    else
        if r2*W1(a,a)<dh
            b=1;
        else
            b=2;
        end
    end
    
    if a==1
        a=2;
    else 
        a=1;
    end
    A(1,m+1)=a-1;B(1,m)=b;
 end
 
  A1=zeros(1,m1);T1=zeros(1,m1);A1(1,m1)=a-1;T1(1,m1)=tf;B1=zeros(1,m1-2);
  for m=1:m1-1
      A1(1,m)=A(1,m);T1(1,m)=T(1,m);
  end
  for m=1:1:m1-2
      B1(1,m)=B(1,m);
  end

fp1=exp(-vpaintegral(W0(a,a),t3,T1(1,m1-1),tf));
fp3=exp(-vpaintegral(W01(a,a),t3,T1(1,m1-1),tf));

for m=1:m1-2
    fp1=fp1*R(1,m);
  fp3=fp3*exp(-vpaintegral(W01(A1(1,m)+1,A1(1,m)+1),t3,T1(1,m),T1(1,m+1)));
  
end
%fp1:y=1 forward probability fp2:y=1 backward probability fp3:y=0 forward probability
%fp8:y=0 backward probability %fp4:y=1 dual forward probability fp5:y=1 dual backward probability 
%fp6:y=0 dual forward probability fp7:y=0 dual backward probability
fp2=fp1;
fp4=fp1;
fp5=fp1;
fp8=fp3;
fp6=fp3;
fp7=fp3;




  BE1=zeros(1,m1-2); BE2=zeros(1,m1-2);BE31=zeros(1,m1-2);BE32=zeros(1,m1-2);
  BE41=zeros(1,m1-2);BE42=zeros(1,m1-2);BE51=zeros(1,m1-2);BE52=zeros(1,m1-2);
  BE8=zeros(1,m1-2);
 %BE61=zeros(1,m1-2);BE62=zeros(1,m1-2);BE71=zeros(1,m1-2);BE72=zeros(1,m1-2);
  
  for k=2:m1-1
      %W2:y=1 transition rate  W21:y=0 transition rate
         %W3: y=1 transition rate  W4:y=0 transition rate
         %BE1:y=1 forward transition probability (not include staying) BE2: y=1 backward transition probability (not include staying)
         %BE51:y=0 forward transition probability (not include staying) BE52: y=0 backward transition probability (not include staying)
         %BE31: y=1 dual forward transition probability    BE32: y=1 dual backward transition probability
         %BE41:y=0 dual forward transition probability     BE42: y=1 dual backward transition probability
      
         u2=[subs(f,t3,T1(1,k)),dh+dc;subs(f,t3,T1(1,k)),dh+dc];
         u21=[ah+ac,bh+bc;ah+ac,bh+bc];
     

          W32=[-u2(1,1),u2(1,2);1,1];
         W33=[0;fi01+fi11];
         W34=W32\W33;

         W42=[-u21(1,1),u21(1,2);1,1];
         W43=[0;fi00+fi10];
         W44=W42\W43;

        

         if B(1,k-1)==1
              W2=[subs(f1,t3,T1(1,k)),dh;subs(f1,t3,T1(1,k)),dh];
              W21=[ah,bh;ah,bh];
              W3=[W2(1,1),W2(2,1)*W34(1,1)/W34(2,1);W2(1,2)*W34(2,1)/W34(1,1),W2(2,2)];
              W4=[W21(1,1),W21(2,1)*W44(1,1)/W44(2,1);W21(1,2)*W44(2,1)/W44(1,1),W21(2,2)];
         else
             W2=[subs(f2,t3,T1(1,k)),dc;subs(f2,t3,T1(1,k)),dc];
             W21=[ac,bc;ac,bc];
              W3=[W2(1,1),W2(2,1)*W34(1,1)/W34(2,1);W2(1,2)*W34(2,1)/W34(1,1),W2(2,2)];
              W4=[W21(1,1),W21(2,1)*W44(1,1)/W44(2,1);W21(1,2)*W44(2,1)/W44(1,1),W21(2,2)];
          end


          
          BE1(1,k-1)=W2(A1(1,k)+1,A1(1,k-1)+1);
          BE2(1,k-1)=W2(A1(1,k-1)+1,A1(1,k)+1);
          BE51(1,k-1)=W21(A1(1,k)+1,A1(1,k-1)+1);
          BE52(1,k-1)=W21(A1(1,k-1)+1,A1(1,k)+1);
         

          BE31(1,k-1)=W3(A1(1,k)+1,A1(1,k-1)+1);
          BE32(1,k-1)=W3(A1(1,k-1)+1,A1(1,k)+1);
          BE41(1,k-1)=W4(A1(1,k)+1,A1(1,k-1)+1);
          BE42(1,k-1)=W4(A1(1,k-1)+1,A1(1,k)+1);
         
  end
fp1=fi11*fp1;fp3=fi10*fp3;fp4=fi11*fp4;fp6=fi10*fp6;
if A(1,m1-1)==0
    fp2=p01*fp2;fp5=fp5*p01;
    fp7=p00*fp7;fp8=fp8*p00;
    
    MI=p01*(fi10+fi11)/((p00+p01)*fi11);
else
    fp2=fp2*p11;fp5=p11*fp5;
    fp7=p10*fp7;fp8=fp8*p10;
    
    MI=p11*(fi10+fi11)/((p10+p11)*fi11);

end
for m=1:m1-2
    fp1=fp1*BE1(1,m);
    fp2=fp2*BE2(1,m);
    fp3=fp3*BE51(1,m);
    fp8=fp8*BE52(1,m);
   
    fp4=fp4*BE31(1,m);
    fp5=fp5*BE32(1,m);
    fp6=fp6*BE41(1,m);
    fp7=fp7*BE42(1,m);
    
end
%ate1: toal           ate2:excess     ate3:housekeeping 
%ate4: marginal excess                ate5: marginal housekeeping  
%ate6: mixed excess                   ate7:mixed housekeeping  
te1=log(fp1/fp2);
te2=log(fp1/fp5);
te3=log(fp1/fp4);
te4=log((fp1+fp3)/(fp5+fp7));
te5=log((fp1+fp3)/(fp4+fp6));
te6=(log(fp1/fp5)-log((fp1+fp3)/(fp5+fp7)));
te7=(log(fp1/fp4)-log((fp1+fp3)/(fp4+fp6)));
te9=log((fp1+fp3)/(fp2+fp8));
te10=(log(fp1/fp2)-log((fp1+fp3)/(fp2+fp8)));


ate1=ate1+1/zgjs*log(fp1/fp2);
ate2=ate2+1/zgjs*log(fp1/fp5);
ate3=ate3+1/zgjs*log(fp1/fp4);
ate4=ate4+(1/zgjs)*log((fp1+fp3)/(fp5+fp7));
ate5=ate5+(1/zgjs)*log((fp1+fp3)/(fp4+fp6));
ate6=ate6+1/zgjs*(log(fp1/fp5)-log((fp1+fp3)/(fp5+fp7)));
ate7=ate7+1/zgjs*(log(fp1/fp4)-log((fp1+fp3)/(fp4+fp6)));
ate9=ate9+(1/zgjs)*log((fp1+fp3)/(fp2+fp8));
ate10=ate10+(1/zgjs)*(log(fp1/fp2)-log((fp1+fp3)/(fp2+fp8)));



fd1=fd1+1/zgjs*exp(-te1);
fd2=fd2+1/zgjs*exp(-te2);
fd3=fd3+1/zgjs*exp(-te3);
fd4=fd4+1/zgjs*exp(-te4);
fd5=fd5+1/zgjs*exp(-te5);
fd6=fd6+1/zgjs*exp(-te6);
fd7=fd7+1/zgjs*exp(-te7);
fd9=fd9+1/zgjs*exp(-te9);
fd10=fd10+1/zgjs*exp(-te10);

te8=log(MI*fp1/fp2);
ate8=ate8+1/zgjs*te8;
fd8=fd8+1/zgjs*exp(-te8);
end
KK1(1,lk)=te1;KK2(1,lk)=te2;KK3(1,lk)=te3;KK4(1,lk)=te4;KK5(1,lk)=te5;KK6(1,lk)=te6;KK7(1,lk)=te7;KK8(1,lk)=te8;
end
ate1=double(ate1);
ate2=double(ate2);
ate3=double(ate3);
ate4=double(ate4);
ate5=double(ate5);
ate6=double(ate6);
ate7=double(ate7);
ate9=double(ate9);
ate10=double(ate10);
ate8=double(ate8);

fd1=double(fd1);
fd2=double(fd2);
fd3=double(fd3);
fd4=double(fd4);
fd5=double(fd5);
fd6=double(fd6);
fd7=double(fd7);
fd9=double(fd9);
fd10=double(fd10);
fd8=double(fd8);
ATE(abc,1)=ate1;ATE(abc,2)=ate2;ATE(abc,3)=ate3;ATE(abc,4)=ate4;ATE(abc,5)=ate5;ATE(abc,6)=ate6;ATE(abc,7)=ate7;ATE(abc,8)=ate8;ATE(abc,9)=ate9;ATE(abc,10)=ate10;
FD(abc,1)=fd1;FD(abc,2)=fd2;FD(abc,3)=fd3;FD(abc,4)=fd4;FD(abc,5)=fd5;FD(abc,6)=fd6;FD(abc,7)=fd7;FD(abc,8)=fd8;FD(abc,9)=fd9;FD(abc,10)=fd10;

end
mui=p00*log(p00/((p00+p01)*(p00+p10)))+p01*log(p01/((p00+p01)*(p01+p11)))+p10*log(p10/((p10+p11)*(p00+p10)))+p11*log(p11/((p10+p11)*(p01+p11)))-(fi00*log(fi00/((fi00+fi01)*(fi00+fi10)))+fi01*log(fi01/((fi00+fi01)*(fi01+fi11)))+fi10*log(fi10/((fi10+fi11)*(fi00+fi10)))+fi11*log(fi11/((fi10+fi11)*(fi01+fi11))));
ate1+mui
