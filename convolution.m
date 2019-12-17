close all;
clear all;
clc;
x=input('enter the value of 1st sequence');
h=input('enter the value of 2nd sequence');
disp('the 1st sequence is-');
disp(x); 
disp('the 2nd sequence is-');
disp(h);
lx=length(x);
lh=length(h);
n=max(lx,lh);
subplot(3,1,1);
stem(x);
title('1st sequence');
subplot(3,1,2);
stem(h);
title('2nd sequence');
hh=[h zeros(1,n-lh)];
xx=zeros(n);
xx(1:lx,1)=x;
for i=2:n
    for j=2:n
        xx(j,i)=xx(j-1,i-1);
      
    end;
end;
for b=1:n-1
    xx(1,b+1)=xx(n,b);
end;
yy=xx*hh';
subplot(3,1,3);
stem(yy);
disp('convoluted o/p->');
disp(yy');
title('y=x*h');