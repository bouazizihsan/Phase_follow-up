function x=Operator(x)
while(any(abs(x)>pi)) 
x=x.*(abs(x)<pi) + (x+2*pi).*(x<(-pi)) + (x-2*pi).*(x>pi); 
end
end