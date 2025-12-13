% Funcao y=rect(t)

function y=rect(t)

a=(sign(t+0.5)-sign(t-0.5))/2;
%l=length(t);
%for cont=1:l
%   if a(cont,1)<=0.5
%      a(cont,1)=0;
%   end;
%   if a(cont,1)>0.5
%      a(cont,1)=1;
%   end;
%   
%end;

y=a;      
%y=round(a);      
