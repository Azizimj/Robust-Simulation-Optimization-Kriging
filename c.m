function c=c(h)
rang=80; sill=100;
c=zeros(size(h,1),1);
for i=1:size(h,1)        
    if h(i)==0
        c(i)=sill;
    elseif (h(i)>0) && (h(i)<=rang)
        c(i)=sill*(1-(1.5*h(i)/rang-.5*(h(i)/rang).^3));
    else
        c(i)=0;
    end
%     c(i)=(1-heaviside(h(i)-rang))*((1-sign(h(i)))*sill + heaviside(rang-h(i))*sign(h(i))*(sill-(sill*(1.5*h(i)/rang-.5*(h(i)/rang).^3))));
% c(i)=sill-(sill*(1.5*h(i)/rang-.5*(h(i)/rang).^3));
% end
end