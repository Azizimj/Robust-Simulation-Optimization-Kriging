function kr = kriging (e_pred)
%% kriging
global n
global x_mat
global w_mean
global e_mat
global x_pred

% syms h1;
% 
% % c(h)= .001*(1-1.5*h/.5+.5*(h/.5).^3);
%     function c=c(h)
%         if h==0
%         c(h)=100;
%         elseif (h>0) && (h<=0.5) 
%         c(h)=100-(100*(1.5*h/.5-.5*(h/.5).^3));
%         else
%         c(h)=0;
%         end
%     end
% c1(h1)= 100*(1-(1.5*h1/.5-.5*(h1/.5).^3));
% c1(h1)= 100*(1-exp(-1*(3*h1)^3/(0.5^2)));
ex_points=horzcat( x_mat,e_mat); %n*(d+d_1)
hh=ex_points-horzcat(repmat(x_pred,n,1),repmat(e_pred,n,1));
% syms h;
% for i=1:n
%     h(i,1)=sqrt( sum( hh(i,:).^2 ) );
% end
h=sum(hh.^2,2).^(1/size(hh,2));
% k=c1(  h  ); %h is n*1 and thus k must be n*1 
k=c(  h  ); %h is n*1 and thus k must be n*1 

K=zeros(n);
for i=1:n
    for j=1:n
    hh= ex_points(i,:)-ex_points(j,:);
%     K(i,j)=c1(sqrt(sumsqr(hh)));
    K(i,j)=c(sumsqr(hh).^(1/size(hh,2)));
    end
end
 
landa = inv(K) *( k + ones(n,1)* (1-ones(1,n)*inv(K)*k ) / (ones(1,n)*inv(K)*ones(n,1)) ) ; 
 
%kr.landa  = landa; % kriging landa n*1
kr   =-1* w_mean'*landa; % krigign prediction 1*1
 
end
