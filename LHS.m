function L= LHS (x_low,x_high,e_low,e_high,n)
x_r=x_high-x_low;
e_r=e_high-e_low;
x_lows=zeros(n,size(x_high,2));
x_highs=zeros(n,size(x_high,2));
for i=1:size(x_high,2)
x_lows(:,i)=x_low(:,i):x_r(:,i)/n:x_high(:,i)-(x_r(:,i)/n);
x_highs(:,i)=x_low(:,i)+x_r(:,i)/n:x_r(:,i)/n:x_high(:,i);
end

e_lows=zeros(n,size(e_high,2));
e_highs=zeros(n,size(e_high,2));
for i=1:size(e_high,2)
e_lows(:,i)=e_low(:,i):e_r(:,i)/n:e_high(:,i)-(e_r(:,i)/n);
e_highs(:,i)=e_low(:,i)+e_r(:,i)/n:e_r(:,i)/n:e_high(:,i);
end

x=unifrnd(x_lows,x_highs);
e=unifrnd(e_lows,e_highs);
xx=x;
ee=e;

for i=1:size(x_high,2)
perm=randperm(n);
x(:,i)=xx(perm);
end
for i=1:size(e_high,2)
perm=randperm(n);
e(:,i)=ee(perm);
end

perm=randperm(n);
perm_e=randperm(n);
L.x=x(perm,:);
L.e=e(perm_e,:);
end
