clear;clc;
ws=zeros(10,1);
for i=1:10
    sim=simulation(25000, 0.5);
    ws(i)=sim.mean;
    
end 
