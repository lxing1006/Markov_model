%% Transition Matrix
%Need t1, t2
clc;clear;
load('times.mat');

%Assign bins after first cell
Bn =11; %number of bins
% N=10000;
% t1 = datasample(t1,N,'Weights',t1); %按t1的权重扩展原t1，增加至N个
% t1 = sort(t1);
% t2 = datasample(t2,N,'Weights',t2); %按t1的权重扩展原t1，增加至N个
% t2 = sort(t2);
v1 = t1;
v2 = t2-t1;
%tbins depends on the max and min of v2
%tbins = linspace(.597,.83,Bn); %case-1
%tbins = linspace(.385,.509,Bn); %case-3
tbins = linspace(0.3560,0.5761,Bn); %case2
%tbins = linspace(min(v2),max(v2),Bn);

tbins(1)=0;
tbins(end)=100;
[~, ind1] = histc(v1,tbins);
[~,ind2] = histc(v2,tbins); 

M = zeros(Bn-1,Bn-1);
for i = 1:length(ind1)
   start = ind1(i);
   stop = ind2(i); 
   M(start,stop) = M(start,stop)+1;
end
for i=1:Bn-1
   M(i,:) = M(i,:)/sum(M(i,:));
end

save('M.mat','M');
figure(3)
title('Transition Matrix')
imagesc(M')
colorbar

CC=1/(Bn-1)*sum(sum(M.^2)) %correlation coefficient: calculate how important correlation is 
