clear; 
clc;
cpu_i=cputime;
load('M.mat');
%M = 1/20*ones(20,20);
steps = 4;
[Bn, Bn] =size(M); 

CDF = M; %cumulative distribution matrix
for j = 1:Bn
    for k = 2:Bn
        CDF(j,k) = CDF(j,k)+CDF(j,k-1);
    end
end
CDF(:,Bn) = ones(Bn,1);

%% load in times sort
load('times.mat');
t1 = [t1; t2-t1; t3-t2; t4-t3]; %case3 CTRW/BERN; case2; case1
%t1 = [t1]; %case3 SMM;
N =10^6;
%%%New
t1 = datasample(t1,N,'Weights',t1); %按t1的权重扩展原t1，增加至N个
t1 = sort(t1);
t1 = reshape(t1,[],Bn);
%tbins = linspace(0.1,5,100); %case1
tbins = linspace(0.01,3,100); %case2
%tbins = linspace(0.01,2,100); %case3
%%%%%%%%

[aa, bb] =size(t1);
ind = ones(aa,bb);
for j = 1:Bn
   ind(:,j) = j*ones(aa,1); 
end

%% Correlated SMM
t = reshape(t1,[],1);
vtd = reshape(ind,[],1);
N =length(t); 
vtd2 = zeros(N,1);

for j =2:steps
    j
   rand_i = randi(10000,N,1)/10000; %由介于 1 和 10000 之间的伪随机整数组成的 N×1 数组
   rand_i(find(rand_i==1)) = .99999;  
    
   for i = 1:length(vtd)
    a = vtd(i); %gives bin that particle was in after first cell. 
                %This tells us the row to choose in the cumulative probaility matrix
    vtd2(i) = find(CDF(a,:)>rand_i(i),1); %vtd2 is vector of bins that particle will be after cell2
    %we assign bins to vtd2 by finding first time our random integer is greater than a cell in cdm matrix
   end
   
   rr=randi([1,aa],N,1); %the row
   rr = sub2ind(size(t1),rr,vtd2); %sub2ind函数将矩阵中指定元素的行列下标转换成存储的序号
   
   ts = t1(rr);
   t= t+ts;
   vtd  =vtd2; 
   
%    if j==1
%      BTC(:,1) = histc(t,tbins);   %histc函数制定数值边界为分界条件    
%    end
%    if j==steps*(1/4)
%      BTC(:,1) = histc(t,tbins);       
%    end
   if j==steps*(2/4)
     BTC(:,2) = histc(t,tbins);       
   end
   if j==steps*(3/4)
     BTC(:,3) = histc(t,tbins);       
   end  
    if j==steps
     BTC(:,4) = histc(t,tbins);       
    end  
end

%% CTRW and Bernoulli model
t = reshape(t1,[],1);
t_ctrw = t(randi([1,length(t)],N,1)); %t1 is the vector of \tau1
t_bern = t_ctrw;
t_temp = t_bern; 
p =0.8;

for j=2:steps
       
        t_ctrw = t_ctrw + t(randi([1,length(t)],N,1));
        
        flip = rand(N,1); 
        ind = flip>p; %the flips where the time changes  
    
        t_temp(ind) = t(randi([1,length(t)],sum(ind),1));
        t_bern= t_bern + t_temp; 
          
        if j==steps*(2/4)
          BTC_CTRW(:,2) = histc(t_ctrw,tbins);   
          BTC_BERN(:,2) = histc(t_bern,tbins);
        end
        if j==steps*(3/4)
          BTC_CTRW(:,3) = histc(t_ctrw,tbins);   
          BTC_BERN(:,3) = histc(t_bern,tbins);       
        end  
        if j==steps
          BTC_CTRW(:,4) = histc(t_ctrw,tbins);   
          BTC_BERN(:,4) = histc(t_bern,tbins);       
        end  
end


%% Measure_data

  load('times_measured2.mat'); %不同的case导入相应的文件

  %BTC_M(:,1) = histc(t11,tbins);    
  BTC_M(:,2) = histc(t22,tbins);
  BTC_M(:,3) = histc(t33,tbins);
  BTC_M(:,4) = histc(t44,tbins);
  
 %% error analysis
 
 bern = log(BTC_BERN);
 ctrw = log(BTC_CTRW);
 smm = log(BTC);
 measure = log(BTC_M); 
 
 for  i = 2:steps
     
      BERN = find(bern(:,i)>0);
      CTRW = find(ctrw(:,i)>0);
      SMM = find(smm(:,i)>0);
      MEASURE = find(measure(:,i)>0);

      %找到三个模型、实验值横坐标的交集
      trans1 = intersect(BERN,CTRW);
      trans2 = intersect(trans1,SMM);
      tag = intersect(trans2,MEASURE);
      
      if i==2 
          tag2=tag;
      end
      if i==3
         tag3=tag;
      end
      if i==4
         tag4=tag;
      end
 end

%用0填充三个模型、实验值的横坐标个数至所有steps的最大值
N=max([length(tag2),length(tag3),length(tag4)]);
Tag=[padarray(tag2,[N-length(tag2) 0],'post') padarray(tag3,[N-length(tag3) 0],'post') padarray(tag4,[N-length(tag4) 0],'post')];

%mape = zeros(3,3); %mean absolute percent error

 for i = 1:steps-1    
   vec= Tag(find(Tag(:,i))+(i-1)*N); %找到对应的模型和实验结果的相同横坐标值（matlab先按照行计数，再按列）
   mape(1,i)=mean(abs((smm(vec,i+1)-measure(vec,i+1))./smm(vec,i+1)));  %SMM & Measure
   %mape(2,i)=mean(abs((bern(vec,i+1)-measure(vec,i+1))./bern(vec,i+1))); %BERN & Measure
   %mape(3,i)=mean(abs((ctrw(vec,i+1)-measure(vec,i+1))./ctrw(vec,i+1)));  %CTRW & Measure
 end

%  stem((1:3),mape(1,:))
%  hold on
%  stem((1:3),mape(2,:))
%  hold on
%  stem((1:3),mape(3,:))

%%
% dt = zeros(length(tbins),1);
% for i =2:length(dt)
%    dt(i) = tbins(i) - tbins(i-1);    
% end
% dt(1) = dt(2);

%% plot CTRW & Measured data
loglog(tbins, BTC_CTRW(:,2), 'k',tbins,BTC_M(:,2),'*k','Linewidth',1)
hold on
loglog(tbins, BTC_CTRW(:,3), '-.k',tbins,BTC_M(:,3),'+k','Linewidth',1)
hold on
loglog(tbins, BTC_CTRW(:,4), '--k',tbins,BTC_M(:,4),'xk','Linewidth',1)
hold on
set(gca,'FontSize',10,'Fontname','Times New Roman');
xlabel('t (s)','FontSize',10,'Fontname','Times New Roman');
ylabel('N','FontSize',10,'Fontname','Times New Roman');
xlim([0.67,2.5]); %case3(0,2) case1(1,3.7) case2(0.67,2.5)
ylim([10^0,10^6]);
annotation('textbox',[0.730867346938775,0.839816325047738,0.189795924223547,0.067346940258386],...
           'LineStyle','none','string','u=0.2 m/s','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
legend('CTRW-step 2', 'measured data-step 2','CTRW-step 3', 'measured data-step 3','CTRW-step 4', 'measured data-step 4',...
       'FontSize',10,'Fontname','Times New Roman','EdgeColor','w');
   
%% plot Bernoulli & Measured data
loglog(tbins, BTC_BERN(:,2), 'k',tbins,BTC_M(:,2),'*k','Linewidth',1)
hold on
loglog(tbins, BTC_BERN(:,3), '-.k',tbins,BTC_M(:,3),'+k','Linewidth',1)
hold on
loglog(tbins, BTC_BERN(:,4), '--k',tbins,BTC_M(:,4),'xk','Linewidth',1)
hold on
set(gca,'FontSize',10,'Fontname','Times New Roman');
xlabel('t (s)','FontSize',10,'Fontname','Times New Roman');
ylabel('N','FontSize',10,'Fontname','Times New Roman');
xlim([0.67,2.5]); %case3(0,2) case1(1,3.7) case2(0.67,2.5)
ylim([10^0,10^6]);
annotation('textbox',[0.730867346938775,0.839816325047738,0.189795924223547,0.067346940258386],...
           'LineStyle','none','string','u=0.2 m/s','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
legend('Bernoulli model-step 2', 'measured data-step 2','Bernoulli model-step 3', 'measured data-step 3','Bernoulli model-step 4', 'measured data-step 4',...
       'FontSize',10,'Fontname','Times New Roman','EdgeColor','w');

%% plot SMM & Measured data
loglog(tbins, BTC(:,2), 'k',tbins,BTC_M(:,2),'*k','Linewidth',1)
hold on
loglog(tbins, BTC(:,3), '-.k',tbins,BTC_M(:,3),'+k','Linewidth',1)
hold on
loglog(tbins, BTC(:,4), '--k',tbins,BTC_M(:,4),'xk','Linewidth',1)
hold on
set(gca,'FontSize',10,'Fontname','Times New Roman');
xlabel('t (s)','FontSize',10,'Fontname','Times New Roman');
ylabel('N','FontSize',10,'Fontname','Times New Roman');
xlim([0,2]); %case3(0,2) case1(1,3.7) case2(0.67,2.5)
ylim([10^0,10^6]);
annotation('textbox',[0.730867346938775,0.839816325047738,0.189795924223547,0.067346940258386],...
           'LineStyle','none','string','u=0.25 m/s','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
legend('SMM-step 2', 'measured data-step 2','SMM-step 3', 'measured data-step 3','SMM-step 4', 'measured data-step 4',...
       'FontSize',10,'Fontname','Times New Roman','EdgeColor','w');

%% calculate cpu time   
cpu = cputime-cpu_i; %seconds

