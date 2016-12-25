%unhedged_sample_v1
%sample code for unhedged portfolio only
%hedging instruments(e.g., trigger,stoploss,option) not included
%assumed that market data and parameters are given from SQL server

close all; clear all;
tic

%%User input
Start_month = 1;
Start_year = 2018;
End_month = 12;
End_year = 2018;
n_obs = 10000; % Number of Monte Carlo simulation

volume = ...
[70	100	80	30
 90	110	90	20
 90	120	110	80
 50	110	120	70
 40	120	100	80
 80	150	120	90
 50	160	110	70
 70	140	100	80
 80	130	90	100
 90	120	110	100
 80	180	120	90
 70	210	130	100];

% User should also choose the products to draw data from SQL server accordingly. 

%%Database(SQL server) input
sigma1 = [0.6 0.4 0.7 0.3]; % Estimated multi-factor parameter (1) from observed market data
sigma2 = [0.2 0.15 0.3 0.15]; % Estimated multi-factor parameter(2)from observed market data
kappa = [3 2.5 3.5 2.7]; % Estimated multi-factor parameter(3)from observed market data

Forward_price = ... % 1-year monthly market forward prices for four products
[40	41.5 42	38
 41	42.5 43	39
 42	43.5 44	40
 43	44.5 45	41
 44	45.5 46	42
 45	46.5 47	43
 46	47.5 48	44
 47	48.5 49	45
 48	49.5 50	46
 49	50.5 51	47
 50	51.5 52	48
 51	52.5 53	49];

corr_matrix = [1 0.6 0.4 0.8;0.6 1 0.2 0.7;0.4 0.2 1 0.1;0.8 0.7 0.1 1];

%%%% End of inputs %%%%

Today = datetime('today');
Forward_price_orig = Forward_price';
volume_orig = volume';
zos = size(Forward_price_orig,1); % number of products
mos = size(Forward_price_orig,2); % number of months

corr_chol = chol(corr_matrix);
positivedefinite = all(eig(corr_matrix) > 0); % Check the PSD of corr_matrix
eig(corr_matrix);

if positivedefinite == 0 
   msgbox('Correlation Matrix is not valid'); 
end 

formatOut1 = 'mm/dd/yy';
Today_Str = datestr(Today,formatOut1);
End_day = eomday(End_year,End_month);

day_count_ttl = (datenum(End_year,End_month,End_day) - datenum(Today_Str));
day_count_start_end = (datenum(End_year,End_month,End_day)- datenum(Start_year,Start_month,1));
today_to_start = day_count_ttl - day_count_start_end -1; 
Time = (datenum(End_year,End_month,End_day) - datenum(Today_Str))./365;
Time_adj = [Time/day_count_ttl:Time/day_count_ttl:Time]'; %T0

num_month = floor((datenum(End_year,End_month,End_day) - datenum(Start_year,Start_month,1))/30);
Month = Start_month:1:Start_month+num_month-1;
Month_adj = mod(Month,12);
Month_adj(Month_adj==0) = 12;
[wid len] = size(Month_adj);

Year = Start_year:End_year;
Year_wid = size(Year,2);
num_12 = find(Month_adj == 12);
ind = [num_12,len];

%Creat Year mapping table
Year_ind = [1:ind(1)];
Year_mapping = zeros(1,len);
for i = 1:size(num_12,2)
     Year_mapping(Year_ind) = Year(i);
     Year_ind = [ind(i)+1:ind(i+1)];
end 
Year_mapping(Year_mapping == 0) = Year(end);

%Calculate numeric time index
if Start_year == End_year
   End_day_month = eomday(End_year,Month_adj);
   for i = 1:12
        Time_numeric_ind(i) = (datenum(End_year,Month_adj(i),End_day_month(i))-datenum(Today_Str))./365;
   end
   End_day_month = [today_to_start End_day_month];
   Time_day_ind = cumsum(End_day_month);
else 
    End_day_month = zeros(1,len);
    Time_numeric_ind = zeros(1,len);
    for i = 1:len
        End_day_month(i) = eomday(Year_mapping(i),Month_adj(i));
        Time_numeric_ind(i) = (datenum(Year_mapping(i),Month_adj(i),End_day_month(i))-datenum(Today_Str))./365;
    end
        End_day_month = [today_to_start End_day_month];
        End_day_month = End_day_month(End_day_month ~=0);
        Time_day_ind = cumsum(End_day_month);
        Time_numeric_ind = Time_numeric_ind(Time_numeric_ind ~=0);
end
[wid2 len2] = size(Time_numeric_ind);

for num_zone = 1:zos
Forward_price = Forward_price_orig;
volume = volume_orig;

%Calculate two factors 
m=2;
factor1 = zeros(Time_day_ind(end),len2);
factor2 = zeros(Time_day_ind(end),len2);
for i = 1:len2
    for j = 1:Time_day_ind(m)
        factor1(j,i) =(sigma1(:,num_zone)^2)/(2*kappa(:,num_zone))*(exp(-2*kappa(:,num_zone)*(Time_numeric_ind(i)...
                       - Time_adj(j)))-exp(-2*kappa(:,num_zone)*Time_numeric_ind(i)));
        factor2(j,i) = sigma2(:,num_zone)^2*Time_adj(j);
        ttl_variance(j,i) = factor1(j,i) + factor2(j,i);
    end
        m = m + 1;
end

%Calculate local varianes
m=2;
local_vol_1 = zeros(Time_day_ind(end),len2);
local_vol_2 = zeros(Time_day_ind(end),len2);
for i = 1:len2
    for j = 1:Time_day_ind(m)-1
        local_vol_1(j+1,i) = factor1(j+1,i) - factor1(j,i);
        local_vol_2(j+1,i) = factor2(j+1,i) - factor2(j,i);
        local_vol_1(1,i) = factor1(1,i);   
        local_vol_2(1,i) = factor2(1,i);
    end
        m = m + 1;
end

%Conduct Monte Carlo simulations
Forward_price = reshape(repmat(Forward_price(num_zone,:),n_obs,1),[],n_obs*len2);
local_vol_1 = reshape(repmat(local_vol_1,n_obs,1),[],n_obs*len2);
local_vol_2 = reshape(repmat(local_vol_2,n_obs,1),[],n_obs*len2);
Time_day_ind_adj1 = Time_day_ind(1:end-1);
Time_day_ind_adj1 = reshape(repmat(Time_day_ind_adj1,n_obs,1),[],n_obs*len2);
Price = zeros(Time_day_ind_adj1(end)+1,n_obs*len2);
Price(1,:) = Forward_price(1,:);
m = 1;
for i = 1:len2*n_obs
    %Generate correlated random numbers
    Random_Sample_1_orig = (corr_chol'*normrnd(0,1,[zos,size(local_vol_1,1)]))';
    Random_Sample_2_orig = (corr_chol'*normrnd(0,1,[zos,size(local_vol_1,1)]))';
    Random_Sample_1 = Random_Sample_1_orig(:,num_zone); 
    Random_Sample_2 = Random_Sample_2_orig(:,num_zone);
    for j = 2:Time_day_ind_adj1(m)+1
            %Apply two factors stochastic model with two random numbers 
            Price(j,i) = Price(j-1,i)*exp(-0.5*(local_vol_1(j-1,i)+local_vol_2(j-1,i))+...
                         sqrt(local_vol_1(j-1,i))*Random_Sample_1(j-1)+ ...
                         sqrt(local_vol_2(j-1,i))*Random_Sample_2(j-1));           
    end
            m = m + 1;
end
Price(1,:) = []; 

Cal_Price = zeros(size(Price,1),n_obs);
P = zeros(size(Price,1),mos);
for i = 1:size(Price,1) 
    for j = 1:n_obs
            P(i,:) = Price(i,j:n_obs:end);
            Cal_Price(i,j) = sum(Price(i,j:n_obs:end))/numel(P(P(i,:)~=0));
    end
end

m = 1;
n = 1;
x = [1:Time_day_ind(n)];
for i = 1:n_obs
    while m < (size(Time_day_ind,2)-1)*n_obs + 1 
          avg_price(n,i) = sum(Price(x(1):x(end),m))/(x(end)-x(1)+1);
          x = [Time_day_ind(n)+1:Time_day_ind(n+1)];
          m = m + n_obs;
          n = n + 1;
    end
          m = i + 1;
          n = 1;
          x = [1:Time_day_ind(n)];
end

for a = 1:n_obs
    unhedged_strategy_drft(num_zone,a) = mean(avg_price(:,a));
end

end

if num_zone == 1
   unhedged_strategy = unhedged_strategy_drft;
else   
   unhedged_strategy = sum(unhedged_strategy_drft)/num_zone;
end

numpts = 40;

%unhedged portfolio 
figure(1)
h_1 = histogram(unhedged_strategy,numpts,'Normalization','pdf');
[cnt_un, xcost_un] = hist(unhedged_strategy,numpts);
 
%normalize
yprob_un = cnt_un./(sum(cnt_un)*mean(diff(xcost_un)));
pdfdist_un = [xcost_un' yprob_un'*10^6];
%also return kernel smoothing function estimate
%min and max are extended up and down by 10% of range
mincost_un = min(unhedged_strategy) - 0.1*range(unhedged_strategy);
maxcost_un = max(unhedged_strategy) + 0.1*range(unhedged_strategy);
xfcost_un = mincost_un:(maxcost_un-mincost_un)/(numpts-1):maxcost_un;
[yfprob_un, xfcost_un] = ksdensity(unhedged_strategy,xfcost_un);
lognparms_un = lognfit(unhedged_strategy);
yfprob_un = lognpdf(xfcost_un,lognparms_un(1),lognparms_un(2));
yfprob_un = yfprob_un ./ (sum(yfprob_un)*(maxcost_un-mincost_un)/(numpts-1));
pdfdist_un = [pdfdist_un xfcost_un' yfprob_un'*10^6];
hold on;
plot(xfcost_un,yfprob_un,'r','LineWidth',1) 
title('Unhedged Portfolio Cost Distribution')
xlabel('Cost')
ylabel('Probability')

toc


