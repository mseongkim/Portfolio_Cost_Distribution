%unhedged_portfolio_cost_distribution
%sample code for unhedged portfolio only
%hedging instruments(e.g., trigger,stoploss,option) not included

function [histogram_xy, curve_xy] = ...
          unhedged_portfolio_distribution(Start_month, Start_year, End_month, End_year, n_obs, ...
          Forward_price, volume, corr_matrix, sigma1, sigma2, kappa, numpts )

%   Inputs
%       Start_month: Scalar. Start month of contract. Must be nonnegative. 
%       Start_year: Scalar. Start year of contract. Must be nonnegative.
%       End_month: Scalar. End month of contract. Must be nonnegative.
%       End_year: Scalar. End year of contract. Must be nonnegative. 
%       n_obs: Scalar. Number of Monte Carlo simulation. 
%              Suggested value = 10000. 
%       Forward_price: [# months] x [# zones]. Forward prices. 
%                       Must be nonnegative.
%       volume: [# months] x [# zones]. Monthly volatilities.
%                Must be nonnegative.
%       corr_matrix: [# zones] x [# zones]. Correlation matrix.
%                     Must be between -1 and 1 inclusive.
%       sigma1: [1] x [# zones]. Estimated multi-factor parameter(1)
%                Must be nonnegative.
%       sigma2: [1] x [# zones]. Estimated multi-factor parameter(2)
%                Must be nonnegative.
%       kappa: [1] x [# zones].  Estimated multi-factor parameter(3)
%       numpts: Scalar. Number of points (or bins) along histogram 
%               and normalized curve.
%
%   Outputs
%       histogram_xy: [numpts] x 2. Matrix of (x,y) coordinates of 
%                      histogram where x is center cost (of each bin) and y is probability.  
%       curve_xy: [numpts] x 2. Matrix of (x,y) coordinates of 
%                  normalized curve where x is cost and y is probability.

Forward_price_orig = Forward_price';
volume_orig = volume';
zos = size(Forward_price_orig,1); % number of products (zones)
mos = size(Forward_price_orig,2); % number of months

%error check 
if any( size(Forward_price) ~= size(volume) )
   error('Dimension mismatch: forward price is %dx%d.; volume is %dx%d.', ...
          mos, zos, size(volume_orig,2), size(volume_orig,1))
elseif ~isscalar(Start_month) || ~isscalar(Start_year) || ~isscalar(End_month) || ~isscalar(End_year)
   error('Month & Year must be scalar.')
elseif Start_month < 0 || Start_year < 0 || End_month < 0 || End_year < 0
   error('Month & Year must be nonnegative.')
elseif Start_year > End_year 
   error('Start year must be less than or equal to End year.')   
elseif n_obs < 1000
   error('Number of Monte Carlo simulations must be greater than 1000.')   
elseif any(any(Forward_price < 0)) 
   error('Forward price must be nonnegative.')     
elseif any(any(volume < 0)) 
   error('Volume must be nonnegative.')  
elseif any(any(corr_matrix < -1)) || any(any(corr_matrix > 1))
   error('Correlation values must be in range [-1,1].')  
elseif size(corr_matrix,1) ~= size(corr_matrix,2)
   error('Correlation matrix must be square %dx%d matrix.', ...
          zos, zos)
elseif any(any(sigma1 < 0)) || any(any(sigma2 < 0))
   error('Sigma 1 & 2 must be nonnegative.')        
elseif ~isscalar(numpts) || numpts < 2
   error('numpts must be scalar integer greater than 1')
end

% Check the positive definite of correlation matrix
corr_chol = chol(corr_matrix);
positivedefinite = all(eig(corr_matrix) > 0); 
eig(corr_matrix);

if positivedefinite == 0 
   error('Correlation matrix is not valid.'); 
end 

Today = datetime('today');
formatOut1 = 'mm/dd/yy';
Today_Str = datestr(Today,formatOut1);
End_day = eomday(End_year,End_month);

day_count_ttl = (datenum(End_year,End_month,End_day) - datenum(Today_Str));
day_count_start_end = (datenum(End_year,End_month,End_day)- datenum(Start_year,Start_month,1));
today_to_start = day_count_ttl - day_count_start_end -1; 
Time = (datenum(End_year,End_month,End_day) - datenum(Today_Str))./365;
Time_adj = [Time/day_count_ttl:Time/day_count_ttl:Time]';

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
Forward_price = reshape(repmat(Forward_price(num_zone,:),n_obs,1),[],n_obs*size(Time_numeric_ind,2));
local_vol_1 = reshape(repmat(local_vol_1,n_obs,1),[],n_obs*size(Time_numeric_ind,2));
local_vol_2 = reshape(repmat(local_vol_2,n_obs,1),[],n_obs*size(Time_numeric_ind,2));
time_to_maturity_ind_adj1 = Time_day_ind(1:end-1);
Price = zeros(time_to_maturity_ind_adj1(end)+1,n_obs*size(Time_numeric_ind,2));
Price(1,:) = Forward_price(1,:);
for i = 1:n_obs
    
    %Generate the correlated random number 1 and 2 using Cholesky
    %decomposition
    %Random number 1 and 2 are independent each other    
    Random_Sample_1_orig = (corr_chol'*normrnd(0,1,[4,size(local_vol_1,1)]))';
    Random_Sample_2_orig = (corr_chol'*normrnd(0,1,[4,size(local_vol_1,1)]))';
    Random_Sample_1 = Random_Sample_1_orig(:,num_zone); 
    Random_Sample_2 = Random_Sample_2_orig(:,num_zone);
    
    %Apply different expiration dates of monthly derivatives
    m = 1;
    x = [2:Time_day_ind(m)+1];
    while m < size(Time_day_ind,2)
          for j = x(1):x(end) 
                Price(j,i:n_obs:end) = Price(j-1,i:n_obs:end).*exp(-0.5.*(local_vol_1(j-1,i:n_obs:end)+local_vol_2(j-1,i:n_obs:end))+...
                                       sqrt(local_vol_1(j-1,i:n_obs:end)).*Random_Sample_1(j)+ ...
                                       sqrt(local_vol_2(j-1,i:n_obs:end)).*Random_Sample_2(j));
          end
          x = [Time_day_ind(m)+2:Time_day_ind(m+1)+1];
          m = m + 1;
    end             
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

%unhedged portfolio 
figure(1)
h_1 = histogram(unhedged_strategy,numpts,'Normalization','pdf');
[cnt_un, xcost_un] = hist(unhedged_strategy,numpts);

%normalize
yprob_un = cnt_un./(sum(cnt_un)*mean(diff(xcost_un)));
histogram_xy = [xcost_un', yprob_un'];
%also return kernel smoothing function estimate
%min and max are extended up and down by 10% of range
mincost_un = min(unhedged_strategy) - 0.1*range(unhedged_strategy);
maxcost_un = max(unhedged_strategy) + 0.1*range(unhedged_strategy);
xfcost_un = mincost_un:(maxcost_un-mincost_un)/(numpts-1):maxcost_un;
[yfprob_un, xfcost_un] = ksdensity(unhedged_strategy,xfcost_un);
lognparms_un = lognfit(unhedged_strategy);
yfprob_un = lognpdf(xfcost_un,lognparms_un(1),lognparms_un(2));

yfprob_un = yfprob_un ./ (sum(yfprob_un)*(maxcost_un-mincost_un)/(numpts-1));
curve_xy = [xfcost_un',yfprob_un'];
hold on;
plot(xfcost_un,yfprob_un,'r','LineWidth',1) 
title('Unhedged Portfolio Cost Distribution')
xlabel('Cost')
ylabel('Probability')

toc

end



