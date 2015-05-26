%n = 5:20;
n = 20;
k = 1:20;
pos = [1 1];
num_sim = 20;
glob_ct_arr = zeros(1,num_sim);
best_ct_arr = zeros(1,num_sim);
avg_ct_arr = zeros(1,num_sim);
worst_ct_arr = zeros(1,num_sim);
%average value for each value of n
glob_ct = zeros(1,length(n));
best_ct = zeros(1,length(n));
avg_ct = zeros(1,length(n));
worst_ct = zeros(1,length(n));

gossiping_time = zeros(1,length(n));

for i = 1:length(k)
    disp(i);
    
    for j =1:num_sim
        [ glob_ct_arr(j), best_ct_arr(j),avg_ct_arr(j), worst_ct_arr(j) ] = RandomWalk( k(i), n, pos );
    end
    glob_ct(i) = mean(glob_ct_arr);
    best_ct(i) = mean(best_ct_arr);
    avg_ct(i) = mean(avg_ct_arr);
    worst_ct(i) = mean(worst_ct_arr);
    
    if (i==1)
        gossiping_time(i) = 0;
    else
        gossiping_time(i) = (n*n) * (log(n)*log(n)) / i;
    end
end
% disp(sprintf('global cover time %d',glob_ct));
% disp(sprintf('best local cover time %d',best_ct));
% disp(sprintf('avarage local cover time %d',avg_ct));
% disp(sprintf('worst local cover time %d',worst_ct));

%theoretical
% k increasing
theo_ct = zeros(1,length(k));
max_bound_theo_ct = zeros(1,length(k));
max_bound_theo_ct2 = zeros(1,length(k));

for i =1:length(k)
    theo_ct(i) = (2/k(i))*(n^2)*2*log(n);
    max_bound_theo_ct(i) = theo_ct(i) + gossiping_time(i);
    
    max_bound_theo_ct2(i) = glob_ct(1)/i + gossiping_time(i);
end

figure(1);
grid on;
hold on; 
plot(k,glob_ct,'k+-');
plot(k,best_ct,'b*-');
plot(k,avg_ct,'rs-');
plot(k,worst_ct,'m^-');
plot(k,theo_ct,'bd--','LineWidth',2);
plot(k,max_bound_theo_ct,'LineWidth',2);
plot(k,max_bound_theo_ct2,'LineWidth',2);
h = legend('glob','best','avg','worst','theo','max theo');
set(h,'Location','NorthEast');
xlabel('Number of random walkers');
ylabel('Cover time');


% n increasing
% theo_ct = zeros(1,length(n));
% for i =1:length(n)
%     theo_ct(i) = (n(i)^2)*2*log(n(i));
% end
% clf('reset');
% figure('name', 'k=2, n increasing');
% hold on; 
% plot(n,glob_ct,'k+-');
% plot(n,best_ct,'b*-');
% plot(n,avg_ct,'rs-');
% plot(n,worst_ct,'m^-');
% plot(n,theo_ct,'bd--','LineWidth',2);
% h = legend('glob','best','avg','worst','theo');
% set(h,'Location','NorthWest');
% xlabel('Length of side of torus');
% ylabel('Cover time');
