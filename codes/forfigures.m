%% plot
p_list = [0,0.01,0.02,0.03,0.05,0.1,0.3,0.5,0.75,1.0];
for i = 1:length(p_list)
    p = p_list(i);
    [CT(1,i,:),CP(1,i,:),C(1,i,:),x(1,i,:)] = Task(p,1,0);
    [CT(2,i,:),CP(2,i,:),C(2,i,:),x(2,i,:)] = Task(p,2,0);
end

for i = 1:length(p_list)
    AAA(i,:) = zeros(1,100);
    AAA(i,1:100) = C(1,i,1:100);
    plot(1:100,AAA(i,1:100),'Linewidth',1.5);
    hold on;
end

legend('pi = 0','pi = 0.01','pi = 0.02','pi = 0.03','pi = 0.05','pi = 0.1','pi = 0.3','pi = 0.5','pi = 0.75','pi = 1','Location','northeast');

% plot(1:100,AAA,1:100,BBB,'LineWidth',1.5);
% xlabel('Time/yr');
% ylabel('Carbon sequestration in products/tCO2');
% legend('CP1','CP2','Location','northwest');
% grid on;
% 
% % plot(1:100,C1(1,1:100),1:100,C2(1,1:100),1:100,C3(1,1:100),1:100,C4(1,1:100),1:100,C5(1,1:100),1:100,C6(1,1:100),1:100,C7(1,1:100),1:100,C8(1,1:100),1:100,C9(1,1:100),'LineWidth',1.5);
% 
% xlabel('Time/yr');
% ylabel('Carbon sequestration/tCO2');
% title('Carbon sequestration at different deforestation rates');
% legend('pi = 0','pi = 0.01','pi = 0.02','pi = 0.03','pi = 0.05','pi = 0.1','pi = 0.3','pi = 0.5','pi = 1','Location','northwest');
% grid on;
% figure(2)
% plot(1:100,CP_5(1,1:100),1:100,CP_6(1,1:100));
% legend('0.5','1');



% AAA = zeros(1,100);BBB = zeros(1,100);
% AAA(1,1:100) = C1(2,2,1:100);BBB(1,1:100) = C(2,2,1:100);
% plot(1:100,C1(1:100),1:100,C2(1:100),'LineWidth',1.5);
% xlabel('Time/yr');
% ylabel('Carbon sequestration/tCO2');
% legend('Pi','Pi+10yr','Location','northwest');
% grid on;