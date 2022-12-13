%% plot
p_list = [0,0.01,0.02,0.03,0.05,0.1,0.3,0.5,0.75,1.0];
for i = 1:length(p_list)
    p = p_list(i);
    [CT(1,i,:),CP(1,i,:),C(1,i,:),x(1,i,:)] = Task(p,1,0);
    [CT(2,i,:),CP(2,i,:),C(2,i,:),x(2,i,:)] = Task(p,2,0);
end

norx = zeros(2,length(p_list),3);
for i = 1:length(p_list)
    for j = 1:2
        for k = 1:3
            if k == 3
                norx(j,i,k) = (max(min(x(:,:,k)))-x(j,i,k)) / (max(max(x(:,:,k)))-min(min(x(:,:,k))))+0.001;
            else
                norx(j,i,k) = (x(j,i,k)-min(min(x(:,:,k)))) / (max(max(x(:,:,k)))-min(min(x(:,:,k))))+0.001;
            end
        end
    end
end

Y = zeros(2,length(p_list),3);
for i = 1:length(p_list)
    for j = 1:2
        for k = 1:3
            Y(j,i,k) = norx(j,i,k)/sum(sum(sum(norx)));
        end
    end
end

E = zeros(1,3);
for k = 1:3
    for j = 1:2
        for i = 1:length(p_list)
            E(1,k) = E(1,k) - 1/18*Y(j,i,k)*log(Y(j,i,k));
        end
    end
end

W = zeros(1,3);
for i = 1:3
    W(1,i) = (1-E(1,i))/(3-sum(E));
end

S = zeros(2,length(p_list));
for j = 1:2
    for i = 1:length(p_list)
        for k = 1:3
            S(j,i) = S(j,i) + W(1,k)*Y(j,i,k);
        end
    end
end
% 
% norx = zeros(2,length(p_list),2);
% for i = 1:length(p_list)
%     for j = 1:2
%         for k = 1:2
%             norx(j,i,k) = (x(j,i,k)-min(min(x(:,:,k)))) / (max(max(x(:,:,k)))-min(min(x(:,:,k))))+0.001;
%         end
%     end
% end
% 
% Y = zeros(2,length(p_list),2);
% for i = 1:length(p_list)
%     for j = 1:2
%         for k = 1:2
%             Y(j,i,k) = norx(j,i,k)/sum(sum(sum(norx)));
%         end
%     end
% end
% 
% E = zeros(1,2);
% for k = 1:2
%     for j = 1:2
%         for i = 1:length(p_list)
%             E(1,k) = E(1,k) - 1/18*Y(j,i,k)*log(Y(j,i,k));
%         end
%     end
% end
% 
% W = zeros(1,2);
% for i = 1:2
%     W(1,i) = (1-E(1,i))/(2-sum(E));
% end
% 
% S = zeros(2,length(p_list));
% for j = 1:2
%     for i = 1:length(p_list)
%         for k = 1:2
%             S(j,i) = S(j,i) + W(1,k)*Y(j,i,k);
%         end
%     end
% end


plot(1:length(p_list),S(1,1:length(p_list)),1:length(p_list),S(2,1:length(p_list)),'LineWidth',1.5);