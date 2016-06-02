% set up indices of events to use

clear; close all; clc;

load chiouYoungsUsed % indices of ground motions used in Chiou and Youngs (2014)
load NGA_W2_corr_meta_data EQID NGA_num

% initialize
compositeUsed = zeros(21539,1);
compositeUsed(chiouYoungsUsed) = 1; % include all records that Chiou and Youngs used

save compositeUsed compositeUsed


%% look at what is being excluded

figure
plot(NGA_num, EQID, 'or')
hold on
% plot(NGA_num(chiouYoungsUsed), EQID(chiouYoungsUsed), '.b')
plot(NGA_num(find(compositeUsed)), EQID(find(compositeUsed)), '.b')

% plot(NGA_num(exclIdx), EQID(exclIdx), '^g')
plot([8163 8163], [0 1500], ':k')
xlabel('NGA number')
ylabel('EQID')
legend('All data', 'Used by Chiou and Youngs')
axis([0 23000 0 1500])

