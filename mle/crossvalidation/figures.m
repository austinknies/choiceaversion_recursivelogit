% Creates cross-validation figures
load('LogL_RLCA.mat')
load('LogL_RL.mat')
load('LogL_RL_LS.mat')


figure
r = plot(LogL_RL);
hold on
s = plot(LL_CA);
hold on
t = plot(LogL_LS);
hold off
r.LineWidth = 2;
s.LineWidth = 2;
t.LineWidth = 2;
xlim ([1 20])
legend({'RL', 'RL-CA','RL-LS'},'Location','east')
xlabel('Training Sample')
ylabel('Estimated LL')


figure
r = plot(PLogL_RL);
hold on
s = plot(PLL_CA);
hold on
t = plot(PLogL_LS);
hold off
r.LineWidth = 2;
s.LineWidth = 2;
t.LineWidth = 2;
xlim ([1 20])
legend({'RL', 'RL-CA','RL-LS'},'Location','southeast')
xlabel('Testing Sample')
ylabel('Predicted LL')

% errors
err_RL = -PLogL_RL/366; % following Mai, Fosgerau, and Frejinger (2017)
err_CA = -PLL_CA/366;
err_LS = -PLogL_LS/366;

figure
r = plot(err_RL);
hold on
s = plot(err_CA);
hold on
t = plot(err_LS);
hold off
r.LineWidth = 2;
s.LineWidth = 2;
t.LineWidth = 2;
xlim ([1 20])
legend({'RL', 'RL-CA','RL-LS'},'Location','southeast')
xlabel('Testing Sample')
ylabel('Test Errors')