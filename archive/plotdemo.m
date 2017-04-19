% This is a demo of plot multi graph in a figure

x = 0:0.1:10;
y1 = sin(2*x);
y2 = cos(2*x);

figure
subplot(2,2,1)
plot(x, y1)
title('Subplot 1')

subplot(2,2,2)
scatter(x, y2)
title('Subplot 2')

subplot(2,2,[3 4])
yyaxis left
plot(x, y1)
yyaxis right
plot(x, y2)
title('Subplot 3')

print('pdemo', '-bestfit', '-dpdf')

%Following is a block comment
%{
load fisheriris
pred = meas(51:end, 1:2);
resp = (1:100)'>50;
mdl = fitglm(pred, resp, 'Distribution', 'binomial', 'Link', 'logit');
scores = mdl.Fitted.Probability;
[X, Y, T, AUC] = perfcurve(species(51:end, :), scores, 'virginica');
AUC
%}
