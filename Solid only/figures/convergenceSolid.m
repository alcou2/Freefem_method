
width = 4; height=2;

% Create figure
figure1 = figure;
figure1.PaperUnits = 'inches';
figure1.PaperPosition = [0 0 2 1];

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

% Create multiple lines using matrix input to plot
plot1 = plot(mhS,EigFreq_bar_P1(1,:),mhS,EigFreq_bar_P2(1,:));
set(plot1(1),'DisplayName','P1 elements','Marker','square','Color',[0 0 0]);
set(plot1(2),'DisplayName','P2 elements','Marker','o','Color',[0 0 1]);

% Create asymptotic value lines
f1 = 3.55*sqrt(D/rho_s);
f2 = 8.11*sqrt(D/rho_s);
plot2 = plot([0 max(mhS)],[f1 f1],[0 max(mhS)],[f2 f2]);

% Create ylabel
ylabel('$$ \bar{\omega} $$');

% Create xlabel
xlabel('Number of element in the height of the plate');

% Set the remaining axes properties
set(axes1,'XTick',[1:1:5],'YTick',[40:20:220]);
% Create legend
legend(axes1,'show');

print -depsc convergenceSolid.eps;
