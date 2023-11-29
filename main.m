clear all
clc
% close all

set(0, 'DefaultLineLineWidth', 2); %set thickness of all the lines = 2

n_pigment=1.196;  %real refractive index of pigment
k_pigment=0;  %imaginary refractive index of pigment
n_medium = 1;
nang=20000;
no_of_f_v=50;
no_of_xs=200;

f_v=logspace(-4,log10(0.75),no_of_f_v);
xs = logspace(-2,3,no_of_xs)';
teta=linspace(eps,pi,nang)';%don't start with zero to avoid division by zero
u=cos(teta);


gridmap=zeros(no_of_f_v*no_of_xs,2);
counter=1;
for i=1:length(xs)
    for j=1:length(f_v)
        gridmap(counter,1)=xs(i);
        gridmap(counter,2)=f_v(j);
        counter=counter+1;
    end
end
figure
s=scatter(gridmap(:,2),gridmap(:,1),3,'filled');
hAx=gca;
s.MarkerFaceColor = [1 0 0];
set(gca,'XScale', 'log')
set(gca,'YScale', 'log')
ylabel('Size parameter, \chi')
xlabel('Volume fraction, f_v')
hAx.XColor = [0 0 0];
hAx.YColor = [0 0 0];
hAx.LineWidth = 1.5;
axis square
xlim([10^-4 0.75])
ylim([10^-2 10^3])
hAx.XAxis.TickValues = [10^-4 10^-3 10^-2 10^-1 0.75];
xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','0.75'})
set(gca,'FontSize',13)
box on
saveas(gcf,['grid.png'])
% saveas(gcf,['grid.emf'])
saveas(gcf,['grid.fig'])


f_v_cr_1 = @(x_s) (0.9047./(pi./(2*x_s)+1)).^3;%drolen tien line

drolen_fv=f_v_cr_1(xs);
drolen_fv=drolen_fv(drolen_fv>0.006);
drolen_xs=xs(end-length(drolen_fv)+1:end);
drolen_fv=[0.006;drolen_fv];
drolen_xs=[10^-2;drolen_xs];

mishchenko_fv=[0.02,0.05,0.1,0.02,0.05,0.1]; %mishchenko's experiments
mishchenko_x=[4.208,4.208,4.208,4.916,4.916,4.916];

Qsca_ind = zeros(length(xs),1);
g_ind = zeros(length(xs),1);
Qsca_dep = zeros(length(xs),length(f_v));
g_dep = zeros(length(xs),length(f_v));

tic
parfor i=1:length(xs)
    x_0=xs(i)/n_medium; %size parameter
    [Qsca_ind(i),g_ind(i),S1,S2] = mie(x_0,n_medium,n_pigment,k_pigment,nang);
    Qsca_dep_i= zeros(1,length(f_v));
    g_dep_i= zeros(1,length(f_v));
    for j=1:length(f_v)
        S11=0.5*(abs(S1).^2+abs(S2).^2)';
        S=SSF_correction(f_v(j), teta, xs(i));
        Qsca_dep_i(j)=2*trapz(teta,sin(teta) .* S11 .* S)/abs(xs(i))^2;
        g_dep_i(j)=trapz(teta,sin(teta) .*cos(teta) .* S11 .* S)/trapz(teta,sin(teta) .* S11 .* S);
    end
    Qsca_dep(i,:)=Qsca_dep_i;
    g_dep(i,:)=g_dep_i;
end
toc

Qsca_dep_transport=Qsca_dep.*(1-g_dep);
Qsca_ind_transport=Qsca_ind.*(1-g_ind);

tra_error_monodisperse=100*abs(Qsca_ind_transport - Qsca_dep_transport)./Qsca_dep_transport;
error_monodisperse=100*abs(Qsca_ind - Qsca_dep)./Qsca_dep;
g_err_monodisperse=(g_ind - g_dep);
% g_err_monodisperse=100*abs(g_ssf-g_ssf(1,:))./g_ssf;

figure('Renderer', 'painters', 'Position', [500 300 528 420]) % starting point and height - width of the frame
contourf(f_v,xs,tra_error_monodisperse,'edgecolor','none','LevelList',[0:0.1:5])
hAx=gca;
caxis([0 5])
a = colorbar;
a.Label.String = 'Error = (\sigma_s''_i_n_d - \sigma_s''_d_e_p) / \sigma_s''_d_e_p  [%]';
a.Label.FontSize = 13;
set(gca,'XScale', 'log')
set(gca,'YScale', 'log')
ylabel('Size parameter')
xlabel('Volume fraction, f_v')
colormap(jet)
hold on
% [C,h]=contour(f_v,xs,error_monodisperse',[1 2 3 5 10 20],"ShowText",true,'edgecolor','black');
[C,h]=contour(f_v,xs,tra_error_monodisperse,[1 3 5],"ShowText",true,'edgecolor','black');
set(gca,'XScale', 'log')
set(gca,'YScale', 'log')
ylabel('Size parameter, \chi')
xlabel('Volume fraction, f_v')
hAx.XColor = [0 0 0];
hAx.YColor = [0 0 0];
hAx.LineWidth = 1.5;
axis square
set(h,'linewidth',1.5)
set(gca,'FontSize',13)
clabel(C,h,'FontSize',13,'Color','Black')
hold on
loglog(drolen_fv,drolen_xs,'--k')
scatter(mishchenko_fv([1,4]),mishchenko_x([1,4]),60,'MarkerEdgeColor',[0 0.8 0],...
          'LineWidth',1.5)
scatter(mishchenko_fv([2,5]),mishchenko_x([2,5]),60,'MarkerEdgeColor',[0.75, 0.75, 0],...
          'LineWidth',1.5)
scatter(mishchenko_fv([3,6]),mishchenko_x([3,6]),60,'MarkerEdgeColor',[1 0 0],...
          'LineWidth',1.5)
xlim([10^-4 0.75])
ylim([10^-2 10^3])
hAx.XAxis.TickValues = [10^-4 10^-3 10^-2 10^-1 0.75];
xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','0.75'})
txt = ['m = ' num2str(n_pigment/n_medium) ' + ' num2str(k_pigment/n_medium) '⋅i'];
t=text(1.25*10^-4,5*10^2,txt);
t.FontSize = 13;
t.Color = [1 1 1];

saveas(gcf,['transport_scat_coeff_ref.png'])
% saveas(gcf,['transport_scat_coeff_ref.emf'])
saveas(gcf,['transport_scat_coeff_ref.fig'])


figure('Renderer', 'painters', 'Position', [500 300 528 420]) % starting point and height - width of the frame

contourf(f_v,xs,error_monodisperse,'edgecolor','none','LevelList',[0:0.1:5])
hAx=gca;
caxis([0 5])
a = colorbar;
a.Label.String = 'Error = (\sigma_s_,_i_n_d - \sigma_s_,_d_e_p) / \sigma_s_,_d_e_p  [%]';
a.Label.FontSize = 13;
set(gca,'XScale', 'log')
set(gca,'YScale', 'log')
ylabel('Size parameter')
xlabel('Volume fraction, f_v')
colormap(jet)
hold on
% [C,h]=contour(f_v,xs,error_monodisperse',[1 2 3 5 10 20],"ShowText",true,'edgecolor','black');
[C,h]=contour(f_v,xs,error_monodisperse,[1 3 5],"ShowText",true,'edgecolor','black');
set(gca,'XScale', 'log')
set(gca,'YScale', 'log')
ylabel('Size parameter, \chi')
xlabel('Volume fraction, f_v')
hAx.XColor = [0 0 0];
hAx.YColor = [0 0 0];
hAx.LineWidth = 1.5;
axis square
set(h,'linewidth',1.5)
set(gca,'FontSize',13)
clabel(C,h,'FontSize',13,'Color','Black')
hold on
loglog(drolen_fv,drolen_xs,'--k')
scatter(mishchenko_fv([1,4]),mishchenko_x([1,4]),60,'MarkerEdgeColor',[0 0.8 0],...
          'LineWidth',1.5)
scatter(mishchenko_fv([2,5]),mishchenko_x([2,5]),60,'MarkerEdgeColor',[0.75, 0.75, 0],...
          'LineWidth',1.5)
scatter(mishchenko_fv([3,6]),mishchenko_x([3,6]),60,'MarkerEdgeColor',[1 0 0],...
          'LineWidth',1.5)
xlim([10^-4 0.75])
ylim([10^-2 10^3])
hAx.XAxis.TickValues = [10^-4 10^-3 10^-2 10^-1 0.75];
xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','0.75'})
txt = ['m = ' num2str(n_pigment/n_medium) ' + ' num2str(k_pigment/n_medium) '⋅i'];
t=text(1.25*10^-4,5*10^2,txt);
t.FontSize = 13;
t.Color = [1 1 1];
saveas(gcf,['scat_coeff_ref.png'])
% saveas(gcf,['scat_coeff_ref.emf'])
saveas(gcf,['scat_coeff_ref.fig'])

% g_err_monodisperse(g_err_monodisperse<0)=0;
figure('Renderer', 'painters', 'Position', [500 300 528 420]) % starting point and height - width of the frame

contourf(f_v,xs,g_err_monodisperse,'edgecolor','none','LevelList',[0:0.0005:0.05])
hAx=gca;
caxis([0 0.05])
a = colorbar;
a.Label.String = 'Error = g_i_n_d - g_d_e_p';
a.Label.FontSize = 13;
set(gca,'XScale', 'log')
set(gca,'YScale', 'log')
ylabel('Size parameter')
xlabel('Volume fraction, f_v')
colormap(jet)
hold on
[C,h]=contour(f_v,xs,g_err_monodisperse,[10 10],"ShowText",true,'edgecolor','black');
set(gca,'XScale', 'log')
set(gca,'YScale', 'log')
ylabel('Size parameter, \chi')
xlabel('Volume fraction, f_v')
hAx.XColor = [0 0 0];
hAx.YColor = [0 0 0];
hAx.LineWidth = 1.5;
axis square
set(h,'linewidth',1.5)
set(gca,'FontSize',13)
clabel(C,h,'FontSize',13,'Color','Black')
hold on
loglog(drolen_fv,drolen_xs,'--k')
scatter(mishchenko_fv([1,4]),mishchenko_x([1,4]),60,'MarkerEdgeColor',[0 0.8 0],...
          'LineWidth',1.5)
scatter(mishchenko_fv([2,5]),mishchenko_x([2,5]),60,'MarkerEdgeColor',[0.75, 0.75, 0],...
          'LineWidth',1.5)
scatter(mishchenko_fv([3,6]),mishchenko_x([3,6]),60,'MarkerEdgeColor',[1 0 0],...
          'LineWidth',1.5)
xlim([10^-4 0.75])
ylim([10^-2 10^3])
hAx.XAxis.TickValues = [10^-4 10^-3 10^-2 10^-1 0.75];
xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','0.75'})
txt = ['m = ' num2str(n_pigment/n_medium) ' + ' num2str(k_pigment/n_medium) '⋅i'];
t=text(1.25*10^-4,5*10^2,txt);
t.FontSize = 13;
t.Color = [1 1 1];
saveas(gcf,['asymm.png'])
% saveas(gcf,['asymm.emf'])
saveas(gcf,['asymm.fig'])

