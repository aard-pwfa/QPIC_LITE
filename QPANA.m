% QuickPIC ANALYSIS SCRIPT
% SPENCER M.F. GESSNER 9/8

clear all;

data_dir = '/Users/sgess/Desktop/FACET/qp_tars/2012/Sep/07/QP_eShort2/';
plot_dir = '/Users/sgess/Desktop/FACET/QuickPIC/PLOTS/2012/Sep/07/QP_eShort2/';

QEB_type = 'QEB-XZ';
QEP_type = 'QEP1-XZ';
FEX_type = 'FEX-XZ';
FEY_type = 'FEY-XZ';
FEZ_type = 'FEZ-XZ';

file_number = 1;

% Load beam charge density
beam_rho = LOAD_DATA(data_dir, QEB_type, file_number);
plas_rho = LOAD_DATA(data_dir, QEP_type, file_number);
field_EX = LOAD_DATA(data_dir, FEX_type, file_number);
field_EY = LOAD_DATA(data_dir, FEY_type, file_number);
field_EZ = LOAD_DATA(data_dir, FEZ_type, file_number);

[x_axis, z_axis] = LOAD_AXIS(data_dir, QEB_type, file_number);

XX = linspace(x_axis(1),x_axis(2),length(beam_rho));
ZZ = linspace(z_axis(1),z_axis(2),length(beam_rho));

% figure;
% %imagesc(ZZ,XX,abs(log(plas_rho)));
% imagesc(ZZ(256)-ZZ,XX(129:256),abs(log(plas_rho(129:256,:))));
% xlabel('c/\omega_p','fontsize',16);
% ylabel('c/\omega_p','fontsize',16);
% colorbar;
% t = colorbar('peer',gca);
% set(get(t,'ylabel'),'String', 'n_0','fontsize',16);
% title('Log of Plasma Density','fontsize',16);
% %saveas(gca,[plot_dir 'log_plas_rho.pdf']);
% 
% figure;
% %imagesc(ZZ,XX,abs(log(beam_rho)));
% imagesc(ZZ(256)-ZZ,XX(129:256),abs(log(beam_rho(129:256,:))));
% xlabel('c/\omega_p','fontsize',16);
% ylabel('c/\omega_p','fontsize',16);
% colorbar;
% t = colorbar('peer',gca);
% set(get(t,'ylabel'),'String', 'n_0','fontsize',16);
% title('Log of Beam Density','fontsize',16);
% %saveas(gca,[plot_dir 'log_beam_rho.pdf']);

% figure;
% %imagesc(ZZ,XX,plas_rho);
% imagesc(ZZ(256)-ZZ,XX(129:256),plas_rho(129:256,:));
% xlabel('c/\omega_p','fontsize',16);
% ylabel('c/\omega_p','fontsize',16);
% colorbar;
% t = colorbar('peer',gca);
% set(get(t,'ylabel'),'String', 'n_0','fontsize',16);
% title('Plasma Density','fontsize',16);
% v =  axis;
% text(3*v(2)/5,5*v(4)/6,'Beam Direction \rightarrow','FontSize',16,'FontWeight','bold');
% saveas(gca,[plot_dir 'plas_rho.pdf']);
% 
% figure;
% %imagesc(ZZ,XX,beam_rho);
% imagesc(ZZ(256)-ZZ,XX(129:256),beam_rho(129:256,:));
% xlabel('c/\omega_p','fontsize',16);
% ylabel('c/\omega_p','fontsize',16);
% colorbar;
% t = colorbar('peer',gca);
% set(get(t,'ylabel'),'String', 'n_0','fontsize',16);
% title('Beam Density','fontsize',16);
% v =  axis;
% text(3*v(2)/5,5*v(4)/6,'Beam Direction \rightarrow','FontSize',16,'FontWeight','bold');
% saveas(gca,[plot_dir 'beam_rho.pdf']);
% 
% 
% figure;
% %imagesc(ZZ,XX,field_EZ);
% imagesc(ZZ(256)-ZZ,XX(129:256),field_EZ(129:256,:));
% xlabel('c/\omega_p','fontsize',16);
% ylabel('c/\omega_p','fontsize',16);
% colorbar;
% t = colorbar('peer',gca);
% set(get(t,'ylabel'),'String', 'm c \omega_p / e','fontsize',16);
% title('Longitudinal E Field','fontsize',16);
% v =  axis;
% text(3*v(2)/5,5*v(4)/6,'Beam Direction \rightarrow','FontSize',16,'FontWeight','bold');
% saveas(gca,[plot_dir 'EZ_2D.pdf']);
% 
% 
% figure;
% plot(ZZ(256)-ZZ,field_EZ(129,:));
% xlabel('c/\omega_p','fontsize',16);
% ylabel('m c \omega_p / e','fontsize',16);
% title('Longitudinal E Field on Axis','fontsize',16);
% axis([0 12.5566 (min(field_EZ(129,:))-.01) (max(field_EZ(129,:))+.01)]);
% v =  axis;
% text(2*v(2)/3,4*v(3)/5,'Beam Direction \rightarrow','FontSize',16,'FontWeight','bold');
% saveas(gca,[plot_dir 'EZ_1D.pdf']);
save('../../../COMPARE/QP_eShort2.mat','ZZ','field_EZ');