SI_consts;
cmap = custom_cmap;
%name = 'redo_posi/';
%name = 'ramp_it_up/';
name = 'why/';

data_root = '/Users/sgess/Desktop/sims/data/qp_tars/2015/';
plot_root = '/Users/sgess/Desktop/plots/QP/2015/Fall/';
data_dir = [data_root name];
plot_dir = [plot_root name];

if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end

% Plot units, natural or real
plot_units = 'real';

% Save data?
save_plot = 10;
save_ext  = '.png';

% Get plasma and beam parameters
%load([data_dir 'param.mat']);
n0 = 8e16;
[omega_p, lambda_p, skin_depth, plasma_time, plasma_period, E0, beta_p] = plasma_parameters(n0);
input_struct.plasma.SD = skin_depth;
input_struct.plasma.density = n0;
input_struct.plasma.field = E0;

sz = 50E-6;
zz_cos = linspace(0,20*sz,1001);
zz_rho = linspace(-10*sz,10*sz,1001);
zz_pad = linspace(-10*sz,30*sz,2001);
dz = zz_rho(2)-zz_rho(1);

N = 0.5E10;

rho_b = N/(sqrt(2*pi)*sz)*exp(-zz_rho.^2/(2*sz^2));
rho_pad = N/(sqrt(2*pi)*sz)*exp(-zz_pad.^2/(2*sz^2));

a = 240/skin_depth;
b = 294/skin_depth;


QEB_type = 'QEB-XZ';
FEZX_type = 'FEZ-XZ';
QEP1_type = 'QEP1-XZ';

files = dir([data_root name FEZX_type]);
nFiles = sum(strncmp('FEZ',{files.name},3));
mean_e = zeros(1,nFiles);
mean_loss = zeros(1,nFiles);
for file_number = 10:10:(nFiles*10)
    ind = file_number/10;
    
    
    field_EZX = LOAD_DATA(data_dir, FEZX_type, file_number);
    [x_axis, z_axis] = LOAD_AXIS(data_dir, FEZX_type, file_number);
    [iter, dt, time] = LOAD_TIME(data_dir, FEZX_type, file_number);
    s = time*skin_depth/1e4;
    
    beam_rho = LOAD_DATA(data_dir, QEB_type, file_number);
    plas_rho = LOAD_DATA(data_dir, QEP1_type, file_number);
    beam = LOAD_BEAM(data_dir,file_number);
    
    mean_e(ind) = mean(SI_eM*1e-3*beam(:,6));
    mean_loss(ind) = 1000*(20.35 - mean_e(ind));
    
    n = -n0*plas_rho(70,1);
    [~, ~, skin_depth2] = plasma_parameters(n);
    xx = linspace(skin_depth*x_axis(1),skin_depth*x_axis(2),size(plas_rho,1));
    zz = linspace(skin_depth*z_axis(1),skin_depth*z_axis(2),size(plas_rho,2));
    zz_ana = 1e6*zz_pad-13;
    [~,c1] = min(abs(zz_ana-zz(1)));
    [~,c2] = min(abs(zz_ana-zz(end)));
    
    rho = -n0*(20*beam_rho+plas_rho);
    EzF = 1000*E0*field_EZX(256,:);
    
    a = 230/skin_depth2;
    b = 290/skin_depth2;
    
    [Ez,Ez0,w0] = holo_wake(n,a,b,zz_cos);
    Wz = -dz*conv(rho_b,Ez,'full');
    
    figure(4);
    
    clf;
    imagesc(zz,xx,rho); colormap(cmap.bwr); axis xy;
    xlabel('Z [\mum]');
    ylabel('Y [\mum]');
    set(gca,'fontsize',16);
    caxis([-8e16 8e16]);
    box off;
    hold on;
    plot(zz,EzF,'k',zz_ana(c1:c2),Wz(c1:c2)/1e6,'k--','linewidth',2);
    a2 = axes('YAxisLocation', 'Right');
    set(a2, 'color', 'none');
    set(a2, 'XTick', []);
    set(a2, 'YLim', [xx(1) xx(end)]);
    box on;
    ylabel('eE_z [MeV/m]');
    title(['Energy Loss = ' num2str(mean_loss(ind),'%0.2f') ' MeV, s = ' num2str(s,'%0.2f') ' cm']);
    set(gca,'fontsize',16);
    %line([zz(1) zz(end)],[xx(end) xx(end)],'k','linewidth',2);
    set(gcf,'color','w');
    hold off;
    %pause;
    F = getframe(gcf);
    figure(5);
    imshow(F.cdata);
    num_str = num2str(file_number,'%04d');
    if save_plot
        saveas(gcf,[plot_dir 'charge_rho2_' num_str save_ext]);
    else
        pause;
    end
end