function Current = getCurrent(beam,N,n0,beam_ax)

SI_consts;
[~, ~, skin_depth] = plasma_parameters(n0);

data = hist(skin_depth*beam(:,3),beam_ax);

d_data_ax = 1e-6*(beam_ax(2)-beam_ax(1));

Current = SI_c*SI_e*data*N/(1000*sum(data)*d_data_ax);