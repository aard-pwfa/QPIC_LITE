function fit_data = BEAM_ANA(beam,dim,proj,nbin,dlim,param_struct)

DataDim = beam(:,dim);

if proj ~= 0
    if dim < 4
        DataDim = param_struct.plasma.SD*beam(:,dim)+proj*beam(:,dim+3)./beam(:,6);
    end
else
    if dim < 4
        DataDim = param_struct.plasma.SD*beam(:,dim);
    end
end

if isempty(dlim)
    d_edge = linspace(min(DataDim),max(DataDim),nbin);
else
    d_edge = linspace(dlim(1),dlim(2),nbin);
end

[prof,axis] = hist(DataDim,d_edge);

c = wm(axis,prof,1);
rms = wm(axis,prof,2);

off_guess = mean(prof(1)+prof(end))/2;
amp_guess = max(prof)-off_guess;
sig_guess = rms;
asy_guess = 0;
cen_guess = c;
[x, fval] = fminsearch(@(x) fun2min(x,prof,axis),[amp_guess,cen_guess,asy_guess,sig_guess,off_guess]);

fit_data.amp = x(1);
fit_data.cen = x(2);
fit_data.asy = x(3);
fit_data.sig = x(4);
fit_data.off = x(5);
fit_data.res = fval;
fit_data.fit = x(1)*exp(-(axis-x(2)).^2./(2*((1+x(3)*sign(axis-x(2)))*x(4)).^2))+x(5);

fit_data.rms = rms;
fit_data.prf = prof;
fit_data.axs = axis;