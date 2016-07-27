function [params, fit] = gFit(data_ax,data,guess)

% inputs are data_ax = x axis, data = gaussian to be fit
% guess(1) = amplitude
% guess(2) = mu
% guess(3) = sigma
% guess(4) = offset
% fit_func = guess(1)*exp(-(data_ax - guess(2))^2/(2*guess(3)^3))+guess(4)


params = fminsearch(@(params) gauss(params,data_ax,data),guess);

fit = params(1)*exp(-(data_ax - params(2)).^2/(2*params(3)^2))+params(4);


%ARW = fminsearch(@(z) F_error(lambda,z(1),z(2),z(3)),[2;500;sqrt(140)]);

function fit_error = gauss(guess,data_ax,data)

fit_func = guess(1)*exp(-(data_ax - guess(2)).^2/(2*guess(3)^2))+guess(4);

fit_error = sum(data - fit_func).^2;