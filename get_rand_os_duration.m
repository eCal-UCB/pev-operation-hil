function [overstay_endtime, duration] = get_rand_os_duration(opt)
% TODO: regularization parameter for overstay must be the actual overstay
% duration without penalty
par = get_glob_par();
prb = get_glob_prb();
% lambda = par.lambda.h_c * prb.user.overstay_duration / opt.tariff.overstay;
lambda = prb.user.overstay_duration / opt.tariff.overstay;
range = 0:100;
pdf = exp(-lambda).*(lambda).^range./factorial(range);
cdf = cumsum(pdf);
r = (1-min(cdf))*rand + min(cdf);
duration = 0;
try
    duration = interp1(cdf,range,r);
catch
    duration = range(find(cdf>=r,1));
end
overstay_endtime = opt.time.start + duration;
% disp(['duration: ' num2str(duration)]);
if isnan(duration)
    error('[ ERROR] nan duration');
end
end