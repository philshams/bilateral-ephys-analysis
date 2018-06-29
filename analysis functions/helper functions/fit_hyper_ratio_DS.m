function [ err, pars ] = fit_hyper_ratio_DS(cs,resps,nn, Rmax, sigma, n, R0)
% FIT_HYPER_RATIO_DS fits hyper_ratio to the data
%
% 	[ err, pars ] = fit_hyper_ratio(cs,resps)
% 	[ err, pars ] = fit_hyper_ratio(cs,resps,nn) uses nn starting points (Default:3)
% 	[ err, pars ] = fit_hyper_ratio(cs,resps,nn,Rmax, sigma, n, R0) imposes each paramters if not empty
%
%  hyper_ratio pars are [ Rmax, sigma, n, R0 ]
%
% 
% 4/4/18 created from fit_hyper_ratio
% 5/4/18 fit/evaluate using bins without nan/inf
%
% part of the dsbox
%
% see also: hyper_ratio


if nargin < 3 || isempty(nn)
   nn = 3;
end



cs = cs(:);
resps = resps(:);

if any(~isfinite(cs)) || any(~isfinite(resps))
   okidx = intersect(find(isfinite(cs)), find(isfinite(resps)));
   cs = cs(okidx);
   resps = resps(okidx);
    %error('yoooooooo');  
end

if size(cs)~=size(resps)
   error('size of contrasts and responses not equal');
end

% -------------- initial and min/max values % added on 4/4/18

Rmax_i =  max(resps)-min(resps); %19/4/18

if nargin < 4 || isempty(Rmax)
    Rmax_range = [0 Rmax_i, 2*Rmax_i]; %[min init max]
else
    Rmax_range = [Rmax Rmax Rmax];
end

if  nargin < 5 || isempty(sigma)
    sigma_range = [eps mean(cs) 10*max(cs)]; %23/4/18
else
    sigma_range = [sigma sigma sigma];
end

if  nargin < 6 || isempty(n)
    n_range = [0 2.5 10];%[0 2.5 5];
else
    n_range = [n n n];
end

if  nargin < 7 || isempty(R0)
    if any(cs==0)
      R0_i = mean(resps(cs==0));
   else 
      R0_i = 0;
   end
    R0_range = [min(resps) R0_i Rmax_i+2*R0_i]; %19/4/18 minimum changed from 0 to min(resps)
else
    R0_range = [R0 R0 R0];
end
[ err, pars ] = fitit('hyper_ratio', resps,...
    [Rmax_range(1) sigma_range(1) n_range(1) R0_range(1)], ...
    [Rmax_range(2) sigma_range(2) n_range(2) R0_range(2)], ...
    [Rmax_range(3) sigma_range(3) n_range(3) R0_range(3)], ...
    [0 1e-4 1e-4 nn], cs );


% figure;
% plot(cs,resps,'o');
% hold on
% cc = linspace(min(cs),max(cs));
% plot(cc,hyper_ratio(pars,cc),'k-');
