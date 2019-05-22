function f = fMultisine_time(t, par)
% par_multisine = struct('A',A,'f0',f0,'har',har,'phase',phase);

phase = par.phase;
har = par.har;
f0 = par.f0;
A = par.A;
N = length(har);

f = har'*A*cos(2*pi*[1:N]'*f0*t + phase) /sqrt(sum(har));

% explicit loop, this is how multisine is defined.
% f = 0;
% for i=1:length(har)
%     f = f + har(i)*A*cos(2*pi*i*f0*t + phase(i));
% end
% f = 1/sqrt(sum(har)) * f;
end