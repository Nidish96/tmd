function Y = ode2_modif(odefun,h,y0,N,varargin)
% https://en.wikipedia.org/wiki/Heun%27s_method
% Heun's method may refer to the improved or modified Euler's method (that
% is, the explicit trapezoidal rule), or a similar two-stage Runge–Kutta
% method.
%
% NB:
% calls odefun with the time index, not the time. Modify odefun accordingly

if ~isnumeric(y0)
  error('Y0 should be a vector of initial conditions.');
end

y0 = y0(:);   % Make a column vector.
neq = length(y0);
Y = zeros(neq,N);
F = zeros(neq,2);

Y(:,1) = y0;
for i = 2:N
  yi = Y(:,i-1);
  F(:,1) = feval(odefun,i-1,yi,varargin{:});
  F(:,2) = feval(odefun,i,yi+h*F(:,1),varargin{:});
  Y(:,i) = yi + (h/2)*(F(:,1) + F(:,2));
end
Y = Y.';

end