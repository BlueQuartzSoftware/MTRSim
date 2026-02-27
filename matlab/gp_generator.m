function gp = gp_generator(hx, hy, hz, rho, theta, nx, ny, nz)

%% Compute covariance matrices x direction

x = [hx:hx:(nx-1)*hx]';

Gamma_x = zeros(nx);
x_vec = rho(x, theta(1));
for k = 1:nx-1
	Gamma_x(k,k+1:end) = x_vec(1:nx-k,1)';
end
Gamma_x = Gamma_x + Gamma_x' + eye(nx);
Gamma_x = Gamma_x + 1e-6*eye(nx);  % add a small multiple of the identity to ensure positive definiteness

Rx = chol(Gamma_x);  % Cholesky factor of Gamma_x

%% Compute covariance matrices y direction

y = [hy:hy:(ny-1)*hy]';

Gamma_y = zeros(ny);
y_vec = rho(y, theta(2));
for k = 1:ny-1
	Gamma_y(k,k+1:end) = y_vec(1:ny-k,1)';
end
Gamma_y = Gamma_y + Gamma_y' + eye(ny);
Gamma_y = Gamma_y + 1e-6*eye(ny);  % add a small multiple of the identity to ensure positive definiteness

Ry = chol(Gamma_y);  % Cholesky factor of Gamma_y


%% Compute covariance matrices z direction

z = [hz:hz:(nz-1)*hz]';

Gamma_z = zeros(nz);
z_vec = rho(z, theta(3));
for k = 1:nz-1
	Gamma_z(k,k+1:end) = z_vec(1:nz-k,1)';
end
Gamma_z = Gamma_z + Gamma_z' + eye(nz);
Gamma_z = Gamma_z + 1e-6*eye(nz);  % add a small multiple of the identity to ensure positive definiteness

Rz = chol(Gamma_z);  % Cholesky factor of Gamma_z


%% Random draw - note that this relies on the property that kron(A,B)*vec(X) = B*X*A'. This property is applied twice.

z = randn(nx*ny*nz,1);
W_1 = reshape(z, nx*ny, nz);
u = W_1*Rz;  % multiplication by the cholesky factor in the z-direction first
W_2 = reshape(u, ny, nx, nz);
x = W_2;
for k = 1:nz
	x(:,:,k) = Ry'*(W_2(:,:,k)*Rx);  % multiplication by x and y cholesky factors, slice by slice
end

gp = x;
