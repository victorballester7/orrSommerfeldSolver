function [L, M] = OS_2D(Ny, Re, kx, uDat)
    % OS_2D - Constructs matrices for the Orr-Sommerfeld problem in 2D.
    % 
    % Inputs:
    %   Ny         - Number of grid points in the y-direction.
    %   Re         - Reynolds number.
    %   kx         - Streamwise wavenumber.
    %   uDat       - Velocity profile and its derivatives: [U, U', U''].
    %
    % Outputs:
    %   L - Matrix defining the Orr-Sommerfeld system 
    %   M - Mass matrix for the Orr-Sommerfeld system
    %   L * q = i * omega * M * q

    % Interior grid points count
    Nint = Ny - 2;

    % Extract velocity profiles and their derivatives at interior points
    U = uDat(2:end-1, 1);
    Upp = uDat(2:end-1, 3);

    % Define squared wavenumber and identity matrix
    k2 = kx^2;
    Iint = eye(Nint);         % Identity matrix for interior points
    Zint = zeros(Nint);       % Zero matrix for interior points

    % Generate Chebyshev differentiation matrices
    [~, D2] = chebdif(Ny, 2); % 1st and 2nd derivatives
    [~, D4] = cheb4c(Ny);     % 4th derivative with clamped BCs

    % Extract submatrices for interior points with Dirichlet BCs
    D1Dir = D2(2:end-1, 2:end-1, 1); % First derivative
    D2Dir = D2(2:end-1, 2:end-1, 2); % Second derivative
    clearvars D2

    % Define Laplacian and bi-Laplacian operators
    LapDir = D2Dir - k2 * Iint;                  % Laplacian with Dirichlet BCs
    Lap2 = D4 - 2 * k2 * D2Dir + k2^2 * Iint;    % Bi-Laplacian with clamped BCs

    % Define imaginary unit
    i = sqrt(-1);

    % Define Orr-Sommerfeld operator
    L = -i * kx * (diag(U) * LapDir - diag(Upp)) + 1 / Re * Lap2;

    % Define mass matrix (M) for the system
    M = -LapDir;
end

