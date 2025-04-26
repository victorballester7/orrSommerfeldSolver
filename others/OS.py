#######Librerias#############
import scipy as sy
from pylab import *
from numpy import *
import scipy.integrate as integ
##############################


def cheb(N):
    if N == 0:
        D = 0.0
        x = 1.0
    else:
        n = arange(0, N + 1)  # genera el vector con paso 1, de 0 hasta N+1
        x = -cos(pi * n / N).reshape(N + 1, 1)  # -1<x<1
        c = (hstack(([2.0], ones(N - 1), [2.0])) * (-1) ** n).reshape(
            N + 1, 1
        )  # reshape(N+1,1) me transpone el vector
        X = tile(
            x, (1, N + 1)
        )  # genera la matriz donde la fila i es el elemento i del vector
        dX = X - X.T
        D = (
            dot(c, 1.0 / c.T) / (dX + eye(N + 1))
        )  # poner dot(c,1./c.T) es lo mismo que  c*(1/c).T (la funcion dot multiplica array)
        D -= diag(sum(D.T, axis=0))
    return D, x.reshape(N + 1)


def Dmat(N):
    D, y = cheb(N)
    D2 = dot(D, D)
    y2 = y * y
    s = hstack(([0.0], 1 / (1 - y2[1:-1]), [0.0]))
    s = diag(s)
    D3 = dot(D2, D)
    D4 = dot(D2, D2)
    D4 = dot((dot(diag(1 - y2), D4) - 8 * dot(diag(y), dot(D2, D)) - 12 * D2), s)
    D = D[1:-1, 1:-1]  # chequear que no cambie nada agregar esto aca
    D2 = D2[1:-1, 1:-1]
    D3 = D3[1:-1, 1:-1]
    D4 = D4[1:-1, 1:-1]

    return D, D2, D3, D4


#############################################
################## Mapeos ###################

## parametros mapeo BL##
z_inf = 25  # donde u es casi 1
zi = 5  # donde queres que esten distribuidos la mitad de los puntos
a = (zi * z_inf) / (z_inf - 2 * zi)
b = 1 + (2 * a) / z_inf


############################################
def mapping(y, Np):
    if Np == 2:  # Boundary Layer
        global a
        global b
        eta = a * ((1.0 + y) / (b - y))
    if Np == 3:
        eta = np.arctanh(y)
    return eta


def d_mapping(y, Np):
    z = mapping(y, Np)
    if Np == 2:
        global a
        global b
        d_eta = a * (b + 1) / (z + a) ** 2
    if Np == 3:
        d_eta = 1 / np.cosh(z) ** 2
    # d_eta=eta_inf/2
    return d_eta


def d2_mapping(y, Np):
    z = mapping(y, Np)
    if Np == 2:
        global a
        global b
        d2_eta = -2 * a * (b + 1) / (z + a) ** 3
    if Np == 3:
        d2_eta = -2 * np.tanh(z) / np.cosh(z) ** 2
    return d2_eta


def d3_mapping(y, Np):
    z = mapping(y, Np)
    if Np == 2:
        global a
        global b
        d3_eta = 6 * a * (b + 1) / (z + a) ** 4
    if Np == 3:
        d3_eta = 4 * np.tanh(z) ** 2 / np.cosh(z) ** 2 - 2 / np.cosh(z) ** 4
    return d3_eta


def d4_mapping(y, Np):
    z = mapping(y, Np)
    if Np == 2:
        global a
        global b
        d4_eta = -24 * a * (b + 1) / (z + a) ** 5
    if Np == 3:
        d4_eta = (
            16 * np.tanh(z) / np.cosh(z) ** 4 - 8 * np.tanh(z) ** 3 / np.cosh(z) ** 2
        )
    return d4_eta


############################################
########### Defino flujos###################
def Poiseuille_Flow(N):
    D, y = cheb(N)
    D2 = dot(D, D)
    y2 = y * y
    U = 1 - y2
    U = U.reshape(N + 1, 1)
    y = y.reshape(N + 1, 1)
    U1 = dot(D, U)
    U2 = dot(D2, U)
    U = U[1:-1]
    U1 = U1[1:-1]
    U2 = U2[1:-1]
    return U, U1, U2


def Couette_Flow(N):
    D, y = cheb(N)
    D2 = dot(D, D)
    U = y
    U = U.reshape(N + 1, 1)
    y = y.reshape(N + 1, 1)
    U1 = dot(D, U)
    U2 = np.zeros(N + 1)
    U = U[1:-1]
    U1 = U1[1:-1]
    U2 = U2[1:-1]
    return U, U1, U2


def Boudary_layer(t, x0):
    f = x0[0]
    f1 = x0[1]
    f2 = x0[2]
    df = f1
    df1 = f2
    df2 = -f * f2 / 2
    return df, df1, df2


def simulate_BL(Boudary_layer, t_span, x0, t):
    global y
    r = integ.solve_ivp(
        Boudary_layer, t_span, x0, t_eval=t, method="Radau"
    )  # Radau es la mejor opcion
    return r.y


def Boundary_Layer_Flow(N, Np):
    _, y = cheb(N)
    y = y[1:-1]
    # print(y)
    t = mapping(y, Np)
    # print(t)

    t_span = (t[0], t[t.size - 1])
    A = 0.332057336
    x0 = [0, 0, A]  # saque ese valor con el metodo del tiro
    f, f1, f2 = simulate_BL(Boudary_layer, t_span, x0, t)
    U = f1
    U1 = f2
    U2 = -f2 * f / 2
    return U, U1, U2


def Jet_Flow(N, Np):
    _, y = cheb(N)
    y = y[1:-1]
    x = mapping(y, Np)
    sech = 1 / np.cosh(x)
    U = sech**2
    U1 = -2 * sech**2 * np.tanh(x)
    U2 = 4 * sech**2 * np.tanh(x) ** 2 - 2 * sech**4
    return U, U1, U2


def flow(N, Np):
    # Np=0 Poiseuille flow;Np=1 Couette flow;Np=2 Boundary layer ;Np=3 Jet flow
    if Np == 0:
        U, U1, U2 = Poiseuille_Flow(N)
    if Np == 1:
        U, U1, U2 = Couette_Flow(N)
    if Np == 2:
        U, U1, U2 = Boundary_Layer_Flow(N, Np)
    if Np == 3:
        U, U1, U2 = Jet_Flow(N, Np)
    return U, U1, U2


##########################################################################
################### Normalización  #######################################
##########################################################################
def normalizacion(u, v, w):
    ur = u.real
    ui = u.imag
    norm = ur * ur + ui * ui
    i = np.where(abs(norm) == np.max(norm))
    ci = 1 / (ur[i] + ui[i] ** 2 / ur[i])
    c = ci - 1j * (ci * ui[i] / ur[i])
    u = c * u
    v = c * v
    w = c * w
    u_check = u.real[: int(u.size / 2)]  ## busco maxima amplitud en -1<x<0
    tol = 1e-2
    if 1 - tol < max(u_check) < 1 + tol:
        return u, v, w
    return -u, -v, -w


##########################################################################
####################### Funciones Principales ############################
##########################################################################
def Orr_Sommerfeld_Temporal(N, R, alp, b, n, Np):
    D, D2, D3, D4 = Dmat(N)
    _, y = cheb(N)
    y = y[1:-1]
    U, U1, U2 = flow(N, Np)
    I = eye(N - 1)
    # print("D")
    # print(D)

    # print("D2")
    # print(D2)

    # print("D4")
    # print(D4)

    D_mapping = D
    D2_mapping = D2
    D4_mapping = D4
    if Np == 2 or Np == 3:
        D_mapping = np.dot(d_mapping(y, Np) * I, D)
        D2_mapping = np.array(
            np.dot(d_mapping(y, Np) ** 2 * I, D2) + np.dot(d2_mapping(y, Np) * I, D)
        )
        D4_mapping = np.array(
            np.dot(d_mapping(y, Np) ** 4 * I, D4)
            + 6 * np.dot(d2_mapping(y, Np) * d_mapping(y, Np) ** 2 * I, D3)
            + 3 * np.dot(d2_mapping(y, Np) ** 2 * I, D2)
            + 4 * np.dot(d3_mapping(y, Np) * d_mapping(y, Np) * I, D2)
            + np.dot(d4_mapping(y, Np) * I, D)
        )
    ra = 1
    if Np == 3:  # jet->rayleigh equation
        ra = 0


    k2 = alp * alp + b * b
    A = np.array(
        -(D4_mapping - 2.0 * k2 * D2_mapping + k2 * k2 * I) / (R) * ra
        - 1j * alp * U2 * I
        + 1j * alp * dot(U * I, D2_mapping)
        - 1j * alp * k2 * U * I
    )
    B = np.array(1j * alp * (D2_mapping - k2 * I))

    # print("U")
    # print(U)

    # print("U2")
    # print(U2)


    print("D2_mapping")
    print(D4_mapping)



    # print("B")
    # print(B)

    # print("their A")
    # print(A)
    # print("their B")
    # print(B)


    H = dot(inv(B), A)
    lam, V = eig(H)

    # print("my lam")
    # print(lam)

    if Np == 2: # Boundary Layer
        # remove all the eigenvalues with real part greater than 0.99
        lam = lam[abs(lam.real) < 0.95]

    ii = argsort(-lam.imag)
    lam = lam[ii]
    print(lam[:10])

    V = V[:, ii]
    ##############calculo perturbaciones##########################
    A = -(D2_mapping - k2 * I) / (R) * ra + 1j * alp * I * U - 1j * lam[n] * alp * I
    v = V[:, n].reshape(N - 1, 1)
    eta = dot(dot(inv(A), U1 * I), -1j * b * v)
    u = (1j * alp * dot(D_mapping, v) - 1j * b * eta) / k2
    w = (1j * b * dot(D_mapping, v) + 1j * alp * eta) / k2
    #############################################################
    u, v, w = normalizacion(u, v, w)
    return lam, u, v, w


def Orr_Sommerfeld_Espacial(N, R, w, b, n, Np):
    D, D2, D3, D4 = Dmat(N)
    _, y = cheb(N)
    y = y[1:-1]
    U, U1, U2 = flow(N, Np)  # perfil(N) #generico

    # k2 = alp * alp + b * b
    b2 = b * b
    b3 = b * b * b
    b4 = b * b * b * b
    I = eye(N - 1)
    O = np.zeros((N - 1, N - 1))

    D_mapping = D
    D2_mapping = D2
    D3_mapping = D3
    D4_mapping = D4
    if Np == 2 or Np == 3:
        D_mapping = np.dot(d_mapping(y, Np) * I, D)
        D2_mapping = np.array(
            np.dot(d_mapping(y, Np) ** 2 * I, D2) + np.dot(d2_mapping(y, Np) * I, D)
        )
        D3_mapping = np.array(
            np.dot(d_mapping(y, Np) ** 3 * I, D3)
            + 3 * np.dot(d2_mapping(y, Np) * d_mapping(y, Np) * I, D2)
            + np.dot(d3_mapping(y, Np) * I, D)
        )
        D4_mapping = np.array(
            np.dot(d_mapping(y, Np) ** 4 * I, D4)
            + 6 * np.dot(d2_mapping(y, Np) * d_mapping(y, Np) ** 2 * I, D3)
            + 3 * np.dot(d2_mapping(y, Np) ** 2 * I, D2)
            + 4 * np.dot(d3_mapping(y, Np) * d_mapping(y, Np) * I, D2)
            + np.dot(d4_mapping(y, Np) * I, D)
        )

    R2 = D2_mapping * 4 / R + 2j * dot(U * I, D_mapping)
    R1 = (
        -2j * w * D_mapping
        - D3_mapping * 4 / R
        + 4 * (b2 / R) * D_mapping
        - 1j * dot(U * I, D2_mapping)
        + (1j * b2) * U * I
        + 1j * U2 * I
    )
    R0 = (
        1j * w * D2_mapping
        - 1j * (w * b2) * I
        + D4_mapping / R
        - 2 * (b2) / R * D2_mapping
        + (b4) / R * I
    )
    T1 = 2 * D_mapping / R + 1j * U * I
    T0 = -1j * w * I - D2_mapping / R + b2 / R * I
    S = 1j * b * U1 * I

    R2 = np.array(R2)
    R1 = np.array(R1)
    R0 = np.array(R0)
    T0 = np.array(T0)
    T0 = np.array(T0)
    S = np.array(S)
    F1 = np.concatenate(([-R1, -R0, O]), axis=1)  # con esto creo las "Filas"
    F2 = np.concatenate(([I, O, O]), axis=1)
    F3 = np.concatenate(([O, -S, -T0]), axis=1)

    # A=[[-R1,-R0,O],[I,O,O],[O,-S,-T0]]#
    A = np.concatenate(([F1, F2, F3]), axis=0)

    F4 = np.concatenate(([R2, O, O]), axis=1)
    F5 = np.concatenate(([O, I, O]), axis=1)
    F6 = np.concatenate(([O, O, T1]), axis=1)

    # B=[[R2,O,O],[O,I,O],[O,O,T1]]
    B = np.concatenate(([F4, F5, F6]), axis=0)

    H = dot(inv(B), A)
    lam, V = eig(H)
    # lam=lam[N-1:2*N]#solo me quedo con el medio
    

    V = V[int(lam.size / 3) : int(2 * lam.size / 3), 0 : int(lam.size / 3)]
    lam = lam[int(lam.size / 3) : int(2 * lam.size / 3)]
    # ii = argsort(-1 / lam.imag)


    # changed by victor (uncomment the above 3 lines if you want to use the original code)
    # filterning the eigenvalues
    # print("before filterning")
    # print(lam)

    lam = lam[lam.imag > 0]

    # print("after filterning")
    # print(lam)

    # lam = lam[lam.real < 1]
    ii = argsort(lam.imag)
    lam = lam[ii]
    # finished changing by victor

    # lam = lam[-ii]
    print(lam[:10])
    # V = V[:, ii]
    ##############calculo perturbaciones##########################
    # A = -(D2_mapping - k2 * I) / (R) + 1j * alp * I * U - 1j * lam[n] * alp * I
    # v = V[:, n].reshape(N - 1, 1)
    # eta = dot(dot(inv(A), U1 * I), -1j * b * v)
    # u = (1j * alp * dot(D_mapping, v) - 1j * b * eta) / k2
    # w = (1j * b * dot(D_mapping, v) + 1j * alp * eta) / k2
    # #############################################################
    # u, v, w = normalizacion(u, v, w)
    # return lam, u, v, w
    return lam
