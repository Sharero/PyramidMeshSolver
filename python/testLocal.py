from scipy.integrate import tplquad

def stiffness00(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = pow(-(1 - nu) * (1 - tetta) / k, 2)
    term2 = pow(-(1 - ksi) * (1 - tetta) / i, 2)
    term3 = pow(-(1 - nu) * (1 - ksi) / j, 2)
    
    return term1 + term2 + term3
def stiffness01(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = (-(1 - nu) * (1 - tetta) / k) * ((1 - nu) * (1 - tetta) / k)
    term2 = (-(1 - ksi) * (1 - tetta) / i) * (-ksi * (1 - tetta) / i)
    term3 = (-(1 - nu) * (1 - ksi) / j) * (-(1 - nu) * ksi / j)
    
    return term1 + term2 + term3
def stiffness02(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = (-(1 - nu) * (1 - tetta) / k) * (-nu * (1 - tetta) / k)
    term2 = (-(1 - ksi) * (1 - tetta) / i) * ((1 - ksi) * (1 - tetta) / i)
    term3 = (-(1 - nu) * (1 - ksi) / j) * (-nu * (1 - ksi) / j)
    
    return term1 + term2 + term3
def stiffness03(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = (-(1 - nu) * (1 - tetta) / k) * (nu * (1 - tetta) / k)
    term2 = (-(1 - ksi) * (1 - tetta) / i) * (ksi * (1 - tetta) / i)
    term3 = (-(1 - nu) * (1 - ksi) / j) * (-nu * ksi / j)
    
    return term1 + term2 + term3
def stiffness04(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = (-(1 - nu) * (1 - tetta) / k) * (0.0)
    term2 = (-(1 - ksi) * (1 - tetta) / i) * (0.0)
    term3 = (-(1 - nu) * (1 - ksi) / j) * (1.0 / j)
    
    return term1 + term2 + term3

def stiffness11(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = pow((1 - nu) * (1 - tetta) / k, 2)
    term2 = pow(-ksi * (1 - tetta) / i, 2)
    term3 = pow(-(1 - nu) * ksi / j, 2)
    
    return term1 + term2 + term3
def stiffness12(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = ((1 - nu) * (1 - tetta) / k) * (-nu * (1 - tetta) / k)
    term2 = (-ksi * (1 - tetta) / i) * ((1 - ksi) * (1 - tetta) / i)
    term3 = (-(1 - nu) * ksi / j) * (-nu * (1 - ksi) / j)
    
    return term1 + term2 + term3
def stiffness13(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = ((1 - nu) * (1 - tetta) / k) * (nu * (1 - tetta) / k)
    term2 = (-ksi * (1 - tetta) / i) * (ksi * (1 - tetta) / i)
    term3 = (-(1 - nu) * ksi / j) * (-nu * ksi / j)
    
    return term1 + term2 + term3
def stiffness14(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = ((1 - nu) * (1 - tetta) / k) * (0.0)
    term2 = (-ksi * (1 - tetta) / i) * (0.0)
    term3 = (-(1 - nu) * ksi / j) * (1.0 / j)
    
    return term1 + term2 + term3

def stiffness22(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = pow(-nu * (1 - tetta) / k, 2)
    term2 = pow((1 - ksi) * (1 - tetta) / i, 2)
    term3 = pow(-nu * (1 - ksi) / j, 2)
    
    return term1 + term2 + term3
def stiffness23(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = (-nu * (1 - tetta) / k) * (nu * (1 - tetta) / k)
    term2 = ((1 - ksi) * (1 - tetta) / i) * (ksi * (1 - tetta) / i)
    term3 = (-nu * (1 - ksi) / j) * (-nu * ksi / j)
    
    return term1 + term2 + term3
def stiffness24(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = (-nu * (1 - tetta) / k) * (0.0)
    term2 = ((1 - ksi) * (1 - tetta) / i) * (0.0)
    term3 = (-nu * (1 - ksi) / j) * (1.0 / j)
    
    return term1 + term2 + term3

def stiffness33(x,y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = pow(nu * (1 - tetta) / k, 2)
    term2 = pow(ksi * (1 - tetta) / i, 2)
    term3 = pow(-nu * ksi / j, 2)
    
    return term1 + term2 + term3
def stiffness34(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = (nu * (1 - tetta) / k) * (0.0)
    term2 = (ksi * (1 - tetta) / i) * (0.0)
    term3 = (-nu * ksi / j) * (1.0 / j)
    
    return term1 + term2 + term3

def stiffness44(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = pow(0.0, 2)
    term2 = pow(0.0, 2)
    term3 = pow(1.0 / j, 2)
    
    return term1 + term2 + term3

def mass00(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = (1 - ksi) * (1 - nu) * (1 - tetta)
    term2 = (1 - ksi) * (1 - nu) * (1 - tetta)
    
    return term1 * term2
def mass01(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = (1 - ksi) * (1 - nu) * (1 - tetta)
    term2 = ksi * (1 - nu) * (1 - tetta)
    
    return term1 * term2
def mass02(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = (1 - ksi) * (1 - nu) * (1 - tetta)
    term2 = (1 - ksi) * nu * (1 - tetta)
    
    return term1 * term2
def mass03(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = (1 - ksi) * (1 - nu) * (1 - tetta)
    term2 = ksi * nu * (1 - tetta)
    
    return term1 * term2
def mass04(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = (1 - ksi) * (1 - nu) * (1 - tetta)
    term2 = tetta
    
    return term1 * term2

def mass11(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = ksi * (1 - nu) * (1 - tetta)
    term2 = ksi * (1 - nu) * (1 - tetta)
    
    return term1 * term2
def mass12(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = ksi * (1 - nu) * (1 - tetta)
    term2 = (1 - ksi) * nu * (1 - tetta)
    
    return term1 * term2
def mass13(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = ksi * (1 - nu) * (1 - tetta)
    term2 = ksi * nu * (1 - tetta)
    
    return term1 * term2
def mass14(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = ksi * (1 - nu) * (1 - tetta)
    term2 = tetta
    
    return term1 * term2

def mass22(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = (1 - ksi) * nu * (1 - tetta)
    term2 = (1 - ksi) * nu * (1 - tetta)
    
    return term1 * term2
def mass23(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = (1 - ksi) * nu * (1 - tetta)
    term2 = ksi * nu * (1 - tetta)
    
    return term1 * term2
def mass24(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = (1 - ksi) * nu * (1 - tetta)
    term2 = tetta
    
    return term1 * term2

def mass33(x,y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = ksi * nu * (1 - tetta)
    term2 = ksi * nu * (1 - tetta)
    
    return term1 * term2
def mass34(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = ksi * nu * (1 - tetta)
    term2 = tetta
    
    return term1 * term2

def mass44(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = tetta
    term2 = tetta
    
    return term1 * term2

def right_part0(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = x
    term2 = (1 - ksi) * (1 - nu) * (1 - tetta)
    
    return term1 * term2
def right_part1(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = x
    term2 = ksi * (1 - nu) * (1 - tetta)
    
    return term1 * term2
def right_part2(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = x
    term2 = (1 - ksi) * nu * (1 - tetta)
    
    return term1 * term2
def right_part3(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = x
    term2 = ksi * nu * (1 - tetta)
    
    return term1 * term2
def right_part4(x, y, z, m, n, l, k, i, j):
    ksi = (x - m) / k
    nu = (y - n) / i
    tetta = (z - l) / j

    term1 = x
    term2 = tetta
    
    return term1 * term2

def triple_integral(a, b, c, d, e, f, m, n, l, k, i, j, func):
    """
    Вычисление тройного интеграла с передаваемой функцией.
    
    Параметры:
    a, b: пределы интегрирования по x
    c, d: пределы интегрирования по y
    e, f: пределы интегрирования по z
    m, n, l: параметры смещения (центры)
    k, i, j: масштабные параметры
    func: подынтегральная функция, принимающая x, y, z и параметры
    
    Возвращает:
    Значение тройного интеграла
    """
    # Подынтегральная функция с фиксированными параметрами
    def integrand(z, y, x):
        return func(x, y, z, m, n, l, k, i, j)
    
    # Вычисление интеграла с помощью tplquad
    result, error = tplquad(
        func=integrand,   # подынтегральная функция
        a=a, b=b,         # пределы интегрирования по x
        gfun=lambda x: c, hfun=lambda x: d,   # пределы по y
        qfun=lambda x, y: e, rfun=lambda x, y: f  # пределы по z
    )
    
    return result

def print_stiffness_local_matrix():
    G00 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness00)
    G01 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness01)
    G02 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness02)
    G03 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness03)
    G04 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness04)

    G10 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness01)
    G11 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness11)
    G12 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness12)
    G13 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness13)
    G14 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness14)

    G20 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness02)
    G21 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness12)
    G22 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness22)
    G23 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness23)
    G24 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness24)

    G30 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness03)
    G31 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness13)
    G32 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness23)
    G33 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness33)
    G34 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness34)

    G40 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness04)
    G41 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness14)
    G42 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness24)
    G43 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness34)
    G44 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, stiffness44)

    print(G00, G01, G02, G03, G04)
    print(G10, G11, G12, G13, G14)
    print(G20, G21, G22, G23, G24)
    print(G30, G31, G32, G33, G34)
    print(G40, G41, G42, G43, G44)

def print_mass_matrix_local_matrix():
    M00 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass00)
    M01 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass01)
    M02 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass02)
    M03 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass03)
    M04 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass04)

    M10 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass01)
    M11 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass11)
    M12 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass12)
    M13 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass13)
    M14 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass14)

    M20 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass02)
    M21 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass12)
    M22 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass22)
    M23 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass23)
    M24 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass24)

    M30 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass03)
    M31 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass13)
    M32 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass23)
    M33 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass33)
    M34 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass34)

    M40 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass04)
    M41 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass14)
    M42 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass24)
    M43 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass34)
    M44 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, mass44)

    print(M00, M01, M02, M03, M04)
    print(M10, M11, M12, M13, M14)
    print(M20, M21, M22, M23, M24)
    print(M30, M31, M32, M33, M34)
    print(M40, M41, M42, M43, M44)

def print_right_parn_local_vector():
    f0 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, right_part0)
    f1 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, right_part1)
    f2 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, right_part2)
    f3 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, right_part3)
    f4 = triple_integral(a, b, c, d, e, f, xk, yk, zk, hx, hy, hz, right_part4)
    
    print(f0, f1, f2, f3, f4)

xk, yk, zk = 0, 0, 0
hx, hy, hz = 2, 2, 4
a, b = 0, 2
c, d = 0, 2
e, f = 0, 4

# print_stiffness_local_matrix()
# print_mass_matrix_local_matrix()
# print_right_parn_local_vector()