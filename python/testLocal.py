from scipy.integrate import tplquad

def quadratic_mass_matrix(x,y,z,index1, index2):
    ksi_0 = (1 - 2 * x) * (1 - x)
    ksi_1 = 4 * x * (1 - x)
    ksi_2 = -x * (1 - 2 * x)

    nu_0 = (1 - 2 * y) * (1 - y)
    nu_1 = 4 * y * (1 - y)
    nu_2 = -y * (1 - 2 * y)

    tetta_0 = (1 - 2 * z) * (1 - z)
    tetta_1 = 4 * z * (1 - z)
    tetta_2 = -z * (1 - 2 * z)

    value1 = 0
    value2 = 0
    
    match index1:
        case 0:
            value1 = ksi_0 * nu_0 * tetta_0
        case 1:
            value1 = ksi_1 * nu_0 * tetta_0
        case 2:
            value1 = ksi_2 * nu_0 * tetta_0
        case 3:
            value1 = ksi_0 * nu_1 * tetta_0
        case 4:
            value1 = ksi_1 * nu_1 * tetta_0
        case 5:
            value1 = ksi_2 * nu_1 * tetta_0
        case 6:
            value1 = ksi_0 * nu_2 * tetta_0
        case 7:
            value1 = ksi_1 * nu_2 * tetta_0
        case 8:
            value1 = ksi_2 * nu_2 * tetta_0
        case 9:
            value1 = 4 * (0.75 - x) * (0.75 - y) * tetta_1
        case 10:
            value1 = -4 * (0.25 - x) * (0.75 - y) * tetta_1
        case 11:
            value1 = -4 * (0.75 - x) * (0.25 - y) * tetta_1
        case 12:
            value1 = 4 * (0.25 - x) * (0.25 - y) * tetta_1
        case 13:
            value1 = tetta_2
            
    match index2:
        case 0:
            value2 = ksi_0 * nu_0 * tetta_0
        case 1:
            value2 = ksi_1 * nu_0 * tetta_0
        case 2:
            value2 = ksi_2 * nu_0 * tetta_0
        case 3:
            value2 = ksi_0 * nu_1 * tetta_0
        case 4:
            value2 = ksi_1 * nu_1 * tetta_0
        case 5:
            value2 = ksi_2 * nu_1 * tetta_0
        case 6:
            value2 = ksi_0 * nu_2 * tetta_0
        case 7:
            value2 = ksi_1 * nu_2 * tetta_0
        case 8:
            value2 = ksi_2 * nu_2 * tetta_0
        case 9:
            value2 = 4 * (0.75 - x) * (0.75 - y) * tetta_1
        case 10:
            value2 = -4 * (0.25 - x) * (0.75 - y) * tetta_1
        case 11:
            value2 = -4 * (0.75 - x) * (0.25 - y) * tetta_1
        case 12:
            value2 = 4 * (0.25 - x) * (0.25 - y) * tetta_1
        case 13:
            value2 = tetta_2
            
    return value1 * value2

def linear_right_part(x,y,z,index1,index2):
    value = 0

    match index1:
        case 0:
            value = (1-x)*(1-y)*(1-z)
        case 1:
            value = x*(1-y)*(1-z)
        case 2:
            value = (1-x)*y*(1-z)
        case 3:
            value = x*y*(1-z)
        case 4:
            value = z
            
    return 5 * value

def quadratic_right_part(x,y,z,index1,index2):
    ksi_0 = (1 - 2 * x) * (1 - x)
    ksi_1 = 4 * x * (1 - x)
    ksi_2 = -x * (1 - 2 * x)

    nu_0 = (1 - 2 * y) * (1 - y)
    nu_1 = 4 * y * (1 - y)
    nu_2 = -y * (1 - 2 * y)

    tetta_0 = (1 - 2 * z) * (1 - z)
    tetta_1 = 4 * z * (1 - z)
    tetta_2 = -z * (1 - 2 * z)

    value = 0

    match index1:
        case 0:
            value = ksi_0 * nu_0 * tetta_0
        case 1:
            value = ksi_1 * nu_0 * tetta_0
        case 2:
            value = ksi_2 * nu_0 * tetta_0
        case 3:
            value = ksi_0 * nu_1 * tetta_0
        case 4:
            value = ksi_1 * nu_1 * tetta_0
        case 5:
            value = ksi_2 * nu_1 * tetta_0
        case 6:
            value = ksi_0 * nu_2 * tetta_0
        case 7:
            value = ksi_1 * nu_2 * tetta_0
        case 8:
            value = ksi_2 * nu_2 * tetta_0
        case 9:
            value = 4 * (0.75 - x) * (0.75 - y) * tetta_1
        case 10:
            value = -4 * (0.25 - x) * (0.75 - y) * tetta_1
        case 11:
            value = -4 * (0.75 - x) * (0.25 - y) * tetta_1
        case 12:
            value = 4 * (0.25 - x) * (0.25 - y) * tetta_1
        case 13:
            value = ksi_1*nu_1*tetta_2
            
    return 5*value
    

def quadratic_mass_matrix(x,y,z,index1, index2):
    ksi_0 = (1 - 2 * x) * (1 - x)
    ksi_1 = 4 * x * (1 - x)
    ksi_2 = -x * (1 - 2 * x)

    nu_0 = (1 - 2 * y) * (1 - y)
    nu_1 = 4 * y * (1 - y)
    nu_2 = -y * (1 - 2 * y)

    tetta_0 = (1 - 2 * z) * (1 - z)
    tetta_1 = 4 * z * (1 - z)
    tetta_2 = -z * (1 - 2 * z)

    value1 = 0
    
    match index1:
        case 0:
            value1 = ksi_0 * nu_0 * tetta_0
        case 1:
            value1 = ksi_1 * nu_0 * tetta_0
        case 2:
            value1 = ksi_2 * nu_0 * tetta_0
        case 3:
            value1 = ksi_0 * nu_1 * tetta_0
        case 4:
            value1 = ksi_1 * nu_1 * tetta_0
        case 5:
            value1 = ksi_2 * nu_1 * tetta_0
        case 6:
            value1 = ksi_0 * nu_2 * tetta_0
        case 7:
            value1 = ksi_1 * nu_2 * tetta_0
        case 8:
            value1 = ksi_2 * nu_2 * tetta_0
        case 9:
            value1 = 4 * (0.75 - x) * (0.75 - y) * tetta_1
        case 10:
            value1 = -4 * (0.25 - x) * (0.75 - y) * tetta_1
        case 11:
            value1 = -4 * (0.75 - x) * (0.25 - y) * tetta_1
        case 12:
            value1 = 4 * (0.25 - x) * (0.25 - y) * tetta_1
        case 13:
            value1 = tetta_2
            
    return value1

def triple_integral(a, b, c, d, e, f, index1, index2, func):
    def integrand(x, y, z):
        return func(x, y, z, index1, index2)
    
    result, error = tplquad(
        func=integrand,
        a=a, b=b,
        gfun=lambda y: c, hfun=lambda y: d,
        qfun=lambda y, z: e, rfun=lambda y, z: f
    )
    return result

def print_quadratic_mass_matrix():
    for i in range(14):
        for j in range(14):
            print(triple_integral(0, 1, 0, 1, 0, 1, i, j, quadratic_mass_matrix), end=' ')
        print()
        
def print_quadratic_right_part():
    for i in range(14):
        print(triple_integral(0, 1, 0, 1, 0, 1, i, 0, quadratic_right_part))

sum = 0
for i in range(14):
    sum += triple_integral(0, 1, 0, 1, 0, 1, i, 0, quadratic_right_part)
# print(sum)
sum = 0
for i in range(5):
    sum += triple_integral(0, 1, 0, 1, 0, 1, i, 0, linear_right_part)
# print(sum)
# print_quadratic_mass_matrix()
print()
# print_quadratic_right_part_local_vector()
print()
# print_quadratic_right_part()

x0, y0, z0 = 0, 0, 0
x1, y1, z1 = 1, 0, 0
x2, y2, z2 = 0, 1, 0
x3, y3, z3 = 1, 1, 0
x4, y4, z4 = 0.5, 0.5, 0.5
xc, yc, zc = (x0 + x1 + x2 + x3) / 4.0, (y0 + y1 + y2 + y3) / 4.0, (z0 + z1 + z2 + z3) / 4.0
D10 = (x1-x0)**2 + (y1-y0)**2 + (z1-z0)**2
D20 = (x2-x0)**2 + (y2-y0)**2 + (z2-z0)**2
D4c = (x4-xc)**2 + (y4-yc)**2 + (z4-zc)**2

dksi_dx = (x1-x0) / D10
dksi_dy = (y1-y0) / D10
dksi_dz = (z1-z0) / D10

dnu_dx = (x2-x0) / D20
dnu_dy = (y2-y0) / D20
dnu_dz = (z2-z0) / D20

dtetta_dx = (x4-xc) / D4c
dtetta_dy = (y4-yc) / D4c
dtetta_dz = (z4-zc) / D4c

def a(ksi,nu,tetta, index1, index2):
    ksi_0 = (1 - 2 * ksi) * (1 - ksi)
    ksi_1 = 4 * ksi * (1 - ksi)
    ksi_2 = -ksi * (1 - 2 * ksi)
    
    dksi_0 = 4 * ksi - 3
    dksi_1 = 4 - 8 * ksi
    dksi_2 = -1 + 4 * ksi

    nu_0 = (1 - 2 * nu) * (1 - nu)
    nu_1 = 4 * nu * (1 - nu)
    nu_2 = -nu * (1 - 2 * nu)

    dnu_0 = 4 * nu - 3
    dnu_1 = 4 - 8 * nu
    dnu_2 = -1 + 4 * nu

    tetta_0 = (1 - 2 * tetta) * (1 - tetta)
    tetta_1 = 4 * tetta * (1 - tetta)
    tetta_2 = -tetta * (1 - 2 * tetta)
    
    dtetta_0 = 4 * tetta - 3
    dtetta_1 = 4 - 8 * tetta
    dtetta_2 = -1 + 4 * tetta
    
    match index1:
        case 0:
            dphi_dksi1   =dksi_0 * nu_0 * tetta_0
            dphi_dnu1    =ksi_0 * dnu_0 * tetta_0
            dphi_dtetta1 =ksi_0 * nu_0 * dtetta_0
        case 1:
            dphi_dksi1   =dksi_1 * nu_0 * tetta_0
            dphi_dnu1    =ksi_1 * dnu_0 * tetta_0
            dphi_dtetta1 =ksi_1 * nu_0 * dtetta_0
        case 2:
            dphi_dksi1   =dksi_2 * nu_0 * tetta_0
            dphi_dnu1    =ksi_2 * dnu_0 * tetta_0
            dphi_dtetta1 =ksi_2 * nu_0 * dtetta_0
        case 3:
            dphi_dksi1   =dksi_0 * nu_1 * tetta_0
            dphi_dnu1    =ksi_0 * dnu_1 * tetta_0
            dphi_dtetta1 =ksi_0 * nu_1 * dtetta_0
        case 4:
            dphi_dksi1   =dksi_1 * nu_1 * tetta_0
            dphi_dnu1    =ksi_1 * dnu_1 * tetta_0
            dphi_dtetta1 =ksi_1 * nu_1 * dtetta_0
        case 5:
            dphi_dksi1   =dksi_2 * nu_1 * tetta_0
            dphi_dnu1    =ksi_2 * dnu_1 * tetta_0
            dphi_dtetta1 =ksi_2 * nu_1 * dtetta_0
        case 6:
            dphi_dksi1   =dksi_0 * nu_2 * tetta_0
            dphi_dnu1    =ksi_0 * dnu_2 * tetta_0
            dphi_dtetta1 =ksi_0 * nu_2 * dtetta_0
        case 7:
            dphi_dksi1   =dksi_1 * nu_2 * tetta_0
            dphi_dnu1    =ksi_1 * dnu_2 * tetta_0
            dphi_dtetta1 =ksi_1 * nu_2 * dtetta_0
        case 8:
            dphi_dksi1   =dksi_2 * nu_2 * tetta_0
            dphi_dnu1    =ksi_2 * dnu_2 * tetta_0
            dphi_dtetta1 =ksi_2 * nu_2 * dtetta_0
        case 9:
            dphi_dksi1   =-4 * (0.75 - nu) * tetta_1
            dphi_dnu1    =-4 * (0.75 - ksi) * tetta_1
            dphi_dtetta1 =4 * (0.75 - ksi) * (0.75 - nu) * dtetta_1
        case 10:
            dphi_dksi1   =4 * (0.75 - nu) * tetta_1
            dphi_dnu1    =4 * (0.25 - ksi) * tetta_1
            dphi_dtetta1 =-4 * (0.25 - ksi) * (0.75 - nu) * dtetta_1
        case 11:
            dphi_dksi1   =4 * (0.25 - nu) * tetta_1
            dphi_dnu1    =4 * (0.75 - ksi) * tetta_1
            dphi_dtetta1 =-4 * (0.75 - ksi) * (0.25 - nu) * dtetta_1
        case 12:
            dphi_dksi1   =-4 * (0.25 - nu) * tetta_1
            dphi_dnu1    =-4 * (0.25 - ksi) * tetta_1
            dphi_dtetta1 =4 * (0.25 - ksi) * (0.25 - nu) * dtetta_1
        case 13:
            dphi_dksi1   =0
            dphi_dnu1    =0
            dphi_dtetta1 =dtetta_2

    match index2:
        case 0:
            dphi_dksi2   =dksi_0 * nu_0 * tetta_0
            dphi_dnu2   =ksi_0 * dnu_0 * tetta_0
            dphi_dtetta2 =ksi_0 * nu_0 * dtetta_0
        case 1:
            dphi_dksi2   =dksi_1 * nu_0 * tetta_0
            dphi_dnu2    =ksi_1 * dnu_0 * tetta_0
            dphi_dtetta2 =ksi_1 * nu_0 * dtetta_0
        case 2:
            dphi_dksi2   =dksi_2 * nu_0 * tetta_0
            dphi_dnu2    =ksi_2 * dnu_0 * tetta_0
            dphi_dtetta2 =ksi_2 * nu_0 * dtetta_0
        case 3:
            dphi_dksi2   =dksi_0 * nu_1 * tetta_0
            dphi_dnu2    =ksi_0 * dnu_1 * tetta_0
            dphi_dtetta2 =ksi_0 * nu_1 * dtetta_0
        case 4:
            dphi_dksi2   =dksi_1 * nu_1 * tetta_0
            dphi_dnu2    =ksi_1 * dnu_1 * tetta_0
            dphi_dtetta2 =ksi_1 * nu_1 * dtetta_0
        case 5:
            dphi_dksi2   =dksi_2 * nu_1 * tetta_0
            dphi_dnu2    =ksi_2 * dnu_1 * tetta_0
            dphi_dtetta2 =ksi_2 * nu_1 * dtetta_0
        case 6:
            dphi_dksi2   =dksi_0 * nu_2 * tetta_0
            dphi_dnu2    =ksi_0 * dnu_2 * tetta_0
            dphi_dtetta2 =ksi_0 * nu_2 * dtetta_0
        case 7:
            dphi_dksi2   =dksi_1 * nu_2 * tetta_0
            dphi_dnu2    =ksi_1 * dnu_2 * tetta_0
            dphi_dtetta2 =ksi_1 * nu_2 * dtetta_0
        case 8:
            dphi_dksi2   =dksi_2 * nu_2 * tetta_0
            dphi_dnu2    =ksi_2 * dnu_2 * tetta_0
            dphi_dtetta2 =ksi_2 * nu_2 * dtetta_0
        case 9:
            dphi_dksi2   =-4 * (0.75 - nu) * tetta_1
            dphi_dnu2    =-4 * (0.75 - ksi) * tetta_1
            dphi_dtetta2 =4 * (0.75 - ksi) * (0.75 - nu) * dtetta_1
        case 10:
            dphi_dksi2   =4 * (0.75 - nu) * tetta_1
            dphi_dnu2    =4 * (0.25 - ksi) * tetta_1
            dphi_dtetta2 =-4 * (0.25 - ksi) * (0.75 - nu) * dtetta_1
        case 11:
            dphi_dksi2   =4 * (0.25 - nu) * tetta_1
            dphi_dnu2    =4 * (0.75 - ksi) * tetta_1
            dphi_dtetta2 =-4 * (0.75 - ksi) * (0.25 - nu) * dtetta_1
        case 12:
            dphi_dksi2   =-4 * (0.25 - nu) * tetta_1
            dphi_dnu2    =-4 * (0.25 - ksi) * tetta_1
            dphi_dtetta2 =4 * (0.25 - ksi) * (0.25 - nu) * dtetta_1
        case 13:
            dphi_dksi2   =0
            dphi_dnu2    =0
            dphi_dtetta2 =dtetta_2

    dphi_dx1 = dphi_dksi1 * dksi_dx + dphi_dnu1 * dnu_dx + dphi_dtetta1 * dtetta_dx
    dphi_dy1 = dphi_dksi1 * dksi_dy + dphi_dnu1 * dnu_dy + dphi_dtetta1 * dtetta_dy
    dphi_dz1 = dphi_dksi1 * dksi_dz + dphi_dnu1 * dnu_dz + dphi_dtetta1 * dtetta_dz
    
    dphi_dx2 = dphi_dksi2 * dksi_dx + dphi_dnu2 * dnu_dx + dphi_dtetta2 * dtetta_dx
    dphi_dy2 = dphi_dksi2 * dksi_dy + dphi_dnu2 * dnu_dy + dphi_dtetta2 * dtetta_dy
    dphi_dz2 = dphi_dksi2 * dksi_dz + dphi_dnu2 * dnu_dz + dphi_dtetta2 * dtetta_dz
  
    return dphi_dx1 * dphi_dx2 + dphi_dy1 * dphi_dy2 + dphi_dz1 * dphi_dz2

def psi00_0(x,y,z,index1,index2):
    return 2*((4*x-3)**2 * (1-2*y)**2 * (1-y)**2 * (1-2*z)**2 * (1-z)**2 + \
        (1-2*x)**2 * (1-x)**2 * (4*y-3)**2 * (1-2*z)**2 * (1-z)**2 + \
        0.25 * (1-2*x)**2 * (1-x)**2 * (1-2*y)**2 * (1-y)**2 * (4*z-3)**2)
def psi01_0(x,y,z,index1,index2):
    return 2*((4*x-3) * (4-8*x) * (1-2*y)**2 * (1-y)**2 * (1-2*z)**2 * (1-z)**2 + \
        (1-2*x) * (1-x)**2 * 4 * x * (4*y-3)**2 * (1-2*z)**2 * (1-z)**2 + \
        0.25 * (1-2*x) * (1-x)**2 * 4 * x * (1-2*y)**2 * (1-y)**2 * (4*z-3)**2)
def psi00_1(x,y,z,index1,index2):
    return 2*((4*x-3)**2 * (1-2*y)**2 * (1-y)**2 * (1-2*z)**2 * (1-z)**2 + \
        0.25 * (1-2*x)**2 * (1-x)**2 * (4*y-3)**2 * (1-2*z)**2 * (1-z)**2 + \
        (1-2*x)**2 * (1-x)**2 * (1-2*y)**2 * (1-y)**2 * (4*z-3)**2)
def psi01_1(x,y,z,index1,index2):
    return 2*((4*x-3) * (4-8*x) * (1-2*y)**2 * (1-y)**2 * (1-2*z)**2 * (1-z)**2 + \
        0.25 * (1-2*x) * (1-x)**2 * 4 * x * (4*y-3)**2 * (1-2*z)**2 * (1-z)**2 + \
        (1-2*x) * (1-x)**2 * 4 * x * (1-2*y)**2 * (1-y)**2 * (4*z-3)**2)
def psi00_2(x,y,z,index1,index2):
    return 2*(0.25 * (4*x-3)**2 * (1-2*y)**2 * (1-y)**2 * (1-2*z)**2 * (1-z)**2 + \
        (1-2*x)**2 * (1-x)**2 * (4*y-3)**2 * (1-2*z)**2 * (1-z)**2 + \
        (1-2*x)**2 * (1-x)**2 * (1-2*y)**2 * (1-y)**2 * (4*z-3)**2)
def psi01_2(x,y,z,index1,index2):
    return 2*(0.25 * (4*x-3) * (4-8*x) * (1-2*y)**2 * (1-y)**2 * (1-2*z)**2 * (1-z)**2 + \
        (1-2*x) * (1-x)**2 * 4 * x * (4*y-3)**2 * (1-2*z)**2 * (1-z)**2 + \
        (1-2*x) * (1-x)**2 * 4 * x * (1-2*y)**2 * (1-y)**2 * (4*z-3)**2)

# for i in range(14):
#     for j in range(14):
#         print(triple_integral(0, 1, 0, 1, 0, 1, i, j, a), end=' ')
#     print()

print(triple_integral(0, 1, 0, 1, 0, 1, 0, 0, psi00_0))
print(triple_integral(0, 1, 0, 1, 0, 1, 0, 0, psi00_1))
print(triple_integral(0, 1, 0, 1, 0, 1, 0, 0, psi00_2))
print()
print(triple_integral(0, 1, 0, 1, 0, 1, 0, 0, psi01_0))
print(triple_integral(0, 1, 0, 1, 0, 1, 0, 0, psi01_1))
print(triple_integral(0, 1, 0, 1, 0, 1, 0, 0, psi01_2))