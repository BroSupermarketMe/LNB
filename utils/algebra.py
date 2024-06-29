import math,cmath

# create three-dimensional coordinates with uniform distribution on the sphere 创建球面上均匀随机分布的三维坐标
def balance(r,a):
    q = 34
    s = 1/(math.sqrt(5)*q)
    pi = (math.sqrt(5)-1)/2
    coord = []
    n= int(4*math.pi*pow(r,2)/a)

    for i in range(1, n):
        z = r*(2*i-1)/(n-1)
        az = abs(z-r)
        y = math.sqrt(r**2 - az**2)*math.cos(2*math.pi*i*pi)
        x = math.sqrt(r**2 - az**2)*math.sin(2*math.pi*i*pi)
        z = z -r
        coord.append([x, y, z])
    return coord

# Convert cartesian coordinates to spherical coordinates 直角坐标转球坐标
def trans1(a, r):
    coordi = []
    for i in a:
        x = i[0]
        y = i[1]
        z = i[2]
        theta = math.acos(z/r)
        phi = math.atan(y/x)
        if x > 0 and y > 0:
            phi = phi
        if x < 0 and y > 0:
            phi = phi + math.pi
        if x < 0 and y < 0:
            phi = phi + math.pi
        if x > 0 and y < 0:
            phi = phi + 2*math.pi
        coordi.append([phi, theta])
    return coordi

# Spherical pole projection: turn spherical coordinates to plane cartesian coordinates
# 球坐标投影到平面：将球坐标转换为平面直角坐标
def trans2(r, a):
    i = (-1)**0.5
    c = []
    for j in a:
        phi = j[0]
        theta = j[1]
        z = r*1/math.tan(theta/2)*cmath.exp(i*phi)
        x = z.real
        y = z.imag
        c.append([x, y])
    return c

# Reverse projection of plane cartesian coordinates into spherical coordinates 平面直角坐标反投影到球坐标
def trans3(z, r):
    zz = cmath.polar(z)
    zr = zz[0]
    phi = zz[1]
    theta = 2*math.atan(r/zr)
    return(phi,theta)


# Create a plural x,y坐标转复数
def trans4(x,y):
    z = complex(x,y)
    return z

# Spherical coordinates to cartesian coordinates 球坐标转直角坐标
def trans5(a,r):

    phi = a[0]
    theta = a[1]

    x = r*math.sin(theta)*math.cos(phi)
    y = r*math.sin(theta)*math.sin(phi)
    z = r*math.cos(theta)

    c = (x,y,z)

    return c

# the shortest distance between the two points 两点间最短距离
def distance(a,b):
    '''

    :param a: [x1,y1,z1]
    :param b: [x2,y2,z2]
    :return:
    '''
    d = math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)
    return d

def vector(v):
    '''

    :param v: str or float
              'x,y,z' or float or 'x'
    :return: [float,float,float] or float
    '''
    if type(v) == str and "," in v:
        return [float(i) for i in v.split(",")]
    return float(v)

def vvadd(a,b):
    if type(b) in (int,float):
        return [i+b for i in a]
    return [i+j for i,j in zip(a,b)]

def vvsub(a,b):
    if type(b) in (int,float):
        return [i-b for i in a]
    return [i-j for i,j in zip(a,b)]

# Mean of deviations from initial value
def meand(v):
    return sum([i-v[0] for i in v])/len(v)

# Sum of squares/crossproducts of deviations 算方差或协方差
def ssd(u,v):
    return sum([(i-u[0])*(j-v[0]) for i,j in zip(u,v)])/(len(u)-1)


## MIJN EIGEN ROUTINE ##

# Quite short piece of code for diagonalizing symmetric 3x3 matrices :)
# 用于对称3x3矩阵对角化的简短代码
# Analytic solution for third order polynomial
# 三阶多项式的解析解
def solve_p3( a, b, c ):
    Q,R,a3 = (3*b-a**2)/9.0, (-27*c+a*(9*b-2*a**2))/54.0, a/3.0
    if Q**3 + R**2:
        t,R13 = math.acos(R/math.sqrt(-Q**3))/3, 2*math.sqrt(-Q)
        u,v,w = math.cos(t), math.sin(t+math.pi/6), math.cos(t+math.pi/3)
        return R13*u-a3, -R13*v-a3, -R13*w-a3
    else:
        R13   = math.sqrt3(R)
        return 2*R13-a3, -R13-a3, -R13-a3

# Normalization of 3-vector
def normalize(a):
    f = 1.0/math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
    return f*a[0],f*a[1],f*a[2]

# Eigenvectors for a symmetric 3x3 matrix:
# For symmetric matrix A the eigenvector v with root r satisfies
#   v.Aw = Av.w = rv.w = v.rw
#   v.(A-rI)w = v.Aw - v.rw = 0 for all w
# This means that for any two vectors p,q the eigenvector v follows from:
#   (A-rI)p x (A-rI)q
# The input is var(x),var(y),var(z),cov(x,y),cov(x,z),cov(y,z)
# The routine has been checked and yields proper eigenvalues/-vectors
def mijn_eigen_sym_3x3(a,d,f,b,c,e):
    a,d,f,b,c,e=1,d/a,f/a,b/a,c/a,e/a
    b2, c2, e2, df = b*b, c*c, e*e, d*f
    roots = list(solve_p3(-a-d-f, df-b2-c2-e2+a*(f+d), a*e2+d*c2+f*b2-a*df-2*b*c*e))
    roots.sort(reverse=True)
    ux, uy, uz = b*e-c*d, b*c-a*e, a*d-b*b
    u = (ux+roots[0]*c,uy+roots[0]*e,uz+roots[0]*(roots[0]-a-d))
    v = (ux+roots[1]*c,uy+roots[1]*e,uz+roots[1]*(roots[1]-a-d))
    w = u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0] # Cross product
    return normalize(u),normalize(v),normalize(w),roots

