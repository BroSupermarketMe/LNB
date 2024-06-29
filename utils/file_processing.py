import math
from .abbreviation import *

# TODO:isPDBAtom
def isPDBAtom(l):
    '''

    :param l: pdb文件中的一行
    :return: bool :判断这一行是否是原子信息
    '''
    return l.startswith("ATOM") or l.startswith("HETATM")

# 处理 ATOM或HETATM行 （PDB文件中的原子信息）
def pdbAtom(a):
    '''

    :param a: str: pdb文件中是原子信息的一行
    :return:
    '''
    ##01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ##ATOM   2155 HH11 ARG C 203     116.140  48.800   6.280  1.00  0.00
    ## ===>   atom name,   res name,     res id, chain,       x,            y,             z
    return (S(a[12:16]),S(a[17:20]),I(a[22:26]),a[21],F(a[30:38])/10,F(a[38:46])/10,F(a[46:54])/10)

d2r = 3.14159265358979323846264338327950288/180

# TODO:pdbBoxRead
# 处理CRYST1行（PDB文件中的晶格信息）
def pdbBoxRead(a):
    '''

    :param a: 'CRYST1  201.225  347.690  255.114  90.00  90.00  90.00 C 2 2 2      24 '
    :return:
    '''
    # Convert a PDB CRYST1 entry to a lattice definition.
    # Convert from Angstrom to nanometer 从埃转换为纳米
    fa, fb, fc, aa, ab, ac = [float(i) for i in a.split()[1:7]]
    ca, cb, cg, sg         = math.cos(d2r*aa), math.cos(d2r*ab), math.cos(d2r*ac) , math.sin(d2r*ac)
    wx, wy                 = 0.1*fc*cb, 0.1*fc*(ca-cb*cg)/sg
    wz                     = math.sqrt(0.01*fc*fc - wx*wx - wy*wy)
    return [0.1*fa, 0, 0, 0.1*fb*cg, 0.1*fb*sg, 0, wx, wy, wz]


def groAtom(a):
    #012345678901234567890123456789012345678901234567890
    #    1PRN      N    1   4.168  11.132   5.291
    ## ===>   atom name,   res name,     res id, chain,       x,          y,          z
    return (S(a[10:15]), S(a[5:10]),   I(a[:5]), " ", F(a[20:28]),F(a[28:36]),F(a[36:44]))

# TODO:groBoxRead
def groBoxRead(a):
    b = [F(i) for i in a.split()] + 6*[0] # Padding for rectangular boxes
    return b[0],b[3],b[4],b[5],b[1],b[6],b[7],b[8],b[2]



# 这个函数没有用过
def readBox(a):
    x = [ float(i) for i in a.split(",") ] + 6*[0]
    if len(x) == 12: # PDB format
        return pdbBoxRead("CRYST1 "+" ".join([str(i) for i in x]))
    else:            # GRO format
        return x[0],x[3],x[4],x[5],x[1],x[6],x[7],x[8],x[2]