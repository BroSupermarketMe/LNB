from .file_processing import *
from .algebra import vvadd,vvsub

class Structure:
    def __init__(self,filename=None):
        self.title   = ""
        self.atoms   = []
        self.coord   = []
        self.rest    = []
        self.box     = []
        self._center = None

        if filename:
            lines = open(filename).readlines()
            # Try extracting PDB atom/hetatm definitions
            self.rest   = []
            # pdbAtom() returns [(atom name, res name, res id, chain, x, y, z),...]
            self.atoms  = [pdbAtom(i) for i in lines if isPDBAtom(i) or self.rest.append(i)]
            if self.atoms:
                # This must be a PDB file
                self.title = "THIS IS INSANE!\n"
                for i in self.rest:
                    if i.startswith("TITLE"):
                        self.title = i
                self.box   = [0,0,0,0,0,0,0,0,0]
                for i in self.rest:
                    if i.startswith("CRYST1"):
                        self.box = pdbBoxRead(i)
            else:
                # This should be a GRO file
                self.atoms = [groAtom(i) for i in lines[2:-1]]
                self.rest  = [lines[0],lines[1],lines[-1]]
                self.box   = groBoxRead(lines[-1])
                self.title = lines[0]
            self.coord = [i[4:7] for i in self.atoms]
            self.center() # Calculate the center of mass

    def __nonzero__(self):
        return bool(self.atoms)

    def __len__(self):
        return len(self.atoms)

    def __iadd__(self,s):
        for i in range(len(self)):
            self.coord[i] = vvadd(self.coord[i],s)
        return self

    def center(self,other=None):
        '''

        :param other: 给定other矢量算相对other的坐标，同时把self.coord中的坐标也改变
        :return:
        '''
        if not self._center:
            self._center = [ sum(i)/len(i) for i in zip(*self.coord)]
        if other:
            s = vvsub(other,self._center)
            for i in range(len(self)):
                self.coord[i] = vvadd(self.coord[i],s)
            self._center = other
            return s # return the shift
        return self._center

    def diam(self):
        if self._center != (0,0,0):
            self.center((0,0,0))
        return 2*math.sqrt(max([i*i+j*j+k*k for i,j,k in self.coord]))

    def diamxy(self):
        if self._center != (0,0,0):
            self.center((0,0,0))
        return 2*math.sqrt(max([i*i+j*j for i,j,k in self.coord]))

    def fun(self,fn):
        return [fn(i) for i in zip(*self.coord)]