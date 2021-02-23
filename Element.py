import numpy as np
import math
from abc import ABC, abstractmethod


class Element(ABC):

    def __init__(self):
        self.localK = []
        self.i = 0
        self.j = 0
        self.nodes = []

    def __iter__(self):
        return self

    def __next__(self):
        self.j += 1
        if self.j == len(self):
            self.i += 1
            if self.i == len(self) - 1:
                raise StopIteration
            else:
                self.j = self.i + 1
                return self.nodes[self.i], self.nodes[self.j]
        else:
            return self.nodes[self.i], self.nodes[self.j]

    def __repr__(self):
        return self.nodes

    @abstractmethod
    def __len__(self):
        pass

    def __getitem__(self, key):
        return self.nodes[key].id

    @abstractmethod
    def getlocalK(self):
        pass

    @abstractmethod
    def getGlobalK(self):
        pass

    @abstractmethod
    def getIntForces(self):
        pass

    @abstractmethod
    def getGlobaldisp(self):
        pass

    @abstractmethod
    def getLocaldisp(self):
        pass

    @abstractmethod
    def getStress(self):
        pass


class TwoNodeElement(Element):

    @abstractmethod
    def getlocalK(self):
        pass

    def getGlobalK(self):
        intermediate = np.matmul(np.transpose(self.Tstar), self.getlocalK())
        return np.matmul(intermediate, self.Tstar)

    def getIntForces(self):
        d_prime = np.matmul(self.Tstar, self.getGlobaldisp())
        force = np.matmul(self.getlocalK(), d_prime)
        return force

    def getGlobaldisp(self):
        disp = np.concatenate((self.nodes[0].getDisps(), self.nodes[1].getDisps()))
        return disp

    def getLocaldisp(self):
        return np.matmul(self.Tstar, self.getGlobaldisp())

    @abstractmethod
    def getStress(self):
        pass

    def __init__(self, node1, node2):
        super().__init__()
        self.nodes = [node1, node2]
        self.Tstar = []
        self.stress = 0
        self.area = None
        self.E = None
        self.G = None
        self.J = None
        self.stiffness = None
        self.moi = None

    def __len__(self):
        return 2


class ThreeNodeElement(Element):

    def __init__(self, nodei, nodej, nodem):
        super().__init__()
        self.nodes = [nodei, nodej, nodem]
        self.elmArea = 0
        self.poisson = 0
        self.stress = 0
        self.E = 0
        self.t = 0
        self.Tstar = []

    @abstractmethod
    def getlocalK(self):
        pass

    def getGlobalK(self):
        print(self.getlocalK()/4*math.pow(10, -6))
        print(self.Tstar)
        intermediate = np.matmul(np.transpose(self.Tstar), self.getlocalK())
        return np.matmul(intermediate, self.Tstar)

    def getIntForces(self):
        d_prime = np.matmul(self.Tstar, self.getGlobaldisp())
        force = np.matmul(self.getlocalK(), d_prime)
        return force

    def getGlobaldisp(self):
        disp = np.concatenate((self.nodes[0].getDisps(), self.nodes[1].getDisps(), self.nodes[2].getDisps()))
        return disp

    def getLocaldisp(self):
        return np.matmul(self.Tstar, self.getGlobaldisp())

    @abstractmethod
    def getStress(self):
        pass

    def __len__(self):
        return 3
