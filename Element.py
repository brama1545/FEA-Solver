import numpy as np
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

    @abstractmethod
    def __len__(self):
        pass

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
        return [self.nodes[0].id, self.nodes[1].id], force

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
        self.stress = None
        self.area = None
        self.E = None
        self.G = None
        self.J = None
        self.stiffness = None
        self.moi = None

    def __len__(self):
        return 2


class ThreeNodeElement(Element):

    def getlocalK(self):
        pass

    def getGlobalK(self):
        pass

    def getIntForces(self):
        pass

    def getGlobaldisp(self):
        pass

    def getLocaldisp(self):
        pass

    def getStress(self):
        pass

    def __len__(self):
        return 3
