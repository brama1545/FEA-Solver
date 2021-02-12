import numpy as np


class Edge:
    def __init__(self, node1, node2):
        self.node1 = node1
        self.node2 = node2
        self.stiffness = None

    def getlocalK(self):
        localK = self.stiffness * np.array([[1, -1], [-1, 1]])  # local stiffness matrix
        return localK

    def getGlobalK(self):
        return self.getlocalK()

    def getGlobalDisp(self):
        disp = np.concatenate((self.node1.getDisps, self.node2.getDisps))
        return disp

    def getStress(self):
        return 0

    def getIntForces(self):
        d_prime = self.getGlobalDisp()
        force = np.matmul(self.getlocalK(), d_prime)
        return [self.node1.id, self.node2.id], force
