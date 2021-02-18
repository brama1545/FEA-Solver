import numpy as np


class Edge:
    def __init__(self, node1, node2):
        self.nodes = [node1, node2]
        self.stiffness = None

    def getlocalK(self):
        localK = self.stiffness * np.array([[1, -1], [-1, 1]])  # local stiffness matrix
        return localK

    def getGlobalK(self):
        return self.getlocalK()

    def getGlobalDisp(self):
        disp = np.concatenate((self.nodes[0].getDisps, self.nodes[1].getDisps))
        return disp

    def getStress(self):
        return 0

    def getIntForces(self):
        d_prime = self.getGlobalDisp()
        force = np.matmul(self.getlocalK(), d_prime)
        return [self.nodes[0].id, self.nodes[1].id], force
