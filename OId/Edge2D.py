import numpy as np
import math

class Edge:
    def __init__(self, node1, node2):
        self.node1 = node1
        self.node2 = node2
        self.length = self.getDist()
        self.theta = math.atan2(self.node1.y - self.node2.y, self.node1.x - self.node2.x)
        self.area = None
        self.E = None
        self.stiffness = None
        self.localK = []
        c = math.cos(self.theta)
        c2 = np.square(c)
        s = math.sin(self.theta)
        s2 = np.square(s)
        lam = np.array([[c2, c*s], [c*s, s2]])
        neglam = lam*-1
        x = np.concatenate((lam, neglam))
        y = np.concatenate((neglam, lam))
        self.transform = np.concatenate((x, y), 1)
        self.stress = None

        self.Tstar = [[c, s, 0, 0], [0, 0, c, s]]

    def setBarStiffness(self, area, E):
        self.area = area
        self.E = E
        self.stiffness = area*E/self.length

    def getDist(self):
        xterm = np.square(self.node1.x - self.node2.x)
        yterm = np.square(self.node1.y - self.node2.y)
        return np.sqrt(xterm + yterm)

    def getlocalK(self):
        self.localK = np.array(self.stiffness*[[1, -1], [-1, 1]])  # local stiffness matrix
        return self.localK

    def getGlobalK(self):
        return self.transform * self.stiffness

    def getGlobalDisp(self):
        disp = [self.node1.getDispX(), self.node1.getDispY(), self.node2.getDispX(), self.node2.getDispY()]
        return disp

    def getStress(self):
        c = np.matmul([-1, 1], self.Tstar) * self.E/self.area
        self.stress = np.matmul(c, self.getGlobalDisp())
        return self.stress

    def getIntForces(self):
        d_prime = np.matmul(self.Tstar, self.getGlobalDisp())
        force = np.matmul(self.getlocalK(), d_prime)
        return [self.node1.id, self.node2.id], force
