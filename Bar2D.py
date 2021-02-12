import numpy as np
import math

class Edge:
    def __init__(self, node1, node2):
        self.node1 = node1
        self.node2 = node2
        self.xdel = self.node1.position[0] - self.node2.position[0]
        self.ydel = self.node1.position[1] - self.node2.position[1]
        self.theta = math.atan2(self.ydel, self.xdel)
        self.length = self.getDist()
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
        self.stress = 0

        self.Tstar = [[c, s, 0, 0], [0, 0, c, s]]

    def setStiffness(self, area, E):
        self.area = area
        self.E = E
        self.stiffness = area*E/self.length

    def getDist(self):
        xterm = np.square(self.xdel)
        yterm = np.square(self.ydel)
        return np.sqrt(xterm + yterm)

    def getlocalK(self):
        self.localK = self.stiffness * np.array([[1, -1], [-1, 1]])  # local stiffness matrix
        return self.localK

    def getGlobalK(self):
        return self.transform * self.stiffness

    def getLocalDisp(self):
        return np.matmul(self.Tstar, self.getGlobalDisp())

    def getGlobalDisp(self):
        disp = np.concatenate((self.node1.getDisps(), self.node2.getDisps()))
        return disp

    def getStress(self):
        if not(self.E is None):
            c = np.matmul([-1, 1], self.Tstar) * self.E/self.length
            self.stress = np.matmul(c, self.getGlobalDisp())
        return self.stress

    def getIntForces(self):
        d_prime = np.matmul(self.Tstar, self.getGlobalDisp())
        force = np.matmul(self.getlocalK(), d_prime)
        return [self.node1.id, self.node2.id], force
