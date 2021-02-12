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

        self.Tstar = [[c, s, 0, 0], [0,0,0,0], [0, 0, c, s], [0,0,0,0]]
        # self.Tstar = [[c, s, 0, 0], [0, 0, c, s]]

    def getDist(self):
        xterm = np.square(self.xdel)
        yterm = np.square(self.ydel)
        return np.sqrt(xterm + yterm)

    def getlocalK(self):
        r1 = np.array([1, 0, -1, 0])
        r2 = np.zeros(4)
        self.localK = self.stiffness * np.array([r1, r2, -1 * r1, r2])  # local stiffness matrix
        # self.localK = self.stiffness * np.array([[1, -1], [-1, 1]])  # local stiffness matrix
        return self.localK

    def getGlobalK(self):
        return self.transform * self.stiffness

    def getLocalDisp(self):
        return np.matmul(self.Tstar, self.getGlobalDisp())

    def getGlobalDisp(self):
        disp = np.concatenate((self.node1.getDisps(), self.node2.getDisps()))
        return disp

    def getStress(self):
        return 0

    def getIntForces(self):
        d_prime = np.matmul(self.Tstar, self.getGlobalDisp())
        force = np.matmul(self.getlocalK(), d_prime)
        return [self.node1.id, self.node2.id], force
