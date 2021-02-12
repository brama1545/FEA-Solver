import numpy as np
import math

class Edge:
    def __init__(self, node1, node2):
        self.node1 = node1
        self.node2 = node2
        self.xdel = self.node2.position[0] - self.node1.position[0]
        self.ydel = self.node2.position[1] - self.node1.position[1]
        self.theta = math.atan2(self.ydel, self.xdel)
        self.length = self.getDist()
        self.moi = None
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
        self.psi = 0

        self.Tstar = [[c, s, 0, 0], [0, 0, c, s]]

    def setStiffness(self, moi, E, psi):
        self.moi = moi
        self.E = E
        self.psi = psi
        l3 = math.pow(self.length, 3)
        self.stiffness = moi*E/l3

    def getDist(self):
        xterm = np.square(self.xdel)
        yterm = np.square(self.ydel)
        return np.sqrt(xterm + yterm)

    def getlocalK(self):
        l = self.length
        p = self.psi
        r1 = np.array([12, 6*l, -12, 6*l])
        r2 = np.array([6*l, (4+p)*l*l, -6*l, (2-p)*l*l])
        r3 = np.array([-12, -6*l, 12, -6*l])
        r4 = np.array([6*l, (2-p)*l*l, -6*l, (4+p)*l*l])
        self.localK = self.stiffness/(1+p) * np.array([r1, r2, r3, r4])  # local stiffness matrix
        return self.localK

    def getGlobalK(self):
        #print(self.transform)
        #return self.transform * self.stiffness
        return self.getlocalK()

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
