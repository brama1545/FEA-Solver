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
        self.G = None
        self.J = None
        self.bstiffness = None
        self.tstiffness = None
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

        self.Tstar = [[1, 0, 0, 0, 0, 0], [0,c,s,0,0,0], [0,-s,c,0,0,0], [0, 0, 0, 1, 0, 0], [0,0,0,0,c,s], [0,0,0,0,-s,c]]

    def setStiffness(self, J, moi, E, G, psi):
        self.moi = moi
        self.E = E
        self.G = G
        self.J = J
        self.psi = psi
        l3 = math.pow(self.length, 3)
        self.bstiffness = moi*E/l3
        self.tstiffness = G*J/self.length

    def getDist(self):
        xterm = np.square(self.xdel)
        yterm = np.square(self.ydel)
        return np.sqrt(xterm + yterm)

    def getlocalK(self):
        l = self.length
        p = self.psi
        r1 = np.array([12, 0, 6*l, -12, 0, 6*l])
        rz = np.zeros(6)
        r3 = np.array([6*l, 0, (4+p)*l*l, -6*l, 0, (2-p)*l*l])
        r4 = np.array([-12, 0, -6*l, 12, 0, -6*l])
        r6 = np.array([6*l, 0, (2-p)*l*l, -6*l, 0, (4+p)*l*l])
        bendingK = self.bstiffness/(1+p) * np.array([r1, rz, r3, r4, rz, r6])  # local stiffness matrix
        torsionRow = self.tstiffness * np.array([0, 1, 0, 0, -1, 0])
        torsionK = np.array([rz, torsionRow, rz, rz, -1 * torsionRow, rz])
        self.localK = torsionK + bendingK
        return self.localK

    def getGlobalK(self):
        intermediate = np.matmul(np.transpose(self.Tstar), self.getlocalK())
        return np.matmul(intermediate, self.Tstar)

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
