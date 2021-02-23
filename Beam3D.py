import numpy as np
import math
import Element


class Edge(Element.TwoNodeElement):
    def __init__(self, node1, node2):
        super().__init__(node1, node2)
        self.xdel = node2.position[0] - node1.position[0]
        self.ydel = node2.position[1] - node1.position[1]
        self.theta = math.atan2(self.ydel, self.xdel)
        self.length = self.getDist()
        c = math.cos(self.theta)
        s = math.sin(self.theta)
        self.psi = 0

        self.Tstar = [[c, s, 0, 0, 0, 0], [-s,c,0,0,0,0], [0,0,1,0,0,0], [0, 0, 0, c, s, 0], [0,0,0,-s,c,0], [0,0,0,0,0,1]]

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
        r1 = np.array([0, 12, 6*l, 0, -12, 6*l])
        r2 = np.array([0, 6*l, (4+p)*l*l, 0, -6*l, (2-p)*l*l])
        r3 = np.array([0, -12, -6*l, 0, 12, -6*l])
        r4 = np.array([0, 6*l, (2-p)*l*l, 0, -6*l, (4+p)*l*l])
        rz = np.zeros(6)
        self.localK = self.stiffness/(1+p) * np.array([rz, r1, r2, rz, r3, r4])  # local stiffness matrix
        return self.localK

    def getStress(self):
        return 0