import numpy as np
import math
import Beam3D
import Bar3DPlanar


class Edge:
    def __init__(self, node1, node2):
        self.node1 = node1
        self.node2 = node2
        self.xdel = self.node2.position[0] - self.node1.position[0]
        self.ydel = self.node2.position[1] - self.node1.position[1]
        self.theta = math.atan2(self.ydel, self.xdel)
        self.length = self.getDist()
        self.area = None
        self.moi = None
        self.E = None
        self.stiffness = None
        self.localK = []
        c = math.cos(self.theta)
        c2 = np.square(c)
        s = math.sin(self.theta)
        s2 = np.square(s)
        self.stress = 0
        self.psi = 0
        self.bar = Bar3DPlanar.Edge(node1, node2)
        self.beam = Beam3D.Edge(node1, node2)

        self.Tstar = [[c, s, 0, 0, 0, 0], [-s,c,0,0,0,0], [0,0,1,0,0,0], [0, 0, 0, c, s, 0], [0,0,0,-s,c,0], [0,0,0,0,0,1]]

    def setStiffness(self, area, moi, E, psi):
        self.bar.setStiffness(area, E)
        self.beam.setStiffness(moi, E, psi)

    def getDist(self):
        xterm = np.square(self.xdel)
        yterm = np.square(self.ydel)
        return np.sqrt(xterm + yterm)

    def getlocalK(self):
        return self.bar.getlocalK() + self.beam.getlocalK()

    def getGlobalK(self):
        intermediate = np.matmul(np.transpose(self.Tstar), self.getlocalK())
        return np.matmul(intermediate, self.Tstar)
        #return self.bar.getGlobalK() + self.beam.getGlobalK()

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
