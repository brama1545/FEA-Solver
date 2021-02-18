import numpy as np
import math
import Beam3D
import Bar3DPlanar
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

    def getStress(self):
        if not(self.E is None):
            c = np.matmul([-1, 1], self.Tstar) * self.E/self.length
            self.stress = np.matmul(c, self.getGlobaldisp())
        return self.stress
