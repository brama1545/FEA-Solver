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
        c2 = np.square(c)
        s = math.sin(self.theta)
        s2 = np.square(s)
        lam = np.array([[c2, c*s], [c*s, s2]])
        neglam = lam*-1
        x = np.concatenate((lam, neglam))
        y = np.concatenate((neglam, lam))
        self.transform = np.concatenate((x, y), 1)

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

    def getStress(self):
        cprime = self.E/self.length * np.matmul([-1, 1], self.Tstar)
        self.stress = np.matmul(self.getGlobaldisp(), cprime)
        return self.stress
