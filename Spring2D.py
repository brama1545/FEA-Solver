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

        self.Tstar = [[c, s, 0, 0], [0,0,0,0], [0, 0, c, s], [0,0,0,0]]

    def getDist(self):
        xterm = np.square(self.xdel)
        yterm = np.square(self.ydel)
        return np.sqrt(xterm + yterm)

    def getlocalK(self):
        r1 = np.array([1, 0, -1, 0])
        r2 = np.zeros(4)
        self.localK = self.stiffness * np.array([r1, r2, -1 * r1, r2])  # local stiffness matrix
        return self.localK

    def getGlobalK(self):
        return self.transform * self.stiffness

    def getStress(self):
        return 0
