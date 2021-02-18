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
        self.stress = 0

        self.Tstar = [[c, s, 0, 0, 0, 0], [-s,c,0,0,0,0], [0,0,1,0,0,0], [0, 0, 0, c, s, 0], [0,0,0,-s,c,0], [0,0,0,0,0,1]]

    def getDist(self):
        xterm = np.square(self.xdel)
        yterm = np.square(self.ydel)
        return np.sqrt(xterm + yterm)

    def setStiffness(self, area, E):
        self.area = area
        self.E = E
        self.stiffness = area*E/self.length

    def getlocalK(self):
        r1 = np.array([1, 0, 0, -1, 0, 0])
        r2 = np.zeros(6)
        self.localK = self.stiffness * np.array([r1, r2, r2, -1 * r1, r2, r2])  # local stiffness matrix
        return self.localK

    def getStress(self):
        return self.stress
