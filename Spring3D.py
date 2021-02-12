import numpy as np
import math


class Edge:
    def __init__(self, node1, node2):
        self.node1 = node1
        self.node2 = node2
        self.xdel = self.node1.position[0] - self.node2.position[0]
        self.ydel = self.node1.position[1] - self.node2.position[1]
        self.zdel = self.node1.position[2] - self.node2.position[2]
        self.length = self.getDist()
        self.area = None
        self.E = None
        self.stiffness = None
        self.localK = []
        cX = self.xdel/self.length
        cY = self.ydel/self.length
        cZ = self.zdel/self.length
        cX2 = np.square(cX)
        cY2 = np.square(cY)
        cZ2 = np.square(cZ)
        lam = np.array([[cX2, cX*cY, cX*cZ], [cX*cY, cY2, cY*cZ], [cX*cZ, cY*cZ, cZ2]])
        neglam = lam*-1
        x = np.concatenate((lam, neglam))
        y = np.concatenate((neglam, lam))
        self.transform = np.concatenate((x, y), 1)

        self.Tstar = [-cX, -cY, -cZ, cX, cY, cZ]

    def getDist(self):
        xterm = np.square(self.xdel)
        yterm = np.square(self.ydel)
        zterm = np.square(self.zdel)
        return np.sqrt(xterm + yterm + zterm)

    def getlocalK(self):
        r1 = np.array([1, 0, 0 -1, 0, 0])
        r2 = np.zeros(6)
        self.localK = self.stiffness * np.array([r1, r2, r2, -1 * r1, r2, r2])  # local stiffness matrix
        return self.localK

    def getGlobalK(self):
        intermediate = np.matmul(np.linalg.inv(self.Tstar), self.getlocalK())
        return np.matmul(intermediate, self.Tstar)
        # return self.transform * self.stiffness

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
