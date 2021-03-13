import numpy as np
import math
import Element


class Edge(Element.ThreeNodeElement):
    def getlocalK(self):
        intermediate = np.matmul(np.transpose(self.B), self.D)
        self.localK = self.t * self.elmArea * np.matmul(intermediate, self.B)
        return self.localK

    def setProperties(self, thickness, E, poisson):
        self.t = thickness
        self.poisson = poisson
        self.E = E
        self.D = np.array([[1, poisson, 0], [poisson, 1, 0], [0, 0, (1-poisson)/2]])
        self.D *= E/(1-math.pow(poisson, 2))

    def setThermalProperties(self, delT, alpha):
        self.delT = delT
        self.CTE = alpha
        self.Etherm = np.array([self.CTE*self.delT, self.CTE*self.delT, 0])

    def getElmArea(self):
        positions = np.array([self.nodes[0].position, self.nodes[1].position, self.nodes[2].position])
        positions = np.concatenate((positions, [[1], [1], [1]]), 1)
        return .5 * np.linalg.det(positions)

    def getStress(self):
        intermediate = np.matmul(self.D, self.B)
        self.stress = np.matmul(intermediate, self.getGlobaldisp())
        return self.stress

    def getPrincipalStress(self):
        self.getStress()
        s1 = self.stress[0]
        s2 = self.stress[1]
        tau = self.stress[2]
        principal = [0, 0]
        principal[0] = (s1 + s2)/2 + math.sqrt(math.pow((s1 - s2)/2, 2) + math.pow(tau, 2))
        principal[1] = (s1 + s2)/2 - math.sqrt(math.pow((s1 - s2)/2, 2) + math.pow(tau, 2))
        return np.array(principal)

    def __init__(self, nodei, nodej, nodem):
        super().__init__(nodei, nodej, nodem)
        self.elmArea = self.getElmArea()
        self.xdel = nodej.position[0] - nodei.position[0]
        self.ydel = nodej.position[1] - nodei.position[1]
        self.theta = math.atan2(self.ydel, self.xdel)
        self.Ftherm = []
        self.CTE = 0
        self.delT = 0
        c = math.cos(self.theta)
        s = math.sin(self.theta)

        xi, yi = nodei.position[0], nodei.position[1]
        xj, yj = nodej.position[0], nodej.position[1]
        xm, ym = nodem.position[0], nodem.position[1]
        alpha_i = xj*ym - yj*xm
        alpha_j = yi*xm - xi*ym
        alpha_m = xi*yj - yi*xj
        beta_i = yj - ym
        beta_j = ym - yi
        beta_m = yi - yj
        gamma_i = xm - xj
        gamma_j = xi - xm
        gamma_m = xj - xi

        B_i = 1/2/self.elmArea * np.array([[beta_i, 0], [0, gamma_i], [gamma_i, beta_i]])
        B_j = 1/2/self.elmArea * np.array([[beta_j, 0], [0, gamma_j], [gamma_j, beta_j]])
        B_m = 1/2/self.elmArea * np.array([[beta_m, 0], [0, gamma_m], [gamma_m, beta_m]])

        self.B = np.concatenate((B_i, B_j, B_m), 1)
        T = np.array([[c, s], [-s, c]])
        zeros = np.zeros((2, 2))
        r1 = np.concatenate((T, zeros, zeros), 1)
        r2 = np.concatenate((zeros, T, zeros), 1)
        r3 = np.concatenate((zeros, zeros, T), 1)

        self.Tstar = np.concatenate((r1, r2, r3))
