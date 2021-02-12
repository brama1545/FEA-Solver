import math
import numpy as np

class Node:
    def __init__(self, idNo, position, dof):
        self.id = idNo
        self.position = position
        self.dof = dof
        self.disps = []
        self.Forces = []
        for dimension in range(dof):
            self.disps.append(Pointer())
            self.Forces.append(0)
        self.skew = []

    def setBC(self, disp, dimension):
        self.disps[dimension].value = disp

    def getSkewTransform(self):
        if len(self.skew) == 1:
            cskew = math.cos(self.skew[0])
            sskew = math.sin(self.skew[0])
            return np.array([[cskew, -sskew], [-sskew, cskew]])
        else:
            return np.identity(self.dof)

    def applyF(self, F, dimension):
        self.Forces[dimension] = F

    def getDisps(self):
        disps = []
        for dimension in range(self.dof):
            disps.append(self.disps[dimension].value)
        disps = np.matmul(self.getSkewTransform().transpose(), disps)
        return disps

    def __repr__(self):
        return '%s (%s)' % (self.id, self.position)


class Pointer:
    def __init__(self):
        self.value = None
