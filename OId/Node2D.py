class Node:
    def __init__(self, idNo, x, y):
        self.id = idNo
        self.position = [x, y]
        self.x = x
        self.y = y
        self.dispX = Pointer()
        self.dispY = Pointer()
        self.disps = [self.dispX, self.dispY]
        self.Fx = Pointer()
        self.Fy = Pointer()
        self.Forces = [self.Fx, self.Fy]

    def setBCX(self, disp):
        self.dispX.value = disp

    def setBCY(self, disp):
        self.dispY.value = disp

    def applyFx(self, F):
        self.Fx.value = F

    def applyFy(self, F):
        self.Fy.value = F

    def applyF(self, F, dimension):
        self.Forces[dimension].value = F

    def getDispXPntr(self):
        return self.dispX

    def getDispYPntr(self):
        return self.dispY

    def getFxPntr(self):
        return self.Fx

    def getFyPntr(self):
        return self.Fy

    def getDispX(self):
        return self.dispX.value

    def getDispY(self):
        return self.dispY.value

    def getFx(self):
        return self.Fx.value

    def getFy(self):
        return self.Fy.value

    def __repr__(self):
        return '%s (%s, %s)' % (self.id, self.x, self.y)


class Pointer:
    def __init__(self):
        self.value = None
