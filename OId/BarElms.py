import numpy as np
import Node2D as Node
import Edge2D as Edge

inFile = "Input1.csv"
dof = 2
sigFigs = 3


def main():
    np.set_printoptions(precision=sigFigs, suppress=True, linewidth=500)
    (nodeList, edgeList) = fileParse()  # parses input file
    numberNodes = len(nodeList)
    K = kludger(numberNodes, edgeList)  # kludges global stiffness matrix together
    print(K)
    (Disp, DispPntrs, Forces) = getSolnSpace(nodeList)  # Sets up f and u matrices depending on applied loads and BCs
    solvedDisps = dispSolver(K, Disp, DispPntrs, Forces)  # Solves for unknown displacements
    (allForces, allDisps) = forceSolver(K, DispPntrs, nodeList, solvedDisps)  # Multiplies K*displacements to find Force vector

    print('Displacements (NodeID, X, Y...)')
    print(formatOutput(allDisps))
    print('Forces (NodeID, X, Y...)')
    print(formatOutput(allForces))


def fileParse():
    file = open(inFile)
    nodeList = [0]
    edgeList = []
    file.readline()  # consume 'Nodes' line

    line = file.readline().split(' ')
    while not (line[0] == 'Edges'):  # read lines until 'Edges' comes up
        idNo = int(line[0])
        x = float(line[1])
        y = float(line[2])
        nodeList.insert(int(line[0]), Node.Node(idNo, x, y))
        line = file.readline().split(' ')  # get next line for iteration

    line = file.readline().split(' ')
    while not (line[0] == 'BCX\n'):  # read lines until 'BCX' comes up
        node1 = nodeList[int(line[0])]
        node2 = nodeList[int(line[1])]
        edgeList.append(Edge.Edge(node1, node2))

        if len(line) == 3:
            spring = float(line[2])
            edgeList[-1].stiffness = spring  # add stiffness of spring
        else:
            area = float(line[2])
            E = float(line[3])
            edgeList[-1].setBarStiffness(area, E)  # add stiffness of edge
        line = file.readline().split(' ')  # get next line for processing

    line = file.readline().split(' ')
    while not (line[0] == 'BCY\n'):  # read lines until 'BCY' comes up
        node = nodeList[int(line[0])]
        node.setBCX(float(line[1]))  # establish X boundary conditions
        line = file.readline().split(' ')

    line = file.readline().split(' ')
    while not (line[0] == 'FX\n'):  # read lines until 'FX' comes up
        node = nodeList[int(line[0])]
        node.setBCY(float(line[1]))  # establish Y boundary conditions
        line = file.readline().split(' ')

    line = file.readline().split(' ')
    while not (line[0] == 'FY\n'):  # read lines until 'FY' comes up
        node = nodeList[int(line[0])]
        node.applyFx(float(line[1]))  # establish applied X loads
        line = file.readline().split(' ')

    line = file.readline().split(' ')
    while not (line[0] == ''):  # read lines until EOF
        node = nodeList[int(line[0])]
        node.applyFy(float(line[1]))  # establish applied Y loads
        line = file.readline().split(' ')

    file.close()
    nodeList.remove(0)
    return nodeList, edgeList


def kludger(numberNodes, edgeList):
    K = np.zeros((dof * numberNodes, dof * numberNodes))  # initial stiffness matrix
    for edge in edgeList:
        globalK = edge.getGlobalK()
        node1index = dof * (edge.node1.id - 1)
        node2index = dof * (edge.node2.id - 1)  # Trust the math here...each node takes up dof # of spaces in the matrix
        for i in range(dof):
            for j in range(dof):
                K[node1index + i, node1index + j] += globalK[i, j]
                K[node2index + i, node2index + j] += globalK[dof + i, dof + j]
                K[node1index + i, node2index + j] += globalK[i, dof + j]
                K[node2index + i, node1index + j] += globalK[dof + i, j]
    return K


def getSolnSpace(nodeList):
    Disp = []  # Global Displacement Vector
    DispPntrs = []  # "Pointers" to each unknown displacement. Keeps track of which variabes disappear when BCs are applied
    Forces = []  # Applied forces
    for node in nodeList:
        if not (node.getDispX() is None):  # If a BC has been applied to this node in the X direction
            DispPntrs.append(None)  # Then don't consider it an unknown
            Disp.append(node.getDispX())  # And set the displacement to that value
        else:
            DispPntrs.append(node.getDispXPntr())  # Else: add it to the list of unknown
            Disp.append(0)
        if not (node.getDispY() is None):  # Same thing as above but for Y
            DispPntrs.append(None)
            Disp.append(node.getDispY())
        else:
            DispPntrs.append(node.getDispYPntr())
            Disp.append(0)

        if not (node.getFx() is None):  # Add any forces applied to each node
            Forces.append(node.getFx())
        else:
            Forces.append(0)
        if not (node.getFy() is None):
            Forces.append(node.getFy())
        else:
            Forces.append(0)

    return Disp, DispPntrs, Forces


def dispSolver(K, Disp, DispPntrs, Forces):
    annihilator = np.matmul(K, Disp)  # Annihilator is K*BCs. For non-homogenous BCs, use this to adjust the force matrix to make problem solveable
    Forces = Forces - annihilator  # Alters Forces vector to account for BCs

    unknowns = []  # indexes of unknown displacements
    for i in range(len(DispPntrs)):
        if not (DispPntrs[i] is None):
            unknowns.append(i)

    nonsingularK = np.zeros((len(unknowns), len(unknowns)))
    appliedForces = np.zeros((len(unknowns)))
    for i in range(len(unknowns)):
        for j in range(len(unknowns)):
            nonsingularK[i, j] = K[unknowns[i], unknowns[j]]
        appliedForces[i] = Forces[unknowns[i]]

    print(nonsingularK)

    return np.linalg.solve(nonsingularK, appliedForces)  # Takes reduced K and F matrix and solves for unknown displacements


def forceSolver(K, DispPntrs, nodeList, solvedDisps):
    i = 0
    for pntr in DispPntrs:
        if not (pntr is None):
            pntr.value = solvedDisps[i]
            i += 1

    Disp = []
    for i in range(len(nodeList)):
        Disp.append(nodeList[i].getDispX())
        Disp.append(nodeList[i].getDispY())
    Forces = np.matmul(K, Disp)

    allDisps = []
    i = 0
    for node in nodeList:
        node.applyFx(Forces[i])
        i += 1
        node.applyFy(Forces[i])
        i += 1
        allDisps.append(node.getDispX())
        allDisps.append(node.getDispY())

    return Forces, allDisps


def formatOutput(arr):
    i = 0
    output = np.zeros((int(len(arr) / dof), dof + 1))
    while i < len(arr) / dof:
        line = arr[i * dof:i * dof + dof]
        output[i] = np.concatenate(([i + 1], line))
        i += 1
    return output


if __name__ == '__main__':
    main()
