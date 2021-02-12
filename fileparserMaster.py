import Spring1D
import Spring2D
import Spring3D
import Spring3DPlanar
import Bar2D
import Bar3D
import Bar3DPlanar
import Beam2D
import Beam3D
import Frame
import Grid
import Node


def fileParse(inFile):
    file = open(inFile)
    nodeList = [0]
    edgeList = []
    file.readline()  # consume 'dof' line
    dofs = file.readline().split(' ')
    kdof = int(dofs[0])
    gdof = int(dofs[1])
    file.readline()  # consume 'Nodes' line

    line = file.readline().split(' ')
    while not (line[0] == 'Springs'):  # read lines until 'Springs' comes up
        idNo = int(line[0])
        position = []
        for dimension in range(gdof):
            position.append(float(line[dimension+1]))
        nodeList.insert(int(line[0]), Node.Node(idNo, position, kdof))
        line = file.readline().split(' ')  # get next line for iteration

    line = file.readline().split(' ')
    while not (line[0] == 'Bars'):  # read lines until 'Bars' comes up
        node1 = nodeList[int(line[0])]
        node2 = nodeList[int(line[1])]
        if kdof == 1:
            edgeList.append(Spring1D.Edge(node1, node2))
        elif kdof == 2:
            edgeList.append(Spring2D.Edge(node1, node2))
        elif kdof == 3:
            if gdof == 2:
                edgeList.append(Spring3DPlanar.Edge(node1, node2))
            if gdof == 3:
                edgeList.append(Spring3D.Edge(node1, node2))

        spring = float(line[2])
        edgeList[-1].stiffness = spring  # add stiffness of spring
        line = file.readline().split(' ')  # get next line for processing

    line = file.readline().split(' ')
    while not (line[0] == 'Beams'):  # read lines until 'Beams' comes up
        node1 = nodeList[int(line[0])]
        node2 = nodeList[int(line[1])]
        if kdof == 2:
            edgeList.append(Bar2D.Edge(node1, node2))
        elif kdof == 3:
            if gdof == 3:
                edgeList.append(Bar3D.Edge(node1, node2))
            else:
                edgeList.append(Bar3DPlanar.Edge(node1, node2))

        area = float(line[2])
        E = float(line[3])
        edgeList[-1].setStiffness(area, E)  # add stiffness of edge
        line = file.readline().split(' ')  # get next line for processing

    line = file.readline().split(' ')
    while not (line[0] == 'Frames'):  # read lines until 'BCX' comes up
        node1 = nodeList[int(line[0])]
        node2 = nodeList[int(line[1])]
        if kdof == 2:
            edgeList.append(Beam2D.Edge(node1, node2))
        if kdof == 3:
            edgeList.append(Beam3D.Edge(node1, node2))

        moi = float(line[2])
        E = float(line[3])
        psi = 0
        if len(line) == 5:
            psi = line[4]
        edgeList[-1].setStiffness(moi, E, psi)  # add stiffness of edge
        line = file.readline().split(' ')  # get next line for processing

    line = file.readline().split(' ')
    while not (line[0] == 'Grid'):  # read lines until 'BCX' comes up
        node1 = nodeList[int(line[0])]
        node2 = nodeList[int(line[1])]
        edgeList.append(Frame.Edge(node1, node2))

        area = float(line[2])
        moi = float(line[3])
        E = float(line[4])
        psi = 0
        if len(line) == 6:
            psi = line[5]
        edgeList[-1].setStiffness(area, moi, E, psi)  # add stiffness of edge
        line = file.readline().split(' ')  # get next line for processing

    line = file.readline().split(' ')
    while not (line[0] == 'BC'):  # read lines until 'BCX' comes up
        node1 = nodeList[int(line[0])]
        node2 = nodeList[int(line[1])]
        edgeList.append(Grid.Edge(node1, node2))

        J = float(line[2])
        moi = float(line[3])
        E = float(line[4])
        G = float(line[5])
        psi = 0
        if len(line) == 7:
            psi = line[6]
        edgeList[-1].setStiffness(J, moi, E, G, psi)  # add stiffness of edge
        line = file.readline().split(' ')  # get next line for processing

    line = file.readline().split(' ')
    while not (line[0] == 'Skew'):  # read lines until 'BCY' comes up "NODE, DOF_INDEX, VALUE"
        node = nodeList[int(line[0])]
        node.setBC(float(line[2]), int(line[1])-1)  # establish X boundary conditions
        line = file.readline().split(' ')

    line = file.readline().split(' ')
    while not (line[0] == 'Forces'):
        node = nodeList[int(line[0])]
        for i in range(1, len(line)):
            node.skew.append(float(line[i]))

    line = file.readline().split(' ')
    while not (line[0] == ''):  # read lines until EOF
        node = nodeList[int(line[0])]
        for dimension in range(kdof):
            node.applyF(float(line[dimension + 1]), dimension)
        line = file.readline().split(' ')

    file.close()
    nodeList.remove(0)
    return nodeList, edgeList, kdof
