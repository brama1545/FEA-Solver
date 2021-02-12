import Edge1D as Edge
import Node as Node
dof = 1


def fileParse(inFile):
    file = open(inFile)
    nodeList = [0]
    edgeList = []
    file.readline()  # consume 'Nodes' line

    line = file.readline().split(' ')
    while not (line[0] == 'Edges'):  # read lines until 'Edges' comes up
        idNo = int(line[0])
        position = []
        nodeList.insert(int(line[0]), Node.Node(idNo, position, dof))
        line = file.readline().split(' ')  # get next line for iteration

    line = file.readline().split(' ')
    while not (line[0] == 'BCX\n'):  # read lines until 'BCX' comes up
        node1 = nodeList[int(line[0])]
        node2 = nodeList[int(line[1])]
        edgeList.append(Edge.Edge(node1, node2))
        spring = float(line[2])
        edgeList[-1].stiffness = spring  # add stiffness of spring
        line = file.readline().split(' ')  # get next line for processing

    line = file.readline().split(' ')
    while not (line[0] == 'F\n'):  # read lines until 'F' comes up
        node = nodeList[int(line[0])]
        node.setBC(float(line[1]), 0)  # establish X boundary conditions
        line = file.readline().split(' ')

    line = file.readline().split(' ')
    while not (line[0] == ''):  # read lines until EOF
        node = nodeList[int(line[0])]
        for dimension in range(dof):
            node.applyF(float(line[dimension+1]), dimension)
        line = file.readline().split(' ')

    file.close()
    nodeList.remove(0)
    return nodeList, edgeList
