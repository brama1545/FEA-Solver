import numpy as np
import fileparserMaster as fp

inFile = "InputFINAL5.txt"
sigFigs = 6
dof = 0

auxillaryForce = []


def main():
    global dof
    global auxillaryForce
    np.set_printoptions(precision=5, suppress=True, linewidth=500)
    (nodeList, edgeList, dof) = fp.fileParse(inFile)  # parses input file
    numberNodes = len(nodeList)
    if len(auxillaryForce) == 0:
        auxillaryForce = np.zeros(numberNodes * dof)
    K = kludger(numberNodes, edgeList, nodeList)  # kludges global stiffness matrix together
    print('Global K matrix X 10^-6')
    print(K/1000000)
    np.set_printoptions(precision=sigFigs, suppress=True, linewidth=500)
    (Disp, DispPntrs, Forces) = getSolnSpace(nodeList, edgeList)  # Sets up f and u matrices depending on applied loads and BCs
    solvedDisps = dispSolver(K, Disp, DispPntrs, Forces)  # Solves for unknown displacements
    (allForces, allDisps) = forceSolver(K, nodeList)  # Multiplies K*displacements to find Force vector

    print('\nDisplacements (NodeID, X, Y...)')
    print(formatOutput(allDisps))
    print('\nForces (NodeID, X, Y...)')
    print(formatOutput(allForces))

    for edge in edgeList:
        print("\nLocal K-Matrix for element connecting nodes %s, %s, %s X 10^-6" % (edge[0], edge[1], edge[2]))
        print(edge.getlocalK()/100000)

    for edge in edgeList:
        print("\nLocal Elemental Forces for element connecting nodes %s, %s, %s" % (edge[0], edge[1], edge[2]))
        print(np.matmul(edge.getlocalK(), edge.getLocaldisp()))

    for edge in edgeList:
        print("\nLocal stress for element connecting nodes %s, %s, %s X in MPa" % (edge[0], edge[1], edge[2]))
        print(edge.getStress()/1000000)
        print("Principal stress for element connecting nodes %s, %s, %s in MPa" % (edge[0], edge[1], edge[2]))
        print(edge.getPrincipalStress()/1000000)


def kludger(numberNodes, edgeList, nodeList):
    K = np.zeros((dof * numberNodes, dof * numberNodes))  # initial stiffness matrix
    for edge in edgeList:
        globalK = edge.getGlobalK()
        for node1, node2 in edge:  # For each node pair in edge (see Element.py iterator)
            node1index = dof * (node1.id - 1)
            node2index = dof * (node2.id - 1)  # Trust the math here...each node takes up dof # of spaces in the matrix
            start_i = dof*edge.i
            end_i = start_i+dof
            start_j = dof*edge.j
            end_j = start_j + dof
            K[node1index:node1index+dof, node2index:node2index+dof] += globalK[start_i:end_i, start_j:end_j]
            K[node2index:node2index+dof, node1index:node1index+dof] += globalK[start_j:end_j, start_i:end_i]
            # for i in range(dof):
            #     for j in range(dof):
            #         #K[node1index + i, node1index + j] += globalK[i, j]
            #         #K[node2index + i, node2index + j] += globalK[dof*edge.i + i, dof*edge.j + j]
            #         K[node1index + i, node2index + j] += globalK[dof*edge.i + i, dof*edge.j + j]
            #         K[node2index + i, node1index + j] += globalK[dof*edge.i + i, dof*edge.j + j]

        for i in range(len(edge)):
            nodeindex = dof*(edge.nodes[i].id - 1)
            K[nodeindex:nodeindex+dof, nodeindex:nodeindex+dof] += globalK[dof*i:dof*(i+1), dof*i:dof*(i+1)]

    skewTransform = np.identity(len(K))  # Takes into account skew boundary conditions
    for node in nodeList:
        if node.skew:
            index = dof * (node.id - 1)
            size = len(node.getSkewTransform())
            skewTransform[index:index+size, index:index+size] = node.getSkewTransform()
    K = np.matmul(K, skewTransform.transpose())
    K = np.matmul(skewTransform, K)

    return K


def getSolnSpace(nodeList, edgeList):
    global auxillaryForce
    Disp = []  # Global Displacement Vector
    DispPntrs = []  # "Pointers" to each unknown displacement. Keeps track of which variables disappear when BCs are applied
    Forces = []  # Applied forces
    for node in nodeList:
        BCs = node.disps
        appliedForces = node.Forces
        for dimension in range(dof):
            if not (BCs[dimension].value is None):  # If a BC has been applied to this node in this direction
                DispPntrs.append(None)  # Then don't consider it an unknown
                Disp.append(BCs[dimension].value)  # And set the displacement to that value
            else:
                DispPntrs.append(BCs[dimension])  # Else: add it to the list of unknown
                Disp.append(0)
            Forces.append(appliedForces[dimension])
    for element in edgeList:
        Ftherm = element.getFtherm()
        i = 0
        for node in element.nodes:
            nodeid = node.id
            auxillaryForce[dof * (nodeid - 1):dof * nodeid] += Ftherm[dof*i: dof*(i+1)]
            i += 1

    Forces = Forces + auxillaryForce
    print('\nForce Matrix For Governing Eqn')
    print(Forces)

    return Disp, DispPntrs, Forces


def dispSolver(K, Disp, DispPntrs, Forces):
    annihilator = np.matmul(K, Disp)  # Annihilator is K*BCs. For non-homogenous BCs, this adjusts the force matrix to make problem solveable
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

    print('\nReduced K matrix with BCs applied')
    print(nonsingularK)

    solvedDisps = np.linalg.solve(nonsingularK, appliedForces)  # Takes reduced K and F matrix and solves for unknown displacements
    i = 0
    for pntr in DispPntrs:
        if not (pntr is None):  # After having solved the displacements, store them in each nodes variable
            pntr.value = solvedDisps[i]
            i += 1

    return solvedDisps


def forceSolver(K, nodeList):
    Disp = []
    for node in nodeList:
        for dimension in range(dof):
            Disp.append(node.getDisps()[dimension])
    Forces = np.matmul(K, Disp) - np.array(auxillaryForce)

    i = 0
    for node in nodeList:
        for dimension in range(dof):
            node.setF(Forces[i], dimension)
            i += 1

    return Forces, Disp


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
