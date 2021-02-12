import numpy as np
import fileparserMaster as fp

inFile = "InputMASTER.txt"
sigFigs = 8
dof = 0

workEq = []


def main():
    global dof
    global workEq
    np.set_printoptions(precision=sigFigs, suppress=True, linewidth=500)
    (nodeList, edgeList, dof) = fp.fileParse(inFile)  # parses input file
    numberNodes = len(nodeList)
    if len(workEq) == 0:
        workEq = np.zeros(numberNodes*dof)
    K = kludger(numberNodes, edgeList, nodeList)  # kludges global stiffness matrix together
    print('Global K matrix')
    print(K)
    (Disp, DispPntrs, Forces) = getSolnSpace(nodeList)  # Sets up f and u matrices depending on applied loads and BCs
    solvedDisps = dispSolver(K, Disp, DispPntrs, Forces)  # Solves for unknown displacements
    (allForces, allDisps) = forceSolver(K, nodeList)  # Multiplies K*displacements to find Force vector

    print('\nDisplacements (NodeID, X, Y...)')
    print(formatOutput(allDisps))
    print('\nForces (NodeID, X, Y...)')
    print(formatOutput(allForces))

    for edge in edgeList:
        print("\nLocal K-Matrix for element connecting nodes %s, %s" % (edge.node1.id, edge.node2.id))
        print(edge.getGlobalK())

    for edge in edgeList:
        print("\nLocal Elemental Forces for element connecting nodes %s, %s" % (edge.node1.id, edge.node2.id))
        print(np.matmul(edge.getlocalK(), edge.getLocalDisp()))
    # for edge in edgeList:
    #     print("\nAxial Stress for element connecting nodes %s, %s" % (edge.node1.id, edge.node2.id))
    #     print("{:e}".format(edge.getStress()))
    #     if not (edge.E is None):
    #         print("Axial Strain for element connecting nodes %s, %s" % (edge.node1.id, edge.node2.id))
    #         print("{:e}".format(edge.getStress()/edge.E))


def kludger(numberNodes, edgeList, nodeList):
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

    skewTransform = np.identity(len(K))  # Takes into account skew boundary conditions
    for node in nodeList:
        if node.skew:
            index = dof * (node.id - 1)
            size = len(node.getSkewTransform())
            skewTransform[index:index+size, index:index+size] = node.getSkewTransform()
    K = np.matmul(K, skewTransform.transpose())
    K = np.matmul(skewTransform, K)

    return K


def getSolnSpace(nodeList):
    Disp = []  # Global Displacement Vector
    DispPntrs = []  # "Pointers" to each unknown displacement. Keeps track of which variables disappear when BCs are applied
    Forces = []  # Applied forces
    weq = np.array(workEq)
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
            Forces.append(appliedForces[dimension] + workEq[dimension])
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
    Forces = np.matmul(K, Disp) - np.array(workEq)

    i = 0
    for node in nodeList:
        for dimension in range(dof):
            node.applyF(Forces[i], dimension)
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
