from sys import argv

def elStr(j, el):
    return '%d 2 2 0 5 ' %(j+1) + ' '.join('%d' %idx for idx in el)


def nodeStr(j, node):
    return '%d ' %(j+1) + ' '.join('%f' %coord for coord in node)


def generateInitialMesh():
    nodes = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
    elems = [[1, 2, 3], [4, 3, 2]]
    return nodes, elems


def readMeshFromFile(filename):
    nodes = []
    elems = []
    readnodes = False
    readelems = False

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('$Nodes'):
                readnodes = True
            elif line.startswith('$EndNodes'):
                readnodes = False
            elif line.startswith('$Elements'):
                readelems = True
            elif line.startswith('$EndElements'):
                readelems = False
            elif readnodes:
                node = [float(s) for s in line.split()]
                if len(node) == 4:
                    nodes.append(node[1:])
            elif readelems:
                elem = [int(s) for s in line.split()]
                if len(elem) == 8:
                    elems.append(elem[-3:])

    return nodes, elems


def writeMeshToFile(filename, nodes, elements):
    with open(filename, 'w') as f:
        f.writelines('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')

        f.writelines('$Nodes\n%d\n' %len(nodes))
        for j, node in enumerate(nodes):
            f.writelines('%s\n' %nodeStr(j, node))
        f.writelines('$EndNodes\n')

        f.writelines('$Elements\n%d\n' %len(elements))
        for j, el in enumerate(elements):
            f.writelines('%s\n' %elStr(j, el))
        f.writelines('$EndElements\n')


def refinement(nodes, elems, refine:list):
    newNodes = {}
    newNodeIdxs = {}

    # find edge midpoints to refine edges
    for el in refine:
        getNodesForRefinement(nodes, newNodes, elems[el])

    # add new nodes
    idx = len(nodes)
    for edge,midpt in newNodes.items():
        idx += 1
        newNodeIdxs[edge] = idx
        nodes.append(midpt)

    # add new elements
    for el in refine:
        elems += getSubElements(elems[el], newNodeIdxs)

    # remove refined elements
    refine.sort(reverse=True)
    for el in refine:
        elems.pop(el)

    return nodes, elems


def getEdges(el):
    return [tuple(sorted((el[j]-1, el[j-1]-1))) for j in range(3)]


def getNodesForRefinement(nodes, newNodes, el):
    for idxa,idxb in getEdges(el):
        a = nodes[idxa]
        b = nodes[idxb]
        newNodes[(idxa, idxb)] = [(a[k] + b[k])/2 for k in [0, 1]] + [0]


def getIntersection(edge1, edge2):
    for k in edge1:
        if k in edge2:
            return k+1


def getSubElements(el, midpts):
    elems = []
    edges = getEdges(el)
    for k in range(3):
        e1 = edges[k]
        e2 = edges[k-1]
        elems.append([getIntersection(e1, e2), midpts[e1], midpts[e2]])
    elems.append([midpts[e] for e in edges])

    return elems


if __name__ == '__main__':

    inputfile = argv[1]
    outputfile = argv[2]
    refine = []
    if len(argv) > 3:
        refine = [int(k) for k in argv[3].split()]

    if inputfile == '0':
        nodes,elems = generateInitialMesh()
    else:
        nodes, elems = readMeshFromFile(inputfile)

    nodes, elems = refinement(nodes, elems, refine)

    writeMeshToFile(outputfile, nodes, elems)