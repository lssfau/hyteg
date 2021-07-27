import re

def elStr(j, el):
    return '%d 2 2 0 5 ' %(j+1) + ' '.join('%d' %(idx+1) for idx in el)


def nodeStr(j, node):
    return '%d ' %(j+1) + ' '.join('%f' %coord for coord in node)

def parseFile(filename):
    EDGES = re.compile(r'.*Edges:.*')
    FACES = re.compile(r'.*Faces:.*')
    VTXES = re.compile(r'.*Vertices:.*')

    EDGE = re.compile(r'.*\s(\d+) \|\s+\d+ \|\s+(\d+) \|\s+(\d+) \|\s+\d+ \|\s+\d+.*')
    FACE = re.compile(r'.*\s\d+ \|\s+\d+ \|\s+(\d+) \|\s+(\d+) \|\s+\d+.*')
    VTX = re.compile(r'.*\s\d+ \|\s+\d+ \|\s+\[(\d\.?\d*), (\d\.?\d*), (\d\.?\d*)\].*')

    vtxes = []
    edges = {}
    faces = []

    read = ''

    with open(filename, 'r') as f:
        for line in f:
            for hdr,dat in zip([VTXES,EDGES,FACES], ['vtx', 'edge', 'face']):
                m = hdr.match(line)
                if m:
                    read = dat

            if (read == 'vtx'):
                m = VTX.match(line)
                if m:
                    vtxes.append([float(m[i]) for i in range(1,4)])
                    # print(vtxes)

            if (read == 'edge'):
                m = EDGE.match(line)
                if m:
                    edges[int(m[1])] = [int(m[2]), int(m[3])]
                    # print(edges)

            if (read == 'face'):
                m = FACE.match(line)
                if m:
                    faces.append(set(edges[int(m[1])] + edges[int(m[2])]))
                    # print(faces)

    return vtxes,faces


def writeFile(filename, nodes, elements):
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


if __name__ == '__main__':
    nodes, els = parseFile('hytegoutput.txt')
    writeFile('tmp.msh', nodes, els)



