class WireReader(object):
    def __init__(self, wire_file):
        self.__wire_file = wire_file;
        self.__initialize();
        self.__parse();
        self.__remove_duplicated_edges();

    def __initialize(self):
        self.__dim = 0;
        self.vertices = [];
        self.edges = [];

    def __parse(self):
        with open(self.__wire_file, 'r') as fin:
            for line in fin:
                if len(line) == 0:
                    continue;
                if line[0] == 'v':
                    self.__parse_vertex(line);
                elif line[0] == 'l':
                    self.__parse_edge(line);
                elif line[0] == 'f':
                    self.__parse_face(line);
                else:
                    pass;

    def __parse_vertex(self, line):
        fields = line.split();
        assert(len(fields) == 3 or len(fields) == 4);
        self.dim = len(fields) - 1;
        v = [float(entry) for entry in fields[1:]];
        self.vertices.append(v);

    def __parse_edge(self, line):
        fields = line.split();
        assert(len(fields) == 3);
        e1 = int(fields[1]) - 1;
        e2 = int(fields[2]) - 1;
        self.edges.append([e1, e2]);

    def __parse_face(self, line):
        fields = line.split()[1:];
        num_vts = len(fields);
        for i in range(num_vts):
            e1 = int(fields[i]) -1;
            e2 = int(fields[(i+1) % num_vts]) -1;
            self.edges.append([e1, e2]);

    def __remove_duplicated_edges(self):
        edges = set([(max(e), min(e)) for e in self.edges]);
        self.edges = list(edges);

    @property
    def dim(self):
        return self.__dim;

    @dim.setter
    def dim(self, value):
        """ Set dim variable.  This can be done only once.
        """
        if self.__dim == 0:
            self.__dim = value;
        elif value != self.__dim:
            raise RuntimeError(
                    "Dim has already been setted to {}, not {}".format(
                        self.__dim, value));

