class WireWriter(object):
    def __init__(self, wire_file):
        self.__wire_file = wire_file;

    def write(self, vertices, edges):
        with open(self.__wire_file, 'w') as fout:
            for v in vertices:
                fout.write("v {v[0]} {v[1]} {v[2]}\n".format(v=v));
            for e in edges:
                fout.write("l {e[0]} {e[1]}\n".format(e=e+1));


