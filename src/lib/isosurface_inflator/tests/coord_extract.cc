#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <string>
#include <MeshFEM/MeshIO.hh>

using namespace std;

int main(int argc, char *argv[])
{
    string path(argv[1]);
    std::vector<MeshIO::IOVertex> vertices;
    std::vector<MeshIO::IOElement> elements;
    MeshIO::load(path, vertices, elements);

    cout << setprecision(20);

    for (const auto &v : vertices)
        cout << v[0] << "\t" << v[1] << "\t" << v[2] << endl;
    return 0;
}
