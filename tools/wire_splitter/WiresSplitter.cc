#include <iostream>
#include <fstream>
#include <vector>
#include <regex>
#include <sstream>
#include <cassert>

using namespace std;

int main(int argc, char ** argv) 
{
    string wirePath;
    string outputPath;
    string line;

    vector< vector<double> > vertices; 
    vector< vector<int> > edges; 
        

    if (argc == 3)
    {
        wirePath = argv[1]; 
        outputPath = argv[2];       
    } 
    else
    {
        cout << "Please, provide input file and output file (.wire, .obj extension)" << endl;
        return 0;
    }

    ifstream wireFile(wirePath);
    if (wireFile.is_open()) 
    {   
        while (getline(wireFile, line))
        {
            if (regex_match(line, regex("v [- ]*[0-9][0-9]*[\\.]*[0-9]* [- ]*[0-9][0-9]*[\\.]*[0-9]* [- ]*[0-9][0-9]*[\\.]*[0-9]*[ ]*$", regex_constants::grep)))
            {
                cout << "Found vertex" << endl;
                
                double x, y, z;
                char c;
                stringstream lineStream(line);
                lineStream >> c >> x >> y >> z;

                vector<double> point(3);
                point[0] = x;
                point[1] = y;
                point[2] = z;

                vertices.push_back(point);
            }
            else if (regex_match(line, regex("l [ ]*[0-9][0-9]* [ ]*[0-9][0-9]*[ ]*", regex_constants::grep)))
            {
                cout << "Found edge" << endl;

                long p1, p2;
                char c;
                stringstream lineStream(line);
                lineStream >> c >> p1 >> p2;

                vector<int> edge(2); 
                edge[0] = p1 - 1;
                edge[1] = p2 - 1;

                edges.push_back(edge);
            }
            else if (line.find_first_not_of("\t\n\v\f\r") != string::npos)
            {
                cout << "Warning: Line containing something different than vertex or edge" << endl;
            }
            cout << line << '\n';
        }
        wireFile.close(); 
    }

    cout << "Number of vertices: " << vertices.size() << endl;
    cout << "Number of edges: " << edges.size() << endl;

    // Now, let's create new graph
    vector< vector<double> > resultingVertices(vertices);
    vector< vector<int> > resultingEdges;

    for (vector< vector<int> >::iterator it = edges.begin() ; it != edges.end(); ++it)
    {
        vector<int> edge = *it;
        assert(edge.size() == 2);
        
        vector<double> v1 = vertices[edge[0]];
        vector<double> v2 = vertices[edge[1]];
        assert(v1.size() == 3);
        assert(v2.size() == 3);

        vector<double> newPoint(3);
        newPoint[0] = (v1[0] + v2[0]) / 2;
        newPoint[1] = (v1[1] + v2[1]) / 2;
        newPoint[2] = (v1[2] + v2[2]) / 2;

        int newVertexIndex = resultingVertices.size();

        vector<int> newEdge1(2);
        vector<int> newEdge2(2);
        newEdge1[0] = edge[0];
        newEdge1[1] = newVertexIndex;
        
        newEdge2[0] = edge[1];
        newEdge2[1] = newVertexIndex;

        resultingVertices.push_back(newPoint);
        resultingEdges.push_back(newEdge1);
        resultingEdges.push_back(newEdge2);
    }

    assert(resultingVertices.size() == (vertices.size() + edges.size()));

    ofstream outputFile(outputPath);

    for (vector< vector<double> >::iterator it = resultingVertices.begin() ; it != resultingVertices.end(); ++it)
    {
        vector<double> v = *it;
        outputFile << "v " << " " << v[0] << " " << v[1] << " " <<  v[2] << endl;
    }
    outputFile << endl;
    for (vector< vector<int> >::iterator it = resultingEdges.begin() ; it != resultingEdges.end(); ++it)
    {
        vector<int> l = *it;
        outputFile << "l " << " " << (l[0] + 1) << " " << (l[1] + 1)<< endl;
    }

    return 0;
}
