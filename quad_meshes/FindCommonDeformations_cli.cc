////////////////////////////////////////////////////////////////////////////////
////  FindCommonDeformations_cli.cc
//////////////////////////////////////////////////////////////////////////////////
///*! @file
////      This reads in a set of .txt files containing deformations and find the union of them along with a map for each file. 
//*/ 
////  Author  : Morteza H Siboni (mhs), m.hakimi.siboni.@gmail.com
////  Company : New York University
////  Created : 02/19/2016
//////////////////////////////////////////////////////////////////////////////////
#include <boost/algorithm/string.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <algorithm>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "table.hh"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace fs;
using namespace std;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: FindCommonDeformations [options] list" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}


po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("list", po::value<string>(), "'comma-separated list of input deformation files'");

    po::positional_options_description p;
    p.add("list", 1);

    po::options_description visible_opts;
   	visible_opts.add_options()("help",			"Produce this help message")
        ("output,o", 	      po::value<string>()->default_value("output.txt"),	"output.txt file containing the common deformations")
        ("maps,m", 		      po::value<string>()->default_value("Map_"),		"prefix for map files")
        ;
	
	po::options_description cli_opts;
	cli_opts.add(visible_opts).add(hidden_opts);

	po::variables_map vm;

	try {
		po::store(po::command_line_parser(argc, argv).
				  options(cli_opts).positional(p).run(), vm);
		po::notify(vm);
	}
	catch (std::exception &e){
		cout << "Error: " << e.what() << endl << endl;
		usage(1, visible_opts);
	}

	bool fail = false;
	if (vm.count("list") == 0) {
		cout << "Error: must specify the comma-separated list of filenames!" << endl;
		fail = true;
	}

	if (fail || vm.count("help"))
		usage(fail, visible_opts);

	return vm;
}



void readTableFromFile(vector<vector<double>> & dataTable, string fileName)
{
	ifstream in;
	in.open(fileName);
	string line;
	
	while (getline(in, line)){
		
		vector<string> numbers;
		boost::trim(line);
		boost::split(numbers, line, boost::is_any_of("\t "),
					 boost::token_compress_on);

		vector<double> row;
		for (size_t i = 0; i < numbers.size(); ++i)
			row.push_back(stod(numbers[i]));
		
		dataTable.push_back(row);
	}
	in.close();
}

void dumpTable(const vector<vector<double>> & table, ostream & os)
{

	for (auto it = table.begin(); it != table.end(); ++it)
	{
		vector<double> row = *it;
		for (size_t i = 0; i < row.size() - 1; ++i)
			os << showpos << scientific <<  row[i] << "\t";
		os << showpos << scientific <<  row[row.size() - 1] << "\n";
	}
}


void dumpMap(const map<int, int> & inMap, ostream & os)
{
	for (auto & p : inMap)
		os << p.first << "\t" << p.second << endl;
}

vector<vector<double>> checkMap(map<int, int> inMap, vector<vector<double>> longTable, vector<vector<double>> shortTable)
{
	vector<vector<double>> outTable;
	outTable = longTable;
	for (auto & p : inMap)
	{
		int idInLong  = p.first;
		int idInShort = p.second;

		for (size_t i = 0; i < outTable[idInLong].size(); ++i)
			outTable[idInLong][i] = outTable[idInLong][i] - shortTable[idInShort][i]; 
	}
	return outTable;
}

////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on success)
*///////////////////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[])
{
	// read in the cmd aurguments
    po::variables_map args = parseCmdLine(argc, argv);


	vector<string> files;
	string fileList = args["list"].as<string>();

	boost::trim(fileList);
	boost::split(files, fileList, boost::is_any_of(", "));


	vector<size_t> sizes;

	// read all the deformation files and concatenate them 
	vector<vector<double>> allDeformaitons;
	for (size_t i = 0; i < files.size(); ++i){
		cout << "processing file " << files[i] << endl;
		vector<vector<double>> currentDeformatons;
		readTableFromFile(currentDeformatons, files[i]);
		sizes.push_back(currentDeformatons.size());
		for (size_t j = 0; j < currentDeformatons.size(); ++j)
			allDeformaitons.push_back(currentDeformatons[j]);
	}

	// generate the table object with all the deformations
	table<double> deformationTable(allDeformaitons);

	// get the common deformation 
	vector<vector<double>> uniqueDeformaoitns = deformationTable.getUniqueTable();
	ofstream out;
	out.open(args["output"].as<string>());
	dumpTable(uniqueDeformaoitns, out);
	out.close();

	// get the global map
	map<int, int> conversionMap = deformationTable.getMap();

	// divide the map based on size of the input files
	vector<map<int, int>> seperatedMaps;
	map<int, int>::iterator it = conversionMap.begin();
	for (size_t i = 0; i < sizes.size(); ++i){
		map<int, int> currentMap;
		for (size_t j = 0; j < sizes[i]; ++j){
			currentMap.insert(pair<int, int>(j, it->second));
			it++;
		}
		seperatedMaps.push_back(currentMap);
	}

	// write down the maps
	for (size_t i = 0; i < seperatedMaps.size(); ++i){
		out.open(args["maps"].as<string>() + to_string(i) + ".txt");
		dumpMap(seperatedMaps[i], out);
		out.close();
	}

/* 	// check the correctness of the maps */
/* 	for (size_t i = 0; i < files.size(); ++i){ */
/* 		cout << "checking the " << i << "th map:" << endl; */
/* 		vector<vector<double>> actualDefs, generatedDefs; */
/* 		readTableFromFile(actualDefs, files[i]); */
/* 		for (auto & p : seperatedMaps[i]) */
/* 			generatedDefs.push_back(uniqueDeformaoitns[p.second]); */

/* 		for (size_t row = 0; row < actualDefs.size(); ++row){ */
/* 			for (size_t ele = 0; ele < actualDefs[row].size(); ++ele) */
/* 				cout << abs(actualDefs[row][ele] - generatedDefs[row][ele]) << "\t"; */
/* 			cout << endl; */
/* 		} */
/* 		cout << endl; */
/* 	} */

/* 	// compute the stretches/jacobians for each quad in the quadMesh */
/* 	getDeformations(pmesh, defs, args["mode"].as<string>()); */

/* 	// truncate the stretches/jacobians */
/* 	defs = truncate(args["truncate"].as<size_t>(), defs); */	

/* 	table<double, 2> defTable(defs); */



/* 	// getting/writing all the quad deformations to file */
/* 	vector<vector<double>> allDefs    = defTable.getTable(); */
/* 	ofstream out; */
/* 	out.open(args["allDefs"].as<string>()); */
/* 	dumpTable(allDefs, out); */
/* 	out.close(); */

/* 	// getting/writing unique quad deformations to file */
/* 	vector<vector<double>> uniqueDefs = defTable.getUniqueTable(); */
/* 	out.open(args["uniqueDefs"].as<string>()); */
/* 	dumpTable(uniqueDefs, out); */
/* 	out.close(); */

/* 	// getting/writing the conversion map to file */
/* 	map<int, int> conversionMap       = defTable.getMap(); */
/* 	out.open(args["conversionMap"].as<string>()); */
/* 	dumpMap(conversionMap, out); */
/* 	out.close(); */

    return 0;
}
