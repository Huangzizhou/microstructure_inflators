////////////////////////////////////////////////////////////////////////////////
////  UnMap_cli.cc
//////////////////////////////////////////////////////////////////////////////////
///*! @file
////      This reads in a map and a set of data and unmaps the data to an output 
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
    cout << "Usage: UnMap_cli [options] outputData.txt" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}


po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("output", po::value<string>(), "a text file containing the unmaped data");

    po::positional_options_description p;
    p.add("output", 1);

    po::options_description visible_opts;
   	visible_opts.add_options()("help",			"Produce this help message")
        ("inputMap,m", 	      po::value<string>(),	"a text file containing the map information")
        ("inputData,d",       po::value<string>(),	"a text file containing the data to be unmapped")
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
	if (vm.count("output") == 0) {
		cout << "Error: must specify the output filename!" << endl;
		fail = true;
	}

	if (vm.count("inputMap") == 0) {
		cout << "Error: must specify the map filename!" << endl;
		fail = true;
	}

	if (vm.count("inputData") == 0) {
		cout << "Error: must specify the data filename!" << endl;
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

void readMapFromFile(map<int, int> & outMap, string fileName)
{
	ifstream in;
	in.open(fileName);
	string line;
	
	while (getline(in, line)){
		
		vector<string> numbers;
		boost::trim(line);
		boost::split(numbers, line, boost::is_any_of("\t "),
					 boost::token_compress_on);

		outMap.insert(pair<int, int> (stoi(numbers[0]), stoi(numbers[1])));
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


	map<int, int> inputMap;
	readMapFromFile(inputMap, args["inputMap"].as<string>());


	vector<vector<double>> mappedData, unMappedData;
	readTableFromFile(mappedData, args["inputData"].as<string>());

	for (auto & p : inputMap)
		unMappedData.push_back(mappedData[p.second]);
	
	ofstream out;
	out.open(args["output"].as<string>());
	dumpTable(unMappedData, out);


    return 0;
}
