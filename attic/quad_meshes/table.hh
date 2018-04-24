////////////////////////////////////////////////////////////////////////////////
// table.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      a class implementing a table with features like 
//      obtaining unique rows
//      finding a map from the table of unique rows to the actual table
*/ 
//  Author:   Morteza H Siboni (mhs), m.hakimi.siboni@gmail.com
//  Company:  New York University
//  Created:  02/18/2016
////////////////////////////////////////////////////////////////////////////////
#ifndef TABLE_HH
#define TABLE_HH

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <string>
#include <algorithm>
#include <map>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace std;

template<class T, size_t _Dim>
using Matrix = Eigen::Matrix<T, _Dim, _Dim>;

template <class T, size_t Dim = 0>
class table {
public:
	table()	{}

	table(const vector<vector<T>> & inTable)
		:m_table(inTable), m_table_unique(inTable)
	{
		this->absoluteValue();
		this->unique();
		this->generateMap();
	}

	table(const vector<Matrix<T, Dim>> & inTable)
	{
		m_table.clear();
		for (size_t m = 0; m < inTable.size(); ++m)
		{
			vector<T> currentRow;
			Matrix<T, Dim> currentMatrix = inTable[m];
			for (size_t i = 0; i < (size_t)currentMatrix.rows(); ++i)
				for (size_t j = 0; j < (size_t)currentMatrix.cols(); ++j)
					currentRow.push_back(currentMatrix(i,j));
			m_table.push_back(currentRow);
			m_table_unique.push_back(currentRow);
		}

		this->absoluteValue();
		this->unique();
		this->generateMap();
	}

	vector<vector<T>> getTable()
	{
		return m_table;
	}


	vector<vector<T>> getUniqueTable()
	{
		return m_table_unique;
	}

	map<int, int> getMap()
	{
		return m_map;
	}

	static bool vecCompare(vector<T> v1, vector<T> v2)
	{
		for (size_t i = 0; i < v1.size(); ++i)
			if (v1[i] < v2[i])
				return true;
			else if (v1[i] > v2[i])
				return false;
			else 
				;

		return false;
	}

	static bool vecEqual(vector<T> v1, vector<T> v2)
	{
		for (size_t i = 0; i < v1.size(); ++i)
			if (v1[i] != v2[i])
				return false;

		return true;
	}



protected:
	vector<vector<T>> m_table;
	vector<vector<T>> m_table_unique;
	map<int, int>     m_map;

	void absoluteValue()
	{
		for (size_t i = 0; i < m_table.size(); ++i)
			for (size_t j = 0; j < m_table[i].size(); ++j)
				m_table[i][j] = abs(m_table[i][j]);
		
		for (size_t i = 0; i < m_table_unique.size(); ++i)
			for (size_t j = 0; j < m_table_unique[i].size(); ++j)
				m_table_unique[i][j] = abs(m_table_unique[i][j]);
	}


	void unique()
	{
		sort(m_table_unique.begin(), m_table_unique.end(), vecCompare);
		
		typename vector<vector<T>>::iterator it;
		it = unique_copy (m_table_unique.begin(), m_table_unique.end(), m_table_unique.begin(), vecEqual);

		//resize the table
		m_table_unique.resize(distance(m_table_unique.begin(), it));
	}

	void generateMap()
	{
		for (size_t i = 0; i < m_table.size(); ++i)
		{
			typename vector<vector<T>>::iterator it;
			it = search(m_table_unique.begin(), m_table_unique.end(), m_table.begin() + i, m_table.begin() + i + 1, vecEqual);
			int idInShortTable = distance(m_table_unique.begin(), it);
			int idInLongTable  = i;
			pair<int, int> longId2shortId(idInLongTable, idInShortTable);
			m_map.insert(longId2shortId);
		}
	}
};
#endif // TABLE_HH
