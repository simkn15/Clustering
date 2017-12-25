#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>
#include <limits>
#include <cmath>
#include <utility>
#include <stdlib.h>
#include "transclust/TriangularMatrix.hpp"

TriangularMatrix::TriangularMatrix(
		const std::string &filename,
		TCC::TransClustParams& tcp,
		unsigned cc_id)
{
	id = cc_id;
	// map <object name> -> <index in matrix>
	std::map<std::string, unsigned> object2index;

	// map for similarity value <<object1 name>-<object2 name>> -> <value>
	std::map<std::pair<std::string, std::string>, float> sim_value;

	// map for findeing one-way sim values
	std::map<std::pair<std::string, std::string>, bool> hasPartner;


	// reads the input file and initializes the 'matrix' array
	readFile(filename,object2index,sim_value,hasPartner);

	if(tcp.file_type == "LEGACY"){

		/* Init 1d matrix  */
		parseLegacySimDataFile(
				object2index,
				sim_value,
				hasPartner,
				tcp);
	}else if(tcp.file_type == "SIMPLE"){

		parseSimpleSimDataFile(
				object2index,
				sim_value,
				hasPartner,
				tcp);
	}
}

TriangularMatrix::TriangularMatrix(
		const TriangularMatrix &m,
		const std::vector<unsigned> &objects,
		unsigned cc_id)
{
	id = cc_id;
	maxValue = std::numeric_limits<float>::lowest();
	minValue = std::numeric_limits<float>::max();

	num_o = objects.size();
	unsigned msize = ((num_o * num_o) - num_o) / 2;

	matrix.resize(msize);
	for (unsigned i = 0; i < num_o; i++)
	{
		// indexes
		index2ObjName.push_back(m.getObjectName(objects.at(i)));
		index2ObjId.push_back(m.getObjectId(objects.at(i)));
		for (unsigned j = i + 1; j < num_o; j++)
		{
			float val = m(objects.at(i), objects.at(j));
			if (maxValue < val){maxValue = val;}
			if (minValue > val){minValue = val;}

			matrix.at(index(i, j)) = val;
		}
	}
}

void TriangularMatrix::parseLegacySimDataFile(
		std::map<std::string, unsigned> &object2index,
		std::map<std::pair<std::string, std::string>, float> & sim_value,
		std::map<std::pair<std::string, std::string>, bool> &hasPartner,
		TCC::TransClustParams& tcp)
{


	// for each pair of object, tjeck if it has a partner and if so assign the
	// minimum similarity value. If a pair only has similarity data in one
	// direction then assign the defautl similarity value
	num_o = index2ObjName.size();
	unsigned msize = ((num_o * num_o) - num_o) / 2;
	matrix.resize(msize);
	for (unsigned i = 0; i < num_o; i++)
	{
		for (unsigned j = i + 1; j < num_o; j++)
		{
			std::pair<std::string, std::string> key;
			key = std::make_pair(index2ObjName[i], index2ObjName[j]);

			// default value if pair is missing
			float val = tcp.sim_fallback;

			// if the pair is in simvalue and it has a partner
			// choose the smallest sim value
			// else
			// use defualt value
			if(sim_value.find(key) != sim_value.end())
			{
				if(hasPartner.find(key) != hasPartner.end())
				{
					std::pair<std::string, std::string> rkey;
					rkey = std::make_pair(index2ObjName[j], index2ObjName[i]);
					val = std::min(sim_value[key],sim_value[rkey]);
				}
			}
			if(maxValue < val){maxValue = val;};
			matrix.at(index(i, j)) = val;
		}
	}

}

void TriangularMatrix::parseSimpleSimDataFile(
		std::map<std::string, unsigned> &object2index,
		std::map<std::pair<std::string, std::string>, float> & sim_value,
		std::map<std::pair<std::string, std::string>, bool> &hasPartner,
		TCC::TransClustParams& tcp)
{

	std::vector<std::pair<unsigned,unsigned>> positive_inf;
	num_o = index2ObjName.size();
	unsigned long msize = ((num_o * num_o) - num_o) / 2;
		matrix.resize(msize);

		for (unsigned i = 0; i < num_o; i++)
		{
			for (unsigned j = i + 1; j < num_o; j++)
			{
				float val = parseSimpleEdge(positive_inf,object2index,sim_value,hasPartner,i,j,tcp);
				matrix.at(index(i, j)) = val;
			}
		}
		for(std::pair<unsigned,unsigned> &p:positive_inf){
			matrix.at(index(p.first,p.second)) = maxValue;
		}
}

float TriangularMatrix::parseSimpleEdge(
		std::vector<std::pair<unsigned,unsigned>>& positive_inf,
		std::map<std::string, unsigned>& object2index,
		std::map<std::pair<std::string, std::string>, float> & sim_value,
		std::map<std::pair<std::string, std::string>, bool> &hasPartner,
		unsigned i,
		unsigned j,
		TCC::TransClustParams& tcp)
{
	std::pair<std::string, std::string> key;
	key = std::make_pair(index2ObjName[i], index2ObjName[j]);

	float val = tcp.sim_fallback;

	if(sim_value.find(key) != sim_value.end())
	{
		float o1_o2 = sim_value[key];

		if(hasPartner[key])
		{
			std::pair<std::string, std::string> rkey;
			rkey = std::make_pair(index2ObjName[j], index2ObjName[i]);

			float o2_o1 = sim_value[rkey];
			if(std::min(o1_o2,o2_o1) == std::numeric_limits<float>::lowest()){
				val = std::max(o1_o2,o2_o1);
			}else{
				val = std::min(o1_o2,o2_o1);
			}
		}else{
			val = o1_o2;
		}
	}

	if(val == std::numeric_limits<float>::max()){
		positive_inf.push_back(std::make_pair(i,j));
	}else if(maxValue < val){
		maxValue = val;
	}
	return val;

}

void TriangularMatrix::readFile(
	const std::string &filename,
	std::map<std::string, unsigned> &object2index,
	std::map<std::pair<std::string, std::string>, float> & sim_value,
	std::map<std::pair<std::string, std::string>, bool> &hasPartner
	)
{
	maxValue = std::numeric_limits<float>::lowest();
	minValue = std::numeric_limits<float>::max();
	/*********************************
	 * Read the input similarity file
	 ********************************/
	std::ifstream fs(filename);
	for (std::string line; std::getline(fs, line); )
	{
		// split line from similarity file <o1,o2,sim val>
		std::istringstream buf(line);
		std::istream_iterator<std::string> beg(buf), end;
		std::vector<std::string> tokens(beg, end);

		// object 1
		std::string o1 = tokens.at(0);
		// object 2
		std::string o2 = tokens.at(1);

		// o1,o2 similarity
		//float value = std::atof(tokens.at(2).c_str());
		float value = std::stod(tokens.at(2));

		// if inf/-inf
		if(std::isinf(value))
		{
			if(value < 0 ){
				value = std::numeric_limits<float>::lowest();
			}else{
				value = std::numeric_limits<float>::max();
			}
		}

		// build objects maps
		if (object2index.find(o1) == object2index.end())
		{
			index2ObjName.push_back(o1);
			unsigned _id = index2ObjName.size()-1;
			object2index[o1] = _id;
			index2ObjId.push_back(_id);
		}

		if (object2index.find(o2) == object2index.end())
		{
			index2ObjName.push_back(o2);
			unsigned _id = index2ObjName.size()-1;
			object2index[o2] = _id;
			index2ObjId.push_back(_id);
		}

		// create key for similarity map (smallest id first)
		std::pair<std::string, std::string> key;
		if (object2index[o1] < object2index[o2])
		{
			key = std::make_pair(o1, o2);
		}
		else
		{
			key = std::make_pair(o2, o1);
		}

		// if a similarity value for the two objects exists
		// choose the smallest
		if (sim_value.find(key) != sim_value.end())
		{
			hasPartner[key] = true;
		}else{
			hasPartner[key] = false;
		}

		if (minValue > value && value != std::numeric_limits<float>::lowest())
		{
			minValue = value;
		}
		sim_value[std::make_pair(o1,o2)] = value;
	}
	fs.close();

}

TriangularMatrix::TriangularMatrix(
		std::vector<float>& sim_matrix_1d,
		unsigned _num_o,
		unsigned cc_id
   )
{
	id = cc_id;
   num_o = _num_o;
	matrix = sim_matrix_1d;

	maxValue = std::numeric_limits<float>::lowest();
	minValue = std::numeric_limits<float>::max();

	index2ObjName.resize(num_o);
	index2ObjId.resize(num_o);
	for(unsigned i = 0; i < index2ObjName.size();i++){
		index2ObjName.at(i) = std::to_string(i);
		index2ObjId.at(i) = i;
	}
	for(unsigned i = 0; i < matrix.size(); i++){
			float val = matrix.at(i);
			if (maxValue < val){maxValue = val;}
			if (minValue > val){minValue = val;}
	}
}
