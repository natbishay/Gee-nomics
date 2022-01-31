#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
#include <fstream>
#include <cassert>
using namespace std;

class GenomeImpl
{
public:
    GenomeImpl(const string& nm, const string& sequence);
    static bool load(istream& genomeSource, vector<Genome>& genomes);
    int length() const;
    string name() const;
    bool extract(int position, int length, string& fragment) const;
private:
	string m_name;
	string m_sequence;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
{
	m_name = nm;
	m_sequence = sequence;
}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes) 
{
	string name;
	string sequence;
	char c;
	bool nameFound = false;
	bool sequenceFound = false;
	bool carrotFound = false;

	while (genomeSource.get(c)) {
		if (!nameFound && !sequenceFound) {
			if (c != '>') {
				if (c == '\n') {
					if (name == "")
						return false;
					nameFound = true;
					continue;
				}
				else if (carrotFound)
					name += c;
			}
			else carrotFound = true;
		}

		if (!sequenceFound && nameFound) {
			if (c != '\n') {
				if (c == '>') {
					sequenceFound = true;
				}
				else {
					if (toupper(c) == 'A' || toupper(c) == 'G' || toupper(c) == 'C' || toupper(c) == 'T' || toupper(c) == 'N')
						sequence += c;
					else return false;
				}
			}
		}

		if (sequenceFound && nameFound) {
			Genome g(name, sequence);
			genomes.push_back(g);

			sequenceFound = false;
			nameFound = false;
			name = "";
			sequence = "";
		}
	}

	if (genomeSource.eof() && nameFound) {
		Genome g(name, sequence);
		genomes.push_back(g);
		sequenceFound = true;
	}

	return nameFound && sequenceFound;
}

int GenomeImpl::length() const
{
    return m_sequence.length(); 
}
\
string GenomeImpl::name() const
{
    return m_name;  
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
	if (position + length <= this->length())
	{
		fragment = m_sequence.substr(position, length);
		return true;
	}
	else return false;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
    m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
    delete m_impl;
}

Genome::Genome(const Genome& other)
{
    m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
    GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
    delete m_impl;
    m_impl = newImpl;
    return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes) 
{
    return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
    return m_impl->length();
}

string Genome::name() const
{
    return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
    return m_impl->extract(position, length, fragment);
}

