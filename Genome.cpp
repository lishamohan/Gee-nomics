#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
#include <fstream>
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
	genomes.clear();
	string s, name, sequence;
	while (getline(genomeSource, s))
	{
		if (s.length() == 0) //checks if line is empty
			return false;
		if (s[0] == '>')
	 	{
			if (!name.empty() && sequence.empty()) //checks for two name lines without bases in between
				return false;
			if (!name.empty() && !sequence.empty()) //we've reached the next genome, so add the previous genome to the vector
			{
				Genome g(name, sequence);
				genomes.push_back(g);
				name.clear();
				sequence.clear();
			}
			if (s.length() == 1) //checks if name line has at least one character other than '>'
				return false;
			name = s.substr(1);
			if (!getline(genomeSource, s)) //checks that there are bases after the name line
				return false;
		}

		if (name.empty()) //checks if there are bases with no name line before them
			return false;
		for (int i = 0; i < s.size(); i++) //checks if bases are valid
		{
			char c = toupper(s[i]);
			if (c != 'A' && c != 'C' && c != 'T' && c != 'G' && c != 'N')
				return false;
			sequence += c;
		}
	}
	Genome g(name, sequence); //adds the last genome to the vector
	genomes.push_back(g);
	return true;
}

int GenomeImpl::length() const
{
	return m_sequence.length();
}

string GenomeImpl::name() const
{
	return m_name;
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
	if (position + length > m_sequence.length())
		return false;
	fragment = m_sequence.substr(position, length);
	return true;
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
