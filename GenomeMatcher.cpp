#include "provided.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include "Trie.h"
using namespace std;

class GenomeMatcherImpl
{
public:
    GenomeMatcherImpl(int minSearchLength);
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
	Trie<pair<int, int>> m_sequences;
	vector<Genome> m_genomes;
	int m_minSearchLength;
};

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
{
	m_minSearchLength = minSearchLength;
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
	m_genomes.push_back(genome);
	for (int pos = 0; pos <= genome.length() - m_minSearchLength; pos++)
	{
		string subsequence;
		genome.extract(pos, m_minSearchLength, subsequence);
		pair<int, int> genomeAndPos;
		genomeAndPos.first = m_genomes.size() - 1;
		genomeAndPos.second = pos;
		m_sequences.insert(subsequence, genomeAndPos);
	}
}

int GenomeMatcherImpl::minimumSearchLength() const
{
	return m_minSearchLength;
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
	if (fragment.length() < minimumLength || minimumLength < m_minSearchLength)
		return false;
	matches.clear();
	vector<pair<int, int>> potentialMatches; 
	string fragOfFragment = fragment.substr(0, m_minSearchLength); //use Trie's find method to find matches for the first minimumSearchLength bases of fragment
	potentialMatches = m_sequences.find(fragOfFragment, exactMatchOnly);
	unordered_map<string, pair<int, int>> matchesMap; //the unordered_map is used and explained below
	for (int i = 0; i < potentialMatches.size(); i++) //go through all the matches returned from Trie::find() and check the remaining bases
	{
		int curGenome = potentialMatches[i].first;
		int origMatchPos = potentialMatches[i].second;
		int mismatchesLeft;
		if (exactMatchOnly) //set a misMatchesLeft counter to keep track of how many mismatches have been encountered and how many are allowed
			mismatchesLeft = 0;
		else
			mismatchesLeft = 1;
		int fragmentPos = 0, lengthOfMatch = 0; 
		int curMatchPos = origMatchPos;
		string nextBase;
		while (fragmentPos < fragment.length()) //extract the rest of the bases of the potential match genome one-by-one, starting at the pos of the match returned from Trie::find(), and compare them to each base of fragment
		{
			if (!m_genomes[curGenome].extract(curMatchPos, 1, nextBase)) //if we can't extract another base because we've reached the end of the current genome, skip rest of this loop (there are no more bases to compare)
				break;
			if (nextBase[0] != fragment[fragmentPos]) //if characters don't match, check if there are any mismatches left. If there are, update the lengthOfMatch and mismatchesLeft conuters and keep going. If there aren't, stop checking this genome
			{
				if (mismatchesLeft == 0)
					break;
				else
					mismatchesLeft--;
			}
			lengthOfMatch++, fragmentPos++, curMatchPos++;
		}
		/*
		if the length of the current match just found is at least minimumLength bases long, we want to add it to the container of matches only if it is the longest 
		and earliest match found for that particular genome. To determine this, we use an unordered map that maps the genome name to a pair containing the length of
		the match and starting position of the match in the genome. For every match found that is at least minimumLength bases long, we will check the unordered map 
		for matches from that particular genome that have already been found. If none have been found, we add the genome name, match length, and match pos to the map; otherwise,
		we use the find() method for the unordered map to find the data stored for that genome's match. We compare its match length and position to the current match length and position
		and keep the longest, earliest match; the other match is deleted (if shorter/later match was in the map) or not added to the map (if shorter/later match is the one we just found).
		*/
		if (lengthOfMatch >= minimumLength) 
		{
			auto it = matchesMap.find(m_genomes[curGenome].name());
			if (it != matchesMap.end())
			{
				if (lengthOfMatch < matchesMap[m_genomes[curGenome].name()].first)
					continue;
				if (lengthOfMatch == matchesMap[m_genomes[curGenome].name()].first && 
					matchesMap[m_genomes[curGenome].name()].second < origMatchPos)
					continue;
				matchesMap.erase(it);
			}
			pair<int, int> lengthAndPos;
			lengthAndPos.first = lengthOfMatch;
			lengthAndPos.second = origMatchPos;
			matchesMap[m_genomes[curGenome].name()] = lengthAndPos;
		}
	}
	if (matchesMap.empty()) //if the the matchesMap is empty at this point, no matches of minimumLength were found 
		return false;
	for (auto it = matchesMap.begin(); it != matchesMap.end(); it++) //add the matches from matchesMap to the matches vector, transfering the data from the map to a DNAMatches struct
	{
		DNAMatch m;
		m.genomeName = (it->first);
		m.length = matchesMap[it->first].first;
		m.position = matchesMap[it->first].second;
		matches.push_back(m);
	}
	return true;
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
	if (fragmentMatchLength < m_minSearchLength)
		return false;
	unordered_map<string, double> potentialMatchesMap; //use a map to keep track of the number of matches sequences for each genome
	int numSequencesConsidered = query.length() / fragmentMatchLength; 
	for (int i = 0; i < numSequencesConsidered; i++) //go through each subsequence of length fragmentMatchLength of query
	{
		string curSequence;
		query.extract(i*fragmentMatchLength, fragmentMatchLength, curSequence);
		vector<DNAMatch> potentialMatches;
		if (!findGenomesWithThisDNA(curSequence, fragmentMatchLength, exactMatchOnly, potentialMatches)) //if there aren't any matches, move on to the next fragment; otherwise, add the match to a map and update the counter
			continue;
		for (int j = 0; j < potentialMatches.size(); j++)
			potentialMatchesMap[potentialMatches[j].genomeName]++;
	}
	vector <pair<int, string>> matches;
	for (auto it = potentialMatchesMap.begin(); it != potentialMatchesMap.end(); it++) //go through the map of matches and compute the percentages. Add those the percentage and genome name of those above the threshold to a vector for sorting
	{
		double p = (it->second / numSequencesConsidered) * 100;
		if (p > matchPercentThreshold)
		{
			pair<double, string> pg;
			pg.first = p;
			pg.second = it->first;
			matches.push_back(pg);
		}
	}
	if (matches.empty()) //at this point, if the matches vector is empty, there were no matches found
		return false;
	results.clear();
	sort(matches.begin(), matches.end(), greater<>()); //sort the vector of matches by descending order according to the percentage
	for (int i = 0; i < matches.size(); i++) //add the data to the results vector. If two genomes have the same percentage match, order them by ascending order of genome name
	{
		GenomeMatch m;
		int addFirst = i;
		int addSecond = -1;
		if (i + 1 < matches.size() && matches[i].first == matches[i + 1].first)
		{
			if (matches[i].second < matches[i + 1].second)
			{
				addSecond = i + 1;
			}
			else
			{
				addFirst = i + 1;
				addSecond = i;
			}
		}
		m.genomeName = matches[addFirst].second;
		m.percentMatch = matches[addFirst].first;
		results.push_back(m);
		if (addSecond != -1)
		{
			GenomeMatch m1;
			m1.genomeName = matches[addSecond].second;
			m1.percentMatch = matches[addSecond].first;
			results.push_back(m1);
			i++;
		}
	}
	return true;
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
    m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
    delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
    m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
    return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}
