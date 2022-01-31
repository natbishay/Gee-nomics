#include "provided.h"
#include <string>
#include <iostream>
#include <fstream>
#include "Trie.h"
#include <unordered_map>
#include <cassert>
#include "provided.h"
#include <iomanip>
#include <vector>
#include <cctype>
#include <cstdlib>
#include <algorithm>
using namespace std;

#if defined(_MSC_VER)  &&  !defined(_DEBUG)
#include <iostream>
#include <windows.h>
#include <conio.h>

struct KeepWindowOpenUntilDismissed
{
	~KeepWindowOpenUntilDismissed()
	{
		DWORD pids[1];
		if (GetConsoleProcessList(pids, 1) == 1)
		{
			std::cout << "Press any key to continue . . . ";
			_getch();
		}
	}
} keepWindowOpenUntilDismissed;
#endif

class GenomeMatcherImpl
{
public:
    GenomeMatcherImpl(int minSearchLength);
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;

private:
	int m_minSearchLength;
	vector<Genome> genomes;

	struct Sequence
	{
		Sequence(unsigned int pos, int positionInGenomeVector) : m_pos(pos), m_positionInGenomeVector(positionInGenomeVector) {}
		unsigned int m_pos;
		int m_positionInGenomeVector;
	};


	Trie <Sequence> trie;
	int lengthOfLongestCommonPrefix(string fragment, string extracted, bool exactMatchOnly) const;
	void hashDNAMatch(DNAMatch d, unordered_map<string, DNAMatch> &hashOfMatches) const;
};

bool sortGenomeMatches(const GenomeMatch& first, const GenomeMatch& second);

 int GenomeMatcherImpl::lengthOfLongestCommonPrefix(string fragment, string extracted, bool exactMatchOnly) const
{
	int length = 0;
	int mismatches = 0;

	int j = 0;
	while (j < fragment.size())
	{
		if (fragment[j] != extracted[j])
		{
			if (j == 0)
				return -1;
			if (exactMatchOnly)
				break;
			else
			{
				mismatches++;
				if (mismatches > 1)
					break;
			}
		}
		j++;
		length++;
	}
	return length;
}

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
{
	m_minSearchLength = minSearchLength;
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
	genomes.push_back(genome); 

	int index = 0;
	string frag;
	
	while (genome.extract(index, minimumSearchLength(), frag)) {
		trie.insert(frag, Sequence(index, genomes.size()-1)); // pass in the address to the genome it references.
		index++;
	}
	
}

int GenomeMatcherImpl::minimumSearchLength() const
{
    return m_minSearchLength;
}

void GenomeMatcherImpl::hashDNAMatch(DNAMatch d, unordered_map<string, DNAMatch> &hashOfMatches) const
{
	auto it = hashOfMatches.find(d.genomeName);													  // and loop through that bucket (should contain VERY few DNAMatch pointers
	if (it == hashOfMatches.end()) // nothing in the bucket
	{
		hashOfMatches.insert({ d.genomeName, d });
		return;
	}
	//for (unordered_map<string, DNAMatch>::const_local_iterator it = hashOfMatches.begin(bucket); 
		//it != hashOfMatches.end(bucket); it++)
	{
		DNAMatch p = it->second;
		// since there could be possible collisions, need to make sure the genome names are the same to make sure we're comparing the right ones!
		if (p.genomeName == d.genomeName) // this IS the genome we want to compare to, and possibly replace in the genome map
		{
			if (p.length < d.length) // if the one that is currently in the hash is smaller than the second, we need to get rid 
			{						// of the one that is in the hash.
				hashOfMatches.erase(it);
				hashOfMatches.insert({ d.genomeName, d }); // put the new item into the hash table
			}
			else if (p.length == d.length) // if the lengths are equal
			{
				if (p.position > d.position) // we want the one with the earlier position
				{
					hashOfMatches.erase(it);
					hashOfMatches.insert({ d.genomeName, d }); // put the new item in the hash table
				}
			}
		}
	}
}
bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
	if (fragment.size() < minimumLength)
		return false;
	if (minimumLength < minimumSearchLength()) // this means the fragment size will always be greater than minimumLength 
		return false;							// and minimumLength will always be greater than minimumSearchLength. 
												// so fragment will always be greater than minimumSearchLength. 
												// so the split up part of fragment of size minimumSearchLength will be smaller than minimumLength.
	
	unordered_map<string, DNAMatch> hashOfMatches;
	
	string frag = fragment.substr(0, minimumSearchLength());
	vector <Sequence> v = trie.find(frag, exactMatchOnly);


	// for each of the extracted sequences in the vector, extract out the fragment size, and check it for the common prefix
	
	for (auto i = 0; i < v.size(); i++)
	{
		
		unsigned int pos = v[i].m_pos;
		const Genome * g = &genomes[v[i].m_positionInGenomeVector];

		string extracted;
		if (g->extract(v[i].m_pos, fragment.size(), extracted)) // there is enough left to extract the whole thing...
		{
			int len = lengthOfLongestCommonPrefix(fragment, extracted, exactMatchOnly);
			if (len < minimumLength)
				continue;
			else
			{
				DNAMatch d;
				d.genomeName = g->name();
				d.position = pos;
				d.length = len;

				hashDNAMatch(d, hashOfMatches);
			}
		}
		else // there isn't enough left to extract the whole thing, like we're at the end of the genome sequence
		{
			int j = 1;
			while (fragment.size() - j >= minimumLength)
			{
				if (g->extract(pos, fragment.size() - j, extracted)) // extract a smaller length that does exist
				{
					int len = lengthOfLongestCommonPrefix(fragment, extracted, exactMatchOnly); // get the matched length
					if (len < minimumLength)
						break; // break out of the while loop, we don't need to check anymore, since the match isn't acceptable anyways
					// i.e. there's no way, when we extract a smaller size, and the len is < minimumLength, that when we extract an even
					// smaller one, we'll get a larger len value.
					else
					{
						DNAMatch d;
						d.genomeName = g->name();
						d.position = pos;
						d.length = len;

						hashDNAMatch(d, hashOfMatches);
						// need to add this somewhere!
						break; // we've found the longest possible length, so let's break
					}
				}
				j++;
			}
		}
	}
	// now we need to get all the items in the hash table...
	// the time complexity is gonna be less than the number of hits!! WOOO! :)
	if (!hashOfMatches.empty())
	{
		for (auto it = hashOfMatches.begin(); it != hashOfMatches.end(); it++) {
			matches.push_back((it->second));
		}
	}

	return !matches.empty(); // if not empty, should return true
							// if empty, should return false.

	// what I need to do:
	// hash table with genome name as a keyword
	// the values of the hash table are a linked list of pointers to DNA matches
	// (which are already in the matches vector)
	// compare to the current DNA match and if it is higher in priority, add it to
	// the hash table and to the matches vector!!!!
	// at the end, the only thing that the matches vector should have are 
	// the highest priority matches.
	// BADDA BING BADDA BOOM
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
	unordered_map<string, int> newHashOfMatches;
	
	int s = query.length() / fragmentMatchLength;
	for (int i = 0; i < query.length(); i += fragmentMatchLength) // O(Q)
	{
		string frag;
		if (query.extract(i, fragmentMatchLength, frag)) // O(1)
		{
			vector<DNAMatch> matches;
			if (findGenomesWithThisDNA(frag, fragmentMatchLength, exactMatchOnly, matches)) // O(X)
			{ // everything else has to be constant time here...
				for (int j = 0; j < matches.size(); j++)
				{ 
					unordered_map<string, int>::iterator it = newHashOfMatches.find(matches[j].genomeName);
					if (it == newHashOfMatches.end()) // no matches yet found for this genomes
					{
						newHashOfMatches.insert({ matches[j].genomeName, 1 }); // so insert the first one
					}
					else
						newHashOfMatches[matches[j].genomeName]++; // increase the matches already found there
				}	
			}
		}
	}

	if (!newHashOfMatches.empty())
	{
		for (unordered_map<string, int>::iterator it = newHashOfMatches.begin(); it != newHashOfMatches.end(); it++) {
			{ 
				double ss = s;
				double val = it->second;
				if (( val / ss) * 100 <= matchPercentThreshold)
					continue;
				GenomeMatch gm; 
				gm.genomeName = it->first;
				gm.percentMatch = (val / ss) * 100;
				results.push_back(gm);
			}
		}
	}
	if (!results.empty())
		sort(results.begin(), results.end(), sortGenomeMatches);

	return !results.empty();
		}

bool sortGenomeMatches(const GenomeMatch& first, const GenomeMatch& second)
{ // ok so we want the first things to be in order of PERCENTAGES first, then order of NAME. 
	if (first.percentMatch != second.percentMatch)
		return first.percentMatch > second.percentMatch;
	// this will return true if the first one's match is LARGER, which is what we want; then first will
	// be earlier in the sort.
	else
		return first.genomeName < second.genomeName;
	// this will return true if the first one's name is smaller, meaning it's alphabetically earlier. 
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



const string PROVIDED_DIR = ".";

const string providedFiles[] = {
	"Ferroplasma_acidarmanus.txt",
	"Halobacterium_jilantaiense.txt",
	"Halorubrum_chaoviator.txt",
	"Halorubrum_californiense.txt",
	"Halorientalis_regularis.txt",
	"Halorientalis_persicus.txt",
	"Ferroglobus_placidus.txt",
	"Desulfurococcus_mucosus.txt"
};

void createNewLibrary(GenomeMatcher*& library)
{
	cout << "Enter minimum search length (3-100): ";
	string line;
	getline(cin, line);
	int len = atoi(line.c_str());
	if (len < 3 || len > 100)
	{
		cout << "Invalid prefix size." << endl;
		return;
	}
	delete library;
	library = new GenomeMatcher(len);
}

void addOneGenomeManually(GenomeMatcher* library)
{
	cout << "Enter name: ";
	string name;
	getline(cin, name);
	if (name.empty())
	{
		cout << "Name must not be empty." << endl;
		return;
	}
	cout << "Enter DNA sequence: ";
	string sequence;
	getline(cin, sequence);
	if (sequence.empty())
	{
		cout << "Sequence must not be empty." << endl;
		return;
	}
	if (sequence.find_first_not_of("ACGTNacgtn") != string::npos)
	{
		cout << "Invalid character in DNA sequence." << endl;
		return;
	}
	for (char ch : sequence)
		ch = toupper(ch);
	library->addGenome(Genome(name, sequence));
}

bool loadFile(string filename, vector<Genome>& genomes)
{
	ifstream inputf(filename);
	if (!inputf)
	{
		cout << "Cannot open file: " << filename << endl;
		return false;
	}
	if (!Genome::load(inputf, genomes))
	{
		cout << "Improperly formatted file: " << filename << endl;
		return false;
	}
	return true;
}

void loadOneDataFile(GenomeMatcher* library)
{
	string filename;
	cout << "Enter file name: ";
	getline(cin, filename);
	if (filename.empty())
	{
		cout << "No file name entered." << endl;
		return;
	}
	vector<Genome> genomes;
	if (!loadFile(filename, genomes))
		return;
	for (const auto& g : genomes)
		library->addGenome(g);
	cout << "Successfully loaded " << genomes.size() << " genomes." << endl;
}

void loadProvidedFiles(GenomeMatcher* library)
{
	for (const string& f : providedFiles)
	{
		vector<Genome> genomes;
		if (loadFile(PROVIDED_DIR + "/" + f, genomes))
		{
			for (const auto& g : genomes)
				library->addGenome(g);
			cout << "Loaded " << genomes.size() << " genomes from " << f << endl;
		}
	}
}

void findGenome(GenomeMatcher* library, bool exactMatch)
{
	if (exactMatch)
		cout << "Enter DNA sequence for which to find exact matches: ";
	else
		cout << "Enter DNA sequence for which to find exact matches and SNiPs: ";
	string sequence;
	getline(cin, sequence);
	int minLength = library->minimumSearchLength();
	if (sequence.size() < minLength)
	{
		cout << "DNA sequence length must be at least " << minLength << endl;
		return;
	}
	cout << "Enter minimum sequence match length: ";
	string line;
	getline(cin, line);
	int minMatchLength = atoi(line.c_str());
	if (minMatchLength > sequence.size())
	{
		cout << "Minimum match length must be at least the sequence length." << endl;
		return;
	}
	vector<DNAMatch> matches;
	if (!library->findGenomesWithThisDNA(sequence, minMatchLength, exactMatch, matches))
	{
		cout << "No ";
		if (exactMatch)
			cout << " matches";
		else
			cout << " matches or SNiPs";
		cout << " of " << sequence << " were found." << endl;
		return;
	}
	cout << matches.size();
	if (exactMatch)
		cout << " matches";
	else
		cout << " matches and/or SNiPs";
	cout << " of " << sequence << " found:" << endl;
	for (const auto& m : matches)
		cout << "  length " << m.length << " position " << m.position << " in " << m.genomeName << endl;
}

bool getFindRelatedParams(double& pct, bool& exactMatchOnly)
{
	cout << "Enter match percentage threshold (0-100): ";
	string line;
	getline(cin, line);
	pct = atof(line.c_str());
	if (pct < 0 || pct > 100)
	{
		cout << "Percentage must be in the range 0 to 100." << endl;
		return false;
	}
	cout << "Require (e)xact match or allow (S)NiPs (e or s): ";
	getline(cin, line);
	if (line.empty() || (line[0] != 'e' && line[0] != 's'))
	{
		cout << "Response must be e or s." << endl;
		return false;
	}
	exactMatchOnly = (line[0] == 'e');
	return true;
}

void findRelatedGenomesManual(GenomeMatcher* library)
{
	cout << "Enter DNA sequence: ";
	string sequence;
	getline(cin, sequence);
	int minLength = library->minimumSearchLength();
	if (sequence.size() < minLength)
	{
		cout << "DNA sequence length must be at least " << minLength << endl;
		return;
	}
	double pctThreshold;
	bool exactMatchOnly;
	if (!getFindRelatedParams(pctThreshold, exactMatchOnly))
		return;

	vector<GenomeMatch> matches;
	library->findRelatedGenomes(Genome("x", sequence), 2 * minLength, exactMatchOnly, pctThreshold, matches);
	if (matches.empty())
	{
		cout << "    No related genomes were found" << endl;
		return;
	}
	cout << "    " << matches.size() << " related genomes were found:" << endl;
	cout.setf(ios::fixed);
	cout.precision(2);
	for (const auto& m : matches)
		cout << " " << setw(6) << m.percentMatch << "%  " << m.genomeName << endl;
}

void findRelatedGenomesFromFile(GenomeMatcher* library)
{
	string filename;
	cout << "Enter name of file containing one or more genomes to find matches for: ";
	getline(cin, filename);
	if (filename.empty())
	{
		cout << "No file name entered." << endl;
		return;
	}
	vector<Genome> genomes;
	if (!loadFile(filename, genomes))
		return;
	double pctThreshold;
	bool exactMatchOnly;
	if (!getFindRelatedParams(pctThreshold, exactMatchOnly))
		return;

	int minLength = library->minimumSearchLength();
	for (const auto& g : genomes)
	{
		vector<GenomeMatch> matches;
		library->findRelatedGenomes(g, 2 * minLength, exactMatchOnly, pctThreshold, matches);
		cout << "  For " << g.name() << endl;
		if (matches.empty())
		{
			cout << "    No related genomes were found" << endl;
			continue;
		}
		cout << "    " << matches.size() << " related genomes were found:" << endl;
		cout.setf(ios::fixed);
		cout.precision(2);
		for (const auto& m : matches)
			cout << "     " << setw(6) << m.percentMatch << "%  " << m.genomeName << endl;
	}
}

void showMenu()
{
	cout << "        Commands:" << endl;
	cout << "         c - create new genome library      s - find matching SNiPs" << endl;
	cout << "         a - add one genome manually        r - find related genomes (manual)" << endl;
	cout << "         l - load one data file             f - find related genomes (file)" << endl;
	cout << "         d - load all provided data files   ? - show this menu" << endl;
	cout << "         e - find matches exactly           q - quit" << endl;
}

int main()
{
	const int defaultMinSearchLength = 10;

	cout << "Welcome to the Gee-nomics test harness!" << endl;
	cout << "The genome library is initially empty, with a default minSearchLength of " << defaultMinSearchLength << endl;
	showMenu();

	GenomeMatcher* library = new GenomeMatcher(defaultMinSearchLength);

	for (;;)
	{
		cout << "Enter command: ";
		string command;
		if (!getline(cin, command))
			break;
		if (command.empty())
			continue;
		switch (tolower(command[0]))
		{
		default:
			cout << "Invalid command " << command << endl;
			break;
		case 'q':
			delete library;
			return 0;
		case '?':
			showMenu();
			break;
		case 'c':
			createNewLibrary(library);
			break;
		case 'a':
			addOneGenomeManually(library);
			break;
		case 'l':
			loadOneDataFile(library);
			break;
		case 'd':
			loadProvidedFiles(library);
			break;
		case 'e':
			findGenome(library, true);
			break;
		case 's':
			findGenome(library, false);
			break;
		case 'r':
			findRelatedGenomesManual(library);
			break;
		case 'f':
			findRelatedGenomesFromFile(library);
			break;
		}
	}
}



