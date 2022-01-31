#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>
#include <map>
using namespace std;

template<typename ValueType>
class Trie
{
public:
    Trie();
    ~Trie();
    void reset();
    void insert(const std::string& key, const ValueType& value);
    std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;

      // C++11 syntax for preventing copying and assignment
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;
private:
	struct Node
	{
		char label;
		vector<ValueType>* values = new vector<ValueType>;
		vector<Node *>* children = new vector<Node*>; // the possible children could be a, c, g, t
	};
	void destroy(Node * root);
	void insertHelper(const string&key, const ValueType& value, Node * root);
	vector<ValueType> findHelper(const string& key, bool exactMatchOnly, vector<ValueType> vector, Node * root, int index) const;
	Node * m_root;
};

template <typename ValueType>
Trie<ValueType>::Trie()
{
	m_root = new Node;
}


template <typename ValueType>
void Trie<ValueType>::destroy(Node * root)
{
	if (root == nullptr)
		return;

	delete root->values;

	if (!root->children->empty()) {
		for (auto it = root->children->begin(); it != root->children->end(); it++) {
			destroy(*it);
		}
	}
	
	if (root != nullptr)
		delete root;
}

template <typename ValueType>
Trie<ValueType>::~Trie()
{
	destroy(m_root);
}

template <typename ValueType>
void Trie<ValueType>::insert(const std::string& key, const ValueType& value)
{
	insertHelper(key, value, m_root);
}

template <typename ValueType>
void Trie<ValueType>::insertHelper(const std::string& key, const ValueType& value, Node * root)
{
	if (key == "") { // if we've gone through all the values in the string
		root->values->push_back(value);
		return;
	}

	if (!root->children->empty()) {
		for (auto it = root->children->begin(); it != root->children->end(); it++) {
			if ((*it)->label == key[0]) {
				if (key.size() == 1)
					return insertHelper("", value, *it);
				else
					return insertHelper(key.substr(1), value, *it); // follow that child pointer
			}
		}
	}
	// if we get here, no child node has the label we want
	Node * p = new Node;
	root->children->push_back(p);
	p->label = key[0];
	if (key.size() == 1)
		insertHelper("", value, p);
	else
		insertHelper(key.substr(1), value, p);
}

template <typename ValueType>
vector<ValueType> Trie<ValueType>::find(const string& key, bool exactMatchOnly) const
{
	vector<ValueType> v;
	return findHelper(key, exactMatchOnly, v, m_root, 0);
}

// this function finds exact and non exact matches
template <typename ValueType>
vector<ValueType> Trie<ValueType>::findHelper(const string& key, bool exactMatchOnly, 
	vector<ValueType> vector, Node * root, int index) const
{
	if (root == nullptr) // this is never going to happen from the recursive call, just need to check the base case.
		return vector;

	if (index == key.size())
	{
		for (auto it = root->values->begin(); it != root->values->end(); it++)
			vector.push_back(*it);
	}
	if (!root->children->empty())
	{
		for (auto it = root->children->begin(); it != root->children->end(); it++)
		{
			if ((*it)->label == key[index]) // matches only
				vector = findHelper(key, exactMatchOnly, vector, *it, index + 1);
			else
			{
				if (!exactMatchOnly) // we've hit our only mismatch, now exactMatchOnly
									// has to be true from now on
					vector = findHelper(key, true, vector, *it, index + 1);
			}
		}
	}
	return vector; // if we get here, that means that it doesnt exist...
}		

template <typename ValueType>
void Trie<ValueType>::reset()
{
	destroy(m_root);
	m_root = new Node;
}
#endif // TRIE_INCLUDED
