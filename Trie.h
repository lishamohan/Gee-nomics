#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>
#include <iostream>

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
	struct Node;
	struct Child
	{
		char m_base;
		Node* m_next;
	};
	struct Node
	{
		std::vector<ValueType> m_values;
		std::vector<Child> m_children;
	};

	Node* m_root;

	void freeTree(Node* cur);
	void findLastMatch(Node*& cur, std::string& sequence);
	void addNodes(Node*& cur, const std::string& sequence);
	void findHelper(Node* cur, const std::string& sequence, bool exactMatchOnly, std::vector<ValueType>& matchValues) const;
};

template<typename ValueType> Trie<ValueType>::Trie()
{
	m_root = new Node;
}

template<typename ValueType> Trie<ValueType>::~Trie()
{
	freeTree(m_root);
}

template<typename ValueType> void Trie<ValueType>::reset()
{
	freeTree(m_root);
	m_root = new Node;
}

template<typename ValueType> void Trie<ValueType>::insert(const std::string& key, const ValueType& value)
{
	Node* cur = m_root;
	std::string remainingBases = key;
	findLastMatch(cur, remainingBases); //sets the cur pointer to that which points to the last node whose base matches the current base in key; this is where we start inserting nodes
	addNodes(cur, remainingBases); //inserts however many nodes are needed to add key to the Trie
	cur->m_values.push_back(value);
}

template<typename ValueType> std::vector<ValueType> Trie<ValueType>::find(const std::string& key, bool exactMatchOnly) const
{
	if (key.empty()) //if key is empty, return values stored in root node
		return m_root->m_values;
	std::vector<ValueType> matchValues;
	int i = 0;
	while (i < m_root->m_children.size() && m_root->m_children[i].m_base != key[0]) //search for the node for the current base in key. If it's not found, return an empty vector. Otherwise, continue searching for the rest of the bases
		i++;
	if (i == m_root->m_children.size())
		return matchValues;
	std::string remainingBases = key.substr(1);
	findHelper(m_root->m_children[i].m_next, remainingBases, exactMatchOnly, matchValues); //search for the rest of the bases
	return matchValues; 
}

template<typename ValueType> void Trie<ValueType>::freeTree(Node* cur)
{
	for (int i = 0; i < cur->m_children.size(); i++)
		freeTree(cur->m_children[i].m_next);
	delete cur;
}

template<typename ValueType> void Trie<ValueType>::findLastMatch(Node*& cur, std::string& sequence)
{
	if (sequence.empty())
		return;
	int i = 0;
	while (i < cur->m_children.size() && cur->m_children[i].m_base != sequence[0])
		i++;
	if (i == cur->m_children.size())
		return;
	cur = cur->m_children[i].m_next;
	sequence = sequence.substr(1);
	findLastMatch(cur, sequence);
}

template<typename ValueType> void Trie<ValueType>::addNodes(Node*& cur, const std::string& sequence)
{
	if (sequence.empty())
		return;
	Child c;
	c.m_base = sequence[0];
	c.m_next = new Node;
	cur->m_children.push_back(c);
	cur = c.m_next;
	addNodes(cur, sequence.substr(1));
}

template<typename ValueType> void Trie<ValueType>::findHelper(Node* cur, const std::string& sequence, bool exactMatchOnly, std::vector<ValueType>& matchValues) const
{
	if (sequence.empty())
	{
		matchValues.insert(matchValues.end(), cur->m_values.begin(), cur->m_values.end());
		return;
	}
	int i = 0;
	std::string remainingBases = sequence.substr(1);
	while (i < cur->m_children.size())
	{
		if (cur->m_children[i].m_base == sequence[0])
			findHelper(cur->m_children[i].m_next, remainingBases, exactMatchOnly, matchValues);
		else if (!exactMatchOnly)
			findHelper(cur->m_children[i].m_next, remainingBases, true, matchValues);
		i++;
	}
}

#endif // TRIE_INCLUDED

