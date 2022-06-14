#ifndef DATASTRUCT_H
#define DATASTRUCT_H
#include<vector>
#include<iostream>
#include<unordered_map>
#include<unordered_set>
#include<string>

#include "global.h"

using namespace std;

struct Node {
	vector<int> pred;
	vector<int> succ;
};

struct Graph {
	vector<Node> data;
	vector<int> nodeMap;
	int nodeNum;
	int edgeNum;
	Graph(int num) {
		nodeNum = 0;
		edgeNum = 0;
		data.resize(num + 1);
		nodeMap.resize(num + 1);
	}
};

struct HashNode {
	unordered_set<int> pred;
	unordered_set<int> succ;
};

struct LRUSet {
	vector<int> nodes;
	unordered_map<int, int> node2Idx;

	bool insert(int val) {
		if (node2Idx.count(val)) return false;
		node2Idx[val] = nodes.size();
		nodes.push_back(val);
		return true;
	}

	void erase(int val) {
		int idx = node2Idx[val];
		node2Idx[nodes.back()] = idx;
		nodes[idx] = nodes.back();
		nodes.pop_back();
		node2Idx.erase(val);
	}

	bool count(int u) {
		return node2Idx.count(u);
	}

	size_t size() {
		return nodes.size();
	}
}; 

struct LRUNode {
	LRUSet pred;
	LRUSet succ;
};

struct LRUGraph {
	vector<int> nodes;
	vector<LRUNode> node2Edge;
	vector<int> node2Idx;
	void resize(int num) {
		nodes.resize(num);
		node2Edge.resize(num + 1);
		node2Idx.resize(num + 1);
		for (int i = 0; i < num; ++i) {
			nodes[i] = i + 1;
			node2Idx[i + 1] = i;
		}
	}

	void erase(int u) {
		int idx = node2Idx[u];
		node2Idx[nodes.back()] = idx;
		nodes[idx] = nodes.back();
		nodes.pop_back();
		node2Idx[u] = -1;
	}

	LRUNode& operator[](int idx) {
		return node2Edge[idx];
	}
};

struct DLinkedNode {
	LL idx;
	int prev;
	int next;
	DLinkedNode(LL _idx = -1) {
		idx = _idx;
		prev = 0;
		next = 0;
	}
};

#endif // !DATASTRUCT_H
