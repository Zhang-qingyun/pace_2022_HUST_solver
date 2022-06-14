#include <fstream>
#include <unordered_set>
#include <queue>
#include <sstream>
#include <cassert>
#include <functional>

#include "SA.h"
#include "LocalSearch.h"

SA::SA() {}

void SA::init()
{
	rand = Random::getInstance();
}

bool SA::reduceFacade(int move)
{
	bool result = false;

	switch (move) {
	case 1:
		result = in0();
		break;
	case 2:
		result = out0();
		break;
	case 3:
		result = in1();
		break;
	case 4:
		result = out1();
		break;
	case 5:
		result = loop();
		break;
	case 6:
		result = pie();
		break;
	case 7:
		result = core();
		break;
	case 8:
		result = dome();
		break;
	/*case 9:
		result = mostPIE();
		break;*/
	}

	return result;
}

void SA::reduceGraph()
{
	int MAXMOVE1 = 5;
	int MAXMOVE2 = 8;
	nodeMap.resize(orgNodeNum + 1, 1);
	
	auto reduce = [&]() {
		bool flag2 = true;
		while (flag2)
		{
			flag2 = false;
			for (int move = 1; move <= MAXMOVE1; ++move) {
				if (reduceFacade(move)) {
					flag2 = true;
				}
			}
		}
	};

	// 先调用 five reduce
	reduce();

	bool flag1 = true;
	while (flag1) {
		flag1 = false;
		for (int move = MAXMOVE1 + 1; move <= MAXMOVE2; ++move) {
			if (reduceFacade(move)) {
				flag1 = true;
				reduce();
			}
		}
	}
}

bool SA::in0()
{
	bool flag = false;
	for (int i = 1; i <= orgNodeNum; ++i) {
		if (nodeMap[i] == 1 && org[i].pred.size() == 0) {
			nodeMap[i] = 0;
			for (auto adj : org[i].succ) {
				org[adj].pred.erase(i);
			}
			flag = true;
		}
	}
	return flag;
}

bool SA::out0()
{
	bool flag = false;
	for (int i = 1; i <= orgNodeNum; ++i) {
		if (nodeMap[i] == 1 && org[i].succ.size() == 0) {
			nodeMap[i] = 0;
			for (auto adj : org[i].pred) {
				org[adj].succ.erase(i);
			}
			flag = true;
		}
	}
	return flag;
}

bool SA::in1()
{
	bool flag = false;
	for (int i = 1; i <= orgNodeNum; ++i) {
		if (nodeMap[i] == 1 && org[i].pred.size() == 1) {
			int u = *org[i].pred.begin();
			if (u == i) continue;
			nodeMap[i] = 0;
			org[u].succ.erase(i);
			for (auto adj : org[i].succ) {
				org[adj].pred.erase(i);
				org[u].succ.insert(adj);
				org[adj].pred.insert(u);
			}
			flag = true;
		}
	}
	return flag;
}

bool SA::out1()
{
	bool flag = false;
	for (int i = 1; i <= orgNodeNum; ++i) {
		if (nodeMap[i] == 1 && org[i].succ.size() == 1) {
			int u = *org[i].succ.begin();
			if (u == i) continue;
			nodeMap[i] = 0;
			org[u].pred.erase(i);
			for (auto adj : org[i].pred) {
				org[adj].succ.erase(i);
				org[u].pred.insert(adj);
				org[adj].succ.insert(u);
			}
			flag = true;
		}
	}
	return flag;
}

bool SA::loop()
{
	bool flag = false;
	for (int i = 1; i <= orgNodeNum; ++i) {
		if (nodeMap[i] == 1 && org[i].succ.count(i)) {
			nodeMap[i] = -1;
			flag = true;
			for (auto adj : org[i].pred) {
				org[adj].succ.erase(i);
			}
			for (auto adj : org[i].succ) {
				org[adj].pred.erase(i);
			}
		}
	}
	return flag;
}

bool SA::pie()
{
	bool flag = false;

	// 标记PIE边
	vector<vector<int>> tmp(orgNodeNum + 1);

	for (int i = 1; i <= orgNodeNum; ++i) {
		if (nodeMap[i] == 1) {
			for (auto adj : org[i].succ) {
				if (!org[i].pred.count(adj)) {
					tmp[i].push_back(adj);
				}
			}
		}
	}

	vector<bool> instack(orgNodeNum + 1, false);
	vector<int> dfn(orgNodeNum + 1, 0);
	vector<int> low(orgNodeNum + 1, 0);
	vector<int> belong(orgNodeNum + 1, 0);
	vector<int> parent(orgNodeNum + 1, 0);
	stack<int> s;
	stack<int> trace;
	int cnt = 0, cntb = 0;

	std::function<void(int)> tarjan = [&](int u)
	{
		++cnt;
		dfn[u] = low[u] = cnt;
		s.push(u);
		instack[u] = true;
		for (auto v : tmp[u]) {
			if (!dfn[v])
			{
				tarjan(v);
				low[u] = min(low[u], low[v]);
			}
			else if (instack[v])
				low[u] = min(low[u], dfn[v]);
		}
		if (dfn[u] == low[u])
		{
			++cntb;
			int node;
			do
			{
				node = s.top();
				s.pop();
				instack[node] = false;
				low[u] = INT_MAX;
				belong[node] = cntb;
			} while (node != u);
		}
	};

	std::function<void(int)> tarjanIteration = [&](int u)
	{
		s.push(u);
		parent[u] = u;
		while (!s.empty()) {
			u = s.top();
			if (dfn[u] == 0) {
				trace.push(u);
				++cnt;
				dfn[u] = low[u] = cnt;
			}
			bool flag = false;
			for (auto v : tmp[u]) {
				parent[v] = u;
				if (dfn[v] == 0)
				{
					flag = true;
					s.push(v);
					/*tarjan(v);
					low[u] = min(low[u], low[v]);*/
				}
				else {
					low[u] = min(low[u], low[v]);
				}	
			}
			if (!flag) {
				u = s.top();
				s.pop();
				low[parent[u]] = min(low[parent[u]], low[u]);
				if (dfn[u] == low[u])
				{
					++cntb;
					int node;
					do
					{
						node = trace.top();
						trace.pop();
						low[node] = INT_MAX;
						belong[node] = cntb;
					} while (node != u);
				}
			}
		}
	};

	for (int i = 1; i <= orgNodeNum; ++i) {
		if (nodeMap[i] == 1 && dfn[i] == 0) {
			//tarjan(i);
			tarjanIteration(i);
		}
	}

	// 遍历图，查找无效边
	for (int i = 1; i <= orgNodeNum; ++i) {
		if (nodeMap[i] == 1) {
			for (auto it = org[i].pred.begin(); it != org[i].pred.end();) {
				if (!org[i].succ.count(*it) && belong[i] != belong[*it]) {
					flag = true;
					org[*it].succ.erase(i);
					it = org[i].pred.erase(it);
				}
				else {
					it++;
				}
			}
		}
	}

	return flag;
}

bool SA::core()
{
	bool flag = false;
	vector<bool> invalid(orgNodeNum + 1, false);
	vector<vector<int>> piv;
	for (int i = 1; i <= orgNodeNum; ++i) {
		if (nodeMap[i] == 1 && isPIV(i)) {
			int dgr = org[i].pred.size() + org[i].succ.size();
			piv.push_back({ i, dgr });
		}
	}
	// 对PIV节点按度进行排序
	sort(piv.begin(), piv.end(), [](const vector<int>& a, const vector<int>& b) {
		return a[1] < b[1];
	});
	vector<int> rm;
	for (auto it : piv) {
		int i = it[0];
		if (invalid[i]) continue;
		vector<int> cand;
		for (auto adj : org[i].pred) {
			cand.push_back(adj);
			invalid[adj] = true;
		}
		if (isClique(cand)) {
			flag = true;
			nodeMap[i] = 0;
			for (auto node : cand) {
				nodeMap[node] = -1;
				for (auto adj : org[node].pred) {
					org[adj].succ.erase(node);
				}
				for (auto adj : org[node].succ) {
					org[adj].pred.erase(node);
				}
			}
		}
	}

	return flag;
}

bool SA::dome()
{
	bool flag = false;

	auto satifyPred = [&](int u, int v) {
		for (auto pu : org[u].pred) {
			if (org[u].succ.count(pu)) continue;
			if (!org[v].pred.count(pu)) {
				return false;
			}
		}
		return true;
	};

	auto satifySucc = [&](int u, int v) {
		for (auto pv : org[v].succ) {
			if (org[v].pred.count(pv)) continue;
			if (!org[u].succ.count(pv)) {
				return false;
			}
		}
		return true;
	};

	vector <pair<int, int>> delEdge;
	for (int i = 1; i <= orgNodeNum; ++i) {
		if (nodeMap[i] == 1) {
			for (auto adj : org[i].succ) {
				if (org[i].pred.count(adj)) continue;
				if (satifyPred(i, adj) || satifySucc(i, adj)) {
					flag = true;
					delEdge.emplace_back(i, adj);
				}
			}
		}
	}

	for (auto it : delEdge) {
		org[it.second].pred.erase(it.first);
		org[it.first].succ.erase(it.second);
	}
	return flag;
}

void SA::restoreGraph(int gidx, vector<DLinkedNode>& sol)
{
	for (int i = 1; i <= subgraphs[gidx].nodeNum; ++i) {
		int idx = subgraphs[gidx].nodeMap[i];
		nodeMap[idx] = sol[i].idx;
	}
	nodeMap[0] += sol[0].idx;
}

bool SA::isPIV(int u)
{
	for (auto adj : org[u].pred) {
		assert(nodeMap[adj] == 1);
		if (!org[u].succ.count(adj)) {
			return false;
		}
	}
	for (auto adj : org[u].succ) {
		if (!org[u].pred.count(adj)) {
			return false;
		}
	}
	return true;
}

bool SA::isClique(vector<int>& cand)
{
	for (size_t i = 0; i < cand.size(); ++i) {
		for (size_t j = i + 1; j < cand.size(); ++j) {
			if (!org[cand[i]].succ.count(cand[j])) {
				return false;
			}
			if (!org[cand[i]].pred.count(cand[j])) {
				return false;
			}
		}
	}
	return true;
}

bool SA::allPIEdgeJudge(int gidx)
{
	auto& subgraph = subgraphs[gidx];
	for (int i = 1; i <= subgraph.nodeNum; ++i) {
		int u = subgraph.nodeMap[i];
		for (auto adj : subgraph.data[i].succ) {
			int v = subgraph.nodeMap[adj];
			if (!org[v].succ.count(u)) {
				return false;
			}
		}
	}
	return true;
}

int SA::calConnect()
{
	UnionFind ufset;
	ufset.init(orgNodeNum);
	// 生成新图
	nodeMap[0] = 0;
	reduceNodeNum = 0;
	reduceEdgeNum = 0;
	for (int i = 1; i <= orgNodeNum; ++i) {
		if (nodeMap[i] > 0) {
			reduceNodeNum++;
			nodeMap[i] = reduceNodeNum;
			for (auto adj : org[i].succ) {
				assert(nodeMap[adj] > 0);
				ufset.merge(i, adj);
				reduceEdgeNum++;
			}
		}
		if (nodeMap[i] == -1) {
			nodeMap[0]++;
		}
	}

	#if !SUBMIT
	printf("Reduce Graph: node %d edge %d loop %d\n", reduceNodeNum, reduceEdgeNum, nodeMap[0]);
	#endif // SUBMIT

	vector<int> nd2GIdx(orgNodeNum + 1, -1);	
	for (int i = 1; i <= orgNodeNum; ++i) {
		if (nodeMap[i] > 0 && ufset.f[i] == i) {
			nd2GIdx[i] = subgraphs.size();
			subgraphs.push_back(Graph(ufset.num[i]));
		}
	}
	for (int i = 1; i <= orgNodeNum; ++i) {
		if (nodeMap[i] > 0) {
			int gidx = nd2GIdx[ufset.find(i)];
			nd2GIdx[i] = gidx;
			subgraphs[gidx].nodeNum++;
			int idx = subgraphs[gidx].nodeNum;
			nodeMap[i] = idx;
			subgraphs[gidx].nodeMap[idx] = i;
		}
	}

	for (int i = 1; i <= orgNodeNum; ++i) {
		if (nodeMap[i] <= 0) continue;
		int gidx = nd2GIdx[i];
		auto& graph = subgraphs[gidx];
		for (auto adj : org[i].succ) {
			assert(nodeMap[adj] > 0);
			assert(nd2GIdx[adj] == gidx);
			graph.data[nodeMap[i]].succ.push_back(nodeMap[adj]);
			graph.data[nodeMap[adj]].pred.push_back(nodeMap[i]);
			graph.edgeNum++;
		}
	}

	return 0;
}

void SA::solve()
{
	init();
	reduceGraph();

	calConnect();
	//return;

	clock_t start = clock();

	vector<LS> sols;
	for (int i = 0; i < subgraphs.size(); ++i) {
		int nodeNum = subgraphs[i].nodeNum;
		int edgeNum = subgraphs[i].edgeNum;
		vector<Node>& graph = subgraphs[i].data;
		sols.push_back(LS(nodeNum, edgeNum, graph));
		sols.back().subgraphInit();
	}

	while (!Config::getInstance()->TLE()) {
		for (int i = 0; i < subgraphs.size(); ++i) {
			if (Config::getInstance()->TLE()) break;
			sols[i].subgraphExe();
		}
	}
	
	int totalIter = 0;
	for (int i = 0; i < subgraphs.size(); ++i) {
		totalIter += sols[i].getTotalIter();
		restoreGraph(i, sols[i].getCurBest());
	}
	#if !SUBMIT
	printf("iter: %d init time: %d ", totalIter, (clock() - start) / CLOCKS_PER_SEC);
	#endif // SUBMIT
}

bool SA::processInstanceFiles(char* dataFile, char* instCode)
{
	ifstream data;
	data.open(dataFile, ios::in);

	if (!data.is_open()) {
		printf("Erro ao abrir o arquivo = %s \n\n", dataFile);
		return false;
	}

	this->instCode = instCode;
	this->instString = dataFile;

	// Problem informations
	data >> this->orgNodeNum >> this->orgEdgeNum;
	printf("inst: %s nodeNum: %d edgeNum: %d\n", instCode, this->orgNodeNum, this->orgEdgeNum);
	org.resize(this->orgNodeNum + 1);

	int x, y;
	for (int i = 0; i < this->orgEdgeNum; ++i) {
		data >> x >> y;
		org[x].succ.insert(y);
		org[y].pred.insert(x);
		assert(x != y);
	}


	data.close();
	return true;
}

bool SA::processContestFiles(char* dataFile, char* instCode)
{
	ifstream data;
	data.open(dataFile, ios::in);

	if (!data.is_open()) {
		printf("Erro ao abrir o arquivo = %s \n\n", dataFile);
		return false;
	}

	this->instCode = instCode;
	this->instString = dataFile;

	int wgt;
	// Problem informations
	data >> this->orgNodeNum >> this->orgEdgeNum >> wgt;
	printf("inst: %s nodeNum: %d edgeNum: %d weight: %d\n", instCode, this->orgNodeNum, this->orgEdgeNum, wgt);
	org.resize(this->orgNodeNum + 1);

	string str;
	getline(data, str);
	int y;
	for (int x = 1; x <= this->orgNodeNum; ++x) {
		getline(data, str);
		while(str[0] == '%')getline(data, str);
		stringstream ss(str);
		while (ss >> y) {
			org[x].succ.insert(y);
			org[y].pred.insert(x);
			assert(x != y);
		}
	}
	data.close();
	return true;
}

bool SA::processCmd()
{
//    freopen("../instances/h_193", "r", stdin);
	this->instCode = "";
	this->instString = "";

	int wgt;
	// Problem informations
	cin >> this->orgNodeNum >> this->orgEdgeNum >> wgt;
	org.resize(this->orgNodeNum + 1);

	string str;
	getline(cin, str);
	int y;
	for (int x = 1; x <= this->orgNodeNum; ++x) {
		getline(cin, str);
		while (str[0] == '%')getline(cin, str);
		stringstream ss(str);
		while (ss >> y) {
			org[x].succ.insert(y);
			org[y].pred.insert(x);
			assert(x != y);
		}
	}
	return true;
}

void SA::printCmd()
{
	for (int i = 1; i <= orgNodeNum; ++i) {
		if (nodeMap[i] == -1) {
			cout << i << endl;
		}
	}
	return;
}

int SA::getResult()
{
	return nodeMap[0];
}

