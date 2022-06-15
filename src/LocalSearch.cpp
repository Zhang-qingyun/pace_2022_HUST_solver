#include <fstream>
#include <unordered_set>
#include <queue>
#include <sstream>
#include <cassert>

#include "LocalSearch.h"
#include "maxVerVC.h"
#include "minVerVC.h"

LS::LS(int _nodeNum, int _edgeNum, vector<Node>& _graph) : nodeNum(_nodeNum), edgeNum(_edgeNum), graph(_graph) {}

void LS::init()
{
	ok.resize(nodeNum + 1, false);
	delta.resize(nodeNum + 1, vector<int>(2, 0));
	pos.resize(nodeNum + 1, vector<int>(2));
	unnumbered.resize(nodeNum + 1);
	unbIdx.resize(nodeNum + 1);
	rm.resize(nodeNum + 1, vector<vector<int>>(2));

	curSol.resize(nodeNum + 1);
	curBest.resize(nodeNum + 1);

	rand = Random::getInstance();
	cfg = Config::getInstance();

	t0 = 1;
	maxMvt = 5 * nodeNum;
	maxFail = 20;
	alph = 0.9999;

	tv = 5;
	nbFail = 0;
}

void LS::insertAfter(int u, int v)
{
	curSol[u].prev = v;
	curSol[u].next = curSol[v].next;
	curSol[v].next = u;
	curSol[curSol[u].next].prev = u;
}

void LS::insertBefore(int u, int v)
{
	curSol[u].prev = curSol[v].prev;
	curSol[u].next = v;
	curSol[v].prev = u;
	curSol[curSol[u].prev].next = u;
}

void LS::unlink(int u)
{
	curSol[curSol[u].next].prev = curSol[u].prev;
	curSol[curSol[u].prev].next = curSol[u].next;
	curSol[u].idx = -1;
}

void LS::updateVertex(int v)
{

	int idx_prev = 0;
	int idx_next = 0;
	LL i_prev = 0;
	LL i_next = INT64_MAX;
	for (auto prev : graph[v].pred) {
		if (curSol[prev].idx != -1) {
			if (curSol[prev].idx > i_prev) {
				i_prev = curSol[prev].idx;
				idx_prev = prev;
			}
		}
	}
	for (auto next : graph[v].succ) {
		if (curSol[next].idx != -1) {
			if (curSol[next].idx < i_next) {
				i_next = curSol[next].idx;
				idx_next = next;
			}
		}
	}

	for (int i = 0; i < 2; ++i) {
		rm[v][i].clear();
	}
	int cv_prev = 0;
	int cv_next = 0;
	for (auto next : graph[v].succ) {
		if (curSol[next].idx != -1 && curSol[next].idx <= i_prev) {
			cv_prev++;
			rm[v][0].push_back(next);
		}
	}

	for (auto prev : graph[v].pred) {
		if (curSol[prev].idx >= i_next) {
			cv_next++;
			rm[v][1].push_back(prev);
		}
	}

	ok[v] = true;
	pos[v][0] = idx_prev;
	pos[v][1] = idx_next;
	delta[v][0] = -1 + cv_prev;
	delta[v][1] = -1 + cv_next;
}

void LS::chooseMove(int& u, int& type)
{
	int idx = rand->pick(curSol[0].idx);
	u = unnumbered[idx];
	type = rand->pick(2);
	updateVertex(u);
}

void LS::applyMove(int u, int type)
{

	int num = unbIdx[u];
	unnumbered[num] = unnumbered[curSol[0].idx - 1];
	unbIdx[unnumbered[num]] = num;
	curSol[0].idx--;


	int v = pos[u][type];
	if (type == 0) {
		insertAfter(u, v);
	}
	else {
		insertBefore(u, v);
	}


	for (auto node : rm[u][type]) {
		unlink(node);
		unnumbered[curSol[0].idx] = node;
		unbIdx[node] = curSol[0].idx++;
	}

	int prev = curSol[u].prev;
	int next = curSol[u].next;
	if (next == 0) {
		curSol[u].idx = curSol[prev].idx + INT_MAX;
	}
	else if (curSol[next].idx > curSol[prev].idx + 1) {
		curSol[u].idx = curSol[prev].idx + (curSol[next].idx - curSol[prev].idx) / 2;
	}
	else {
		// 重新排序
		int idx = curSol[0].next;
		LL base = INT_MAX;
		base <<= 10;
		LL cnt = base;
		while (idx != 0) {
			curSol[idx].idx = cnt;
			cnt += base;
			idx = curSol[idx].next;
		}
		//cout << totalIter << endl;
	}
}

bool LS::check(vector<DLinkedNode>& lst)
{
	unordered_set<int> us;
	vector<int> indegree(nodeNum + 1, 0);
	int idx = lst[0].next;
	while (idx != 0) {
		us.insert(idx);
		idx = curSol[idx].next;
	}
	queue<int> que;
	for (int i = 1; i <= nodeNum; ++i) {
		if (!us.count(i)) continue;
		for (auto adj : graph[i].pred) {
			if (us.count(adj)) {
				indegree[i]++;
			}
		}
		if (indegree[i] == 0) {
			que.push(i);
		}
	}
	int cnt = 0;
	while (!que.empty()) {
		int cur = que.front();
		que.pop();
		cnt++;
		for (auto adj : graph[cur].succ) {
			if (!us.count(adj)) continue;
			if (--indegree[adj] == 0) {
				que.push(adj);
			}
		}
	}
	return cnt != us.size();
}

void LS::subgraphInit()
{
	init();
	//return;
	for (int i = 0; i < nodeNum; ++i) {
		unnumbered[i] = i + 1;
		unbIdx[i + 1] = i;
	}
	curSol[0].idx = nodeNum;
	curBest[0].idx = nodeNum;
	t = t0;
	nbFail = 0;
	totalIter = 0;


	constructInitSol();
}

void LS::RTS()
{
	if (cfg->TLE()) return;
	if (curSol[0].idx > 100000) maxTv = 18;
	else if (curSol[0].idx > 10000) maxTv = 15;
	else if (curSol[0].idx > 1000) maxTv = 12;
	else if (curSol[0].idx > 100) maxTv = 9;
	else maxTv = 5;
	// t 0.1 -> 0.01
	// TBE
	LL tvT = tv + curBest[0].idx;
//	printf("cur:%lld t:%f tvT:%lld\n", curBest[0].idx, t, tvT);
	int nbMvtTBE = 0;
	int maxMvtTBE = 5;
	int nbMvt = 0;
	while (nbMvt < maxMvtTBE && !cfg->TLE()) {
		int len = curSol[0].idx;
		while (len > 0 && !cfg->TLE()) {
		// 选择动作
			int idx = rand->pick(len--);
			int u = unnumbered[idx];
			int type = rand->pick(2);
			updateVertex(u);

			if (delta[u][type] + curSol[0].idx <= tvT) {
				applyMove(u, type);
				if (curSol[0].idx < curBest[0].idx) {
					curBest = curSol;
					//printf("TBE iter:%d f:%lld tvT:%lld t:%f time:%d\n", totalIter, curBest[0].idx, tvT, t, Config::getInstance()->getTime());
				}
			}
			else if (delta[u][1 - type] + curSol[0].idx <= tvT) {
				applyMove(u, 1 - type);
				if (curSol[0].idx < curBest[0].idx) {
					curBest = curSol;
					//printf("TBE iter:%d f:%lld tvT:%lld t:%f time:%d\n", totalIter, curBest[0].idx, tvT, t, Config::getInstance()->getTime());
				}
			}
			totalIter++;
		}
		nbMvt += 1;
	}
	// DBI
	nbMvt = 0;
	int maxMvtDBI = 3 * curSol[0].idx;
	while (nbMvt < maxMvtDBI && !cfg->TLE()) {
		// 选择动作
		int u, type;
		chooseMove(u, type);

		if (delta[u][type] <= 0) {
			applyMove(u, type);

			if (curSol[0].idx < curBest[0].idx) {
				curBest = curSol;
				nbMvt = 0;
				nbFail = 0;
				tv = 1;
				//printf("DBI iter:%d f:%lld tvT:%lld t:%f time:%d\n", totalIter, curBest[0].idx, tvT, t, Config::getInstance()->getTime());
			}
		}
		nbMvt += 1;
		totalIter++;
	}
	nbFail++;
	if (nbFail == maxFail) {
		tv = min(tv + 1, maxTv);
		nbFail = 0;
	}
}

void LS::subgraphExe()
{
	RTS();
	return;
	if (cfg->TLE()) return;
	int nbMvt = 0;
	bool failure = true;
	while (nbMvt < maxMvt && !cfg->TLE()) {
		// 选择动作
		int u, type;
		chooseMove(u, type);

		if (delta[u][type] <= 0 || exp(-delta[u][type] / t) >= rand->pick(10000) * 1.0 / 10000) {
			applyMove(u, type);

			if (curSol[0].idx < curBest[0].idx) {
				failure = false;
				curBest = curSol;
				printf("iter:%d f:%lld time:%d\n", totalIter, curBest[0].idx, Config::getInstance()->getTime());
			}

			nbMvt += 1;
		}
	}
	if (failure == true) nbFail++;
	else nbFail = 0;
	t = t * alph;

	t = max(t, 0.15);

	totalIter += nbMvt;
	//if (curBest[0].idx <= BESTF) break;
}

vector<DLinkedNode>& LS::getCurBest()
{
	return this->curBest;
}

int LS::getTotalIter()
{
	return totalIter;
}

void LS::constructInitSol()
{
	ini.resize(nodeNum);
	nodeMap.resize(nodeNum + 1, 1);
	indgr.resize(nodeNum + 1, 0);
	outdgr.resize(nodeNum + 1, 0);
	remain = nodeNum;

	const static double singleRatio = 0.3;
	int singleEdgeNum = 0;

	for (int i = 1; i <= nodeNum; ++i) {
		for (auto adj : graph[i].pred) {
			ini[i].pred.insert(adj);
		}
		for (auto adj : graph[i].succ) {
			ini[i].succ.insert(adj);
		}

		indgr[i] = ini[i].pred.size();
		outdgr[i] = ini[i].succ.size();
		// 计算单边数量
		for (auto adj : graph[i].succ) {
			if (!ini[i].pred.count(adj)) {
				singleEdgeNum++;
			}
		}
	}

	double sRatio = singleEdgeNum * 1.0 / edgeNum;
	if (sRatio > singleRatio) {
		randomInitSol();
	}
	else {
		MVCInitSolution(sRatio);
	}

	generateSolution();
}

void LS::randomInitSol()
{
	auto cmp = [](const pair<int, int>& p1, const pair<int, int>& p2) {
		return p1.second > p2.second;
	};

	int len = max(nodeNum / 10000, 1) * 2;
	priority_queue<PII, vector<PII>, decltype(cmp)> pq(cmp);

	// 循环构造，移除节点
	int cycle = 0;
	while (remain > 0) {
		int cnt = 1;
		for (int i = 0; i < len; ++i) {
			pq.emplace(-1, -1);
		}
		for (int i : ini.nodes) {
			assert(nodeMap[i] == 1);
			assert(indgr[i] == ini[i].pred.size());
			assert(outdgr[i] == ini[i].succ.size());
			int curf = indgr[i] * outdgr[i];

			if (curf > pq.top().second) {
				pq.emplace(i, curf);
				pq.pop();
				cnt = 1;
			}
			else if (curf == pq.top().second) {
				cnt++;
				int d = rand->pick(cnt);
				if (d < 1) {
					pq.emplace(i, curf);
					pq.pop();
				}
			}
		}

		while (!pq.empty()) {
			auto it = pq.top();
			pq.pop();
			int ru = it.first;
			if (ru == -1) continue;
			assert(nodeMap[ru] == 1);
			nodeMap[ru] = -1;
			remain--;
			for (auto adj : ini[ru].succ.nodes) {
				ini[adj].pred.erase(ru);
				indgr[adj]--;
			}
			for (auto adj : ini[ru].pred.nodes) {
				ini[adj].succ.erase(ru);
				outdgr[adj]--;
			}
			ini.erase(ru);
		}
		reduceGraph();
		cycle++;
	}

	//loop();
	#if !SUBMIT
	//printf("cycle: %d ", cycle);
	#endif // SUBMIT
}

void LS::MVCInitSolution(double sRatio)
{
	int InitRunTime = 0;
	if (sRatio == 0.0) {
		InitRunTime = Config::getInstance()->InitTime() * 0.4 * 1000;
	}
	else if (sRatio < 0.05) {
		InitRunTime = Config::getInstance()->InitTime() * 0.2 * 1000;
	}
	else if (sRatio < 0.1) {
		if (nodeNum < 40000) {
			InitRunTime = 20000;
		}
		else {
			InitRunTime = 30000;
		}
	}
	else {
		if (nodeNum < 40000) {
			InitRunTime = 10000;
		}
		else {
			InitRunTime = 20000;
		}
	}

	InitRunTime = nodeNum * 1.0 / Config::getInstance()->getNodeNum() * InitRunTime;

	#if !SUBMIT
	cout << "MVC:" << InitRunTime << endl;
	#endif // SUBMIT

	vector<int> result;
	if (nodeNum < 40000) {
		vector<minVerVC::Edge> MVCEdge;
		for (int i = 1; i <= nodeNum; ++i) {
			for (auto adj : graph[i].succ) {
				if (i >= adj) continue;
				MVCEdge.push_back({ i - 1, adj - 1 });
			}
		}
		
		minVerVC m(0, std::chrono::milliseconds(InitRunTime));
		bool flag = m.NuMVC_minVertexCovSolver(MVCEdge, nodeNum, result);
	}
	else {
		vector<maxVerVC::Edge> MVCEdge;
		for (int i = 1; i <= nodeNum; ++i) {
			for (auto adj : graph[i].succ) {
				if (i >= adj) continue;
				MVCEdge.push_back({ i - 1, adj - 1 });
			}
		}
		maxVerVC f(0, std::chrono::milliseconds(InitRunTime));
		bool flag = f.FastVC_minVertexCovSolver(MVCEdge, nodeNum, result);
	}

	for (auto node : result) {
		nodeMap[node + 1] = -1;
	}
}

void LS::generateSolution()
{
	vector<int> indegree(nodeNum + 1, 0);
	deque<int> dq;
	int cnt = 0;
	for (int i = 1; i <= nodeNum; ++i) {
		if (nodeMap[i] == -1) {
			unnumbered[cnt] = i;
			unbIdx[i] = cnt;
			cnt++;
		}
		else {
			for (auto adj : graph[i].pred) {
				if (nodeMap[adj] == -1) continue;
				indegree[i]++;
			}
			if (indegree[i] == 0) {
				dq.push_back(i);
			}
		}
	}

	curSol[0].idx = cnt;

	int test = 0;
	LL idx = INT_MAX;
	while (!dq.empty()) {
		int cur = dq.front();
		dq.pop_front();
		insertBefore(cur, 0);
		curSol[cur].idx = idx;
		idx += INT_MAX;
		test++;
		for (auto adj : graph[cur].succ) {
			if (nodeMap[adj] == -1) continue;
			if (--indegree[adj] == 0) {
				dq.push_back(adj);
			}
		}
	}

	curBest = curSol;
	assert(test + cnt == nodeNum);
}

bool LS::reduceFacade(int move)
{
	bool result = false;

	switch (move) {
	case 1:
		result = in0_out0();
		break;
	case 2:
		result = in1_out1();
		break;
	case 3:
		result = loop();
		break;
	}

	return result;
}

void LS::reduceGraph()
{
	int MAXMOVE = 2;
	
	for (int move = 1; move <= MAXMOVE; ++move) {
		reduceFacade(move);
	}
}

bool LS::in0_out0()
{
	bool flag = false;
	vector<int> rm;
	for (int i : ini.nodes) {
		assert(nodeMap[i] == 1);
		assert(ini[i].pred.size() == indgr[i]);
		assert(ini[i].succ.size() == outdgr[i]);
		if (indgr[i] == 0) {
			nodeMap[i] = 0;
			rm.push_back(i);
			remain--;
			for (auto adj : ini[i].succ.nodes) {
				ini[adj].pred.erase(i);
				indgr[adj]--;
			}
			flag = true;
		}
		else if (outdgr[i] == 0) {
			nodeMap[i] = 0;
			rm.push_back(i);
			remain--;
			for (auto adj : ini[i].pred.nodes) {
				ini[adj].succ.erase(i);
				outdgr[adj]--;
			}
			flag = true;
		}
	}
	for (auto node : rm) {
		ini.erase(node);
	}
	return flag;
}

bool LS::in1_out1()
{
	bool flag = false;
	vector<int> rm;
	unordered_set<int> loop;
	for (int i : ini.nodes) {
		assert(nodeMap[i] == 1);
		assert(ini[i].pred.size() == indgr[i]);
		assert(ini[i].succ.size() == outdgr[i]);
		if (indgr[i] == 1) {
			int u = ini[i].pred.nodes[0];
			if (u == i) continue;
			nodeMap[i] = 0;
			rm.push_back(i);
			remain--;
			ini[u].succ.erase(i);
			outdgr[u]--;
			for (auto adj : ini[i].succ.nodes) {
				ini[adj].pred.erase(i);
				indgr[adj]--;
				if (ini[u].succ.insert(adj)) {
					outdgr[u]++;
				}
				if (ini[adj].pred.insert(u)) {
					indgr[adj]++;
				}
				if (adj == u) {
					loop.insert(u);
				}
			}
			flag = true;
		}
		else if (outdgr[i] == 1) {
			int u = ini[i].succ.nodes[0];
			if (u == i) continue;
			nodeMap[i] = 0;
			rm.push_back(i);
			remain--;
			ini[u].pred.erase(i);
			indgr[u]--;
			for (auto adj : ini[i].pred.nodes) {
				ini[adj].succ.erase(i);
				outdgr[adj]--;
				if (ini[u].pred.insert(adj)) {
					indgr[u]++;
				}
				if (ini[adj].succ.insert(u)) {
					outdgr[adj]++;
				}
				if (adj == u) {
					loop.insert(u);
				}
			}
			flag = true;
		}
	}
	for (auto node : rm) {
		ini.erase(node);
	}

	for (auto ru : loop) {
		nodeMap[ru] = -1;
		remain--;
		for (auto adj : ini[ru].succ.nodes) {
			ini[adj].pred.erase(ru);
			indgr[adj]--;
		}
		for (auto adj : ini[ru].pred.nodes) {
			ini[adj].succ.erase(ru);
			outdgr[adj]--;
		}
		ini.erase(ru);
	}
	return flag;
}

bool LS::loop()
{
	bool flag = false;
	vector<int> rm;
	for (int i : ini.nodes) {
		assert(nodeMap[i] == 1);
		if (ini[i].succ.count(i)) {
			nodeMap[i] = -1;
			rm.push_back(i);
			remain--;
			flag = true;
			for (auto adj : ini[i].pred.nodes) {
				ini[adj].succ.erase(i);
			}
			for (auto adj : ini[i].succ.nodes) {
				ini[adj].pred.erase(i);
			}
		}
	}
	for (auto node : rm) {
		ini.erase(node);
	}
	return flag;
}
