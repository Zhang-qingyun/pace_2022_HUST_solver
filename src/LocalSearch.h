#ifndef LS_H
#define LS_H

#include "datastuct.h"
#include "Utility.h"
#include "config.h"
using namespace std;

class LS {
private:
	int nodeNum;				
	int edgeNum;				
	vector<Node>& graph;		

	vector<DLinkedNode> curSol;			
	vector<DLinkedNode> curBest;		


	vector<int> unnumbered;		
	vector<int> unbIdx;			
	vector<bool> ok;			
	vector<vector<int>> delta;	
	vector<vector<int>> pos;	
	vector<vector<vector<int>>> rm;	


	double t0;
	int maxMvt;
	double alph;
	int maxFail;


	double t;
	int nbFail;
	int totalIter;


	LL tv;
	LL maxTv;

	Random* rand;
	Config* cfg;
public:
	LS(int _nodeNum, int _edgeNum, vector<Node>& _graph);


	void init();


	void insertAfter(int u, int v);


	void insertBefore(int u, int v);


	void unlink(int u);


	void updateVertex(int v);


	void chooseMove(int& u, int& type);


	void applyMove(int u, int type);


	void subgraphInit();
	void subgraphExe();

	void RTS();


	bool check(vector<DLinkedNode>& lst);


	vector<DLinkedNode>& getCurBest();

	int getTotalIter();

public:

	int remain;					

	vector<int> nodeMap;		
	vector<int> indgr;
	vector<int> outdgr;
	
	LRUGraph ini;
	void iniErase(int u);
	void generateSolution();	
	void constructInitSol();	
	void randomInitSol();		
	void MVCInitSolution(double sRatio);		
	bool reduceFacade(int move);
	void reduceGraph();
	bool in0_out0();
	bool in1_out1();
	bool loop();

	LL hash(int u, int v) {
		return (LL)u * nodeNum + v;
	}
};

#endif