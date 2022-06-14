#ifndef SA_H
#define SA_H

#include "datastuct.h"
#include "Utility.h"

class SA {
public:
	string instString;
	string instCode;

	int orgNodeNum;
	int orgEdgeNum;
	vector<HashNode> org;
	vector<int> nodeMap;		
	int reduceNodeNum;				
	int reduceEdgeNum;				

	vector<Graph> subgraphs;

	Random* rand;
public:
	SA();
	void init();
	bool reduceFacade(int move);
	void reduceGraph();
	bool in0();
	bool out0();
	bool in1();
	bool out1();
	bool loop();
	bool pie();
	bool core();
	bool dome();
	void restoreGraph(int gidx, vector<DLinkedNode>& sol);

	bool isPIV(int u);					
	bool isClique(vector<int>& cand);	

	bool allPIEdgeJudge(int gidx);

	int calConnect();					

	void solve();

	bool processInstanceFiles(char* dataFile, char* instCode);
	bool processContestFiles(char* dataFile, char* instCode);
	bool processCmd();

	void printCmd();

	int getResult();
};

#endif
