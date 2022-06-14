#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>

#include "global.h"

using namespace std;

class Config {
private:
	Enum_Local_Search_Type localSearchType = Enum_Local_Search_Type::RANDOM;   

	Config();
	static Config* instance;

	clock_t start_time;
	int MAXTIME;
    bool flag;
    int nodeNum;
public:
	static Config* getInstance();

	void init();

	int getTime();

    int InitTime();

    void setFlag();

    void setNodeNum(int _nodeNum);

    int getNodeNum();

    bool TLE();

	Enum_Local_Search_Type getLocalSearchType() const;

};

#endif // !CONFIG_H


