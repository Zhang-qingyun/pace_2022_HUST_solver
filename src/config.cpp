#include "config.h"


Config* Config::instance = new Config();

Config::Config()
{
    MAXTIME = 100;
    flag = false;
}

Config* Config::getInstance()
{
    return instance;
}

void Config::init()
{
    start_time = clock();
}

int Config::getTime()
{
    return (clock() - start_time) / CLOCKS_PER_SEC;
}

int Config::InitTime()
{
    return 500;
}

void Config::setNodeNum(int _nodeNum)
{
    nodeNum = _nodeNum;
}

int Config::getNodeNum()
{
    return nodeNum;
}

bool Config::TLE()
{
//    return (clock() - start_time) / CLOCKS_PER_SEC >= MAXTIME;
    return flag;
}

Enum_Local_Search_Type Config::getLocalSearchType() const
{
    return localSearchType;
}

void Config::setFlag() {
    flag = true;
}
