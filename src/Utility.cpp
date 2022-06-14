#include "Utility.h"
//静态变量初始化
Random* Random::instance = new Random();

Random* Random::getInstance()
{
    return instance;
}