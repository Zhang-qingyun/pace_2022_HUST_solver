#include "Utility.h"
//��̬������ʼ��
Random* Random::instance = new Random();

Random* Random::getInstance()
{
    return instance;
}