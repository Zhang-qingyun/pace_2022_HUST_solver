#ifndef GLOBAL_H
#define GLOBAL_H

#include <vector>
#include <list>
#include <iterator>
#include <string>
#include <time.h>

using namespace std;

#define SUBMIT true
#define SOURCE 1
#define LL long long
#define PII pair<int, int>

#if SOURCE==1
// LOCAL
#define BASE_DIR_DAT "../instances/"
#define BASE_DIR_SOL "../solutions/"
#define LOG_RUN_FILE ""

#endif

typedef struct {
	int depot;
	vector<int>::iterator position;
	float cost;
} typedef_location;

enum class Enum_Local_Search_Type 
{
	RANDOM,
	SEQUENTIAL,
	NOT_APPLIED
};

typedef struct {
	int index;
	float cost;
} typedef_order;

typedef struct {
	int i; // id, indice, ...
	float x;
	float y;
} typedef_point;


#endif // !GLOBAL_H

#pragma once
