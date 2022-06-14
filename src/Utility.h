#ifndef UTILITY_H
#define UTILITY_H

#include <random>
#include <ctime>
#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include <stack>
#include <algorithm>
#include <climits>

#include "global.h"
using namespace std;

class Random {
private:
	Random() : rgen(0) {}
	static Random* instance;
public:
	static Random* getInstance();

	using Generator = std::mt19937;
	using Generator_64 = std::mt19937_64;


	Random(int seed) : rgen(seed), rgen_64(seed){}


	static int generateSeed() {
		return static_cast<int>(std::time(nullptr) + std::clock());
	}

	Generator::result_type operator()() { return rgen(); }

	// pick with probability of (numerator / denominator).
	bool isPicked(unsigned numerator, unsigned denominator) {
		return ((rgen() % denominator) < numerator);
	}

	// pick from [min, max).
	int pick(int min, int max) {
		return ((rgen() % (max - min)) + min);
	}
	// pick from [0, max).
	int pick(int max) {
		return (rgen() % max);
	}

	LL pick(LL max) {
		return (rgen_64() % max);
	}

	static void randPermutation(int size, int min, vector<int>& elements) {

		int i;

		// inizialize
		for (i = 0; i < size; ++i) {
			elements.push_back(min);
			min++;
		}

		random_shuffle(elements.begin(), elements.end());

	}

	Generator rgen;
	Generator_64 rgen_64;
};

struct UnionFind {
	vector<int> f;
	vector<int> num;

	void init(int n) {
		f.resize(n + 1);
		num.resize(n + 1, 1);
		for (int i = 0; i <= n; ++i) {
			f[i] = i;
		}
	}

	int find(int node) {
		int i, j = node;
		while (f[node] != node)node = f[node];
		while (j != f[j]) {   //Â·¾¶Ñ¹Ëõ
			i = j;
			j = f[j];
			f[i] = node;
		}
		return node;
	}

	void merge(int a, int b) {
		int fa = find(a);
		int fb = find(b);
		if (fa != fb) {
			f[fb] = fa;
			num[fa] += num[fb];
		}
	}
};

static string num2str(int i)
{
	char ss[10];
	sprintf(ss, "%03d", i);
	return ss;
}

#endif // !UTILITY_H