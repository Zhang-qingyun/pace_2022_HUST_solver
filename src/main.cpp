#include <iostream>
#include <fstream>
#include <string>
#include <csignal>

#include "config.h"
#include "SA.h"
#include "Utility.h"
using namespace std;

#define debug(x) cout << #x << " : " << x << endl

int submit() {

	Config::getInstance()->init();

	SA problem;

	if (!problem.processCmd()) {
		return EXIT_FAILURE;
	}

	problem.solve();

	problem.printCmd();


	return 0;
}

void handle(int sigNum) {
    Config::getInstance()->setFlag();
}

int main() {
    signal(SIGTERM, handle);

	#if SUBMIT
	submit();
	return 0;
	#endif // SUBMIT


}
