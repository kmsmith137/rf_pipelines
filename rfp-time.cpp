#include <fstream>
#include "rf_pipelines_internals.hpp"

using namespace std;
using namespace rf_pipelines;


// -------------------------------------------------------------------------------------------------


static void usage(const char *msg = nullptr)
{
    cerr << "Usage: rfp-time [-rP] [-t NTHREADS] file.json [file2.json file3.json ...]\n"
	 << "   -r: enable recursive timing of all transforms in pipeline\n"
	 << "   -P: don't pin threads to cores (default is to pin threads, this should be done on an otherwise idle machine)\n";

    if (msg)
	cerr << "Error: " << msg << "\n";

    exit(2);
}


static void usage(const string &msg)
{
    usage(msg.c_str());
}


int main(int argc, char **argv)
{
    bool rflag = false;
    bool tflag = false;
    bool Pflag = false;
    int nthreads = 1;

    vector<char *> json_filenames;
    int iarg = 1;

    while (iarg < argc) {
	char *arg = argv[iarg];
	int arglen = strlen(arg);
	iarg++;

	if (arg[0] != '-')
	    json_filenames.push_back(arg);
	else if (!strcmp(argv[iarg], "-t")) {
	    // Note "iarg++" in this line
	    if ((iarg == argc-1) || !lexical_cast(argv[iarg++], nthreads))
		usage("couldn't parse [-t NTHREADS] argument");
	    if (nthreads < 1)
		usage("invalid 'nthreads': " + to_string(nthreads));
	    if (nthreads > 40)
		usage("very large 'nthreads', presumably unintentional: " + to_string(nthreads));
	    if (tflag)
		usage("double [-t NTHREADS] argument specified");
	    tflag = true;
	}
	else {
	    if (arglen == 1)
		usage();

	    for (int j = 1; j < arglen; j++) {
		if (arg[j] == 'r')
		    rflag = true;
		else if (arg[j] == 'P')
		    Pflag = true;
		else
		    usage("unrecognized flag '-" + string(1,arg[j]) + "'");
	    }
	}

	// advance token
	iarg++;
    }

    int njson = json_filenames.size();
    vector<Json::Value> json_values(njson);
	
    if (njson == 0)
	usage();

    for (size_t i = 0; i < json_filenames.size(); i++) {
	std::ifstream f(json_filenames[i]);
	if (f.fail()) {
	    cerr << json_filenames[i] << ": couldn't open file\n";
	    exit(1);
	}

	f >> json_values[i];
    }

    return 0;
}
