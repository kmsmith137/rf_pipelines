#include <thread>
#include <fstream>
#include "rf_pipelines_internals.hpp"

using namespace std;
using namespace rf_pipelines;


// -------------------------------------------------------------------------------------------------


static run_params make_run_params()
{
    run_params rp;
    rp.outdir = "";
    rp.verbosity = 0;

    return rp;
}


static shared_ptr<pipeline_object> make_pipeline(const vector<Json::Value> &json_v)
{
    rf_assert(json_v.size() > 0);
    
    if (json_v.size() == 1)
	return pipeline_object::from_json(json_v[0]);

    auto p = make_shared<pipeline> ();
    
    for (size_t i = 0; i < json_v.size(); i++)
	p->add(pipeline_object::from_json(json_v[i]));

    return p;
}


// -------------------------------------------------------------------------------------------------
//
// Note: we use bare pointers and don't worry about memory leaks.


static void worker_thread_main(const vector<Json::Value> &json_v, Json::Value *json_output)
{
    auto rp = make_run_params();
    auto p = make_pipeline(json_v);

    p->bind(rp);
    p->allocate();
    *json_output = p->run(rp);
}


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
	else if (!strcmp(arg, "-t")) {
	    if (iarg >= argc)
		usage("couldn't parse [-t NTHREADS] argument");

	    char *arg_t = argv[iarg];
	    iarg++;

	    if (!lexical_cast(arg_t, nthreads))
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
    }

    int njson = json_filenames.size();
    vector<Json::Value> json_v(njson);
	
    if (njson == 0)
	usage();

    // Read json files
    for (int i = 0; i < njson; i++) {
	std::ifstream f(json_filenames[i]);
	if (f.fail()) {
	    cerr << json_filenames[i] << ": couldn't open file\n";
	    exit(1);
	}
	f >> json_v[i];
    }

    // Construct and bind throwaway pipeline, so that some error checking
    // happens before spawning worker threads.
    auto p = make_pipeline(json_v);
    p->bind(make_run_params());
    p.reset();

    vector<std::thread> threads;
    vector<Json::Value> json_output(nthreads);

    for (int i = 0; i < nthreads; i++)
	threads.push_back(std::thread(worker_thread_main, json_v, &json_output[i]));

    for (int i = 0; i < nthreads; i++)
	threads[i].join();

    cout << json_output[0];
    return 0;
}
