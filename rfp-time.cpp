#include <thread>
#include <fstream>
#include "rf_pipelines_internals.hpp"

using namespace std;
using namespace rf_pipelines;


static void usage(const char *msg = nullptr)
{
    cerr << "Usage: rfp-time [-rP] [-t NTHREADS] file.json [file2.json file3.json ...]\n"
	 << "   -t: change number of worker threads (default 1)\n"
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


// -------------------------------------------------------------------------------------------------


struct global_context {
    bool rflag = false;   // -r: enable recursive timing of all transforms in pipeline
    bool tflag = false;   // -t: change number of worker threads (default 1)
    bool Pflag = false;   // -P: don't pin threads to cores (default is to pin threads)
    int nthreads = 1;

    run_params rp;

    int ninputs = 0;
    vector<string> input_filenames;
    vector<Json::Value> input_json;

    global_context(int argc, char **argv);

    shared_ptr<pipeline_object> make_pipeline() const;
};


global_context::global_context(int argc, char **argv)
{
    this->rp.outdir = "";
    this->rp.verbosity = 0;

    int iarg = 1;

    while (iarg < argc) {
	char *arg = argv[iarg];
	int arglen = strlen(arg);
	iarg++;

	if (arg[0] != '-')
	    this->input_filenames.push_back(arg);
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

	    this->tflag = true;
	}
	else {
	    if (arglen == 1)
		usage();

	    for (int j = 1; j < arglen; j++) {
		if (arg[j] == 'r')
		    this->rflag = true;
		else if (arg[j] == 'P')
		    this->Pflag = true;
		else
		    usage("unrecognized flag '-" + string(1,arg[j]) + "'");
	    }
	}
    }

    this->ninputs = input_filenames.size();
    this->input_json.resize(ninputs);
	
    if (ninputs == 0)
	usage();

    // Read json files
    for (int i = 0; i < ninputs; i++) {
	std::ifstream f(input_filenames[i]);
	if (f.fail()) {
	    cerr << input_filenames[i] << ": couldn't open file\n";
	    exit(1);
	}
	f >> input_json[i];
    }
}


// Calls bind() but not allocate()
shared_ptr<pipeline_object> global_context::make_pipeline() const
{
    shared_ptr<pipeline_object> ret;

    if (ninputs == 1)
	ret = pipeline_object::from_json(input_json[0]);
    else {
	auto p = make_shared<pipeline> ();
	for (int i = 0; i < ninputs; i++)
	    p->add(pipeline_object::from_json(input_json[i]));
	ret = p;
    }

    ret->bind(this->rp);
    return ret;
}


// -------------------------------------------------------------------------------------------------
//
// Note: we use bare pointers and don't worry about memory leaks.


static void worker_thread_main(const global_context &c, Json::Value *json_output)
{
    auto p = c.make_pipeline();

    p->allocate();
    *json_output = p->run(c.rp);
}


// -------------------------------------------------------------------------------------------------



int main(int argc, char **argv)
{
    global_context c(argc, argv);

    // Construct and bind throwaway pipeline, so that some error checking
    // happens before spawning worker threads.
    c.make_pipeline();

    vector<std::thread> threads;
    vector<Json::Value> json_output(c.nthreads);

    for (int i = 0; i < c.nthreads; i++)
	threads.push_back(std::thread(worker_thread_main, c, &json_output[i]));

    for (int i = 0; i < c.nthreads; i++)
	threads[i].join();

    cout << json_output[0];
    return 0;
}
