#include <mutex>
#include <thread>
#include <fstream>
#include <condition_variable>
#include "rf_pipelines_internals.hpp"

using namespace std;
using namespace rf_pipelines;


static void usage(const char *msg = nullptr)
{
    cerr << "Usage: rfp-time [-rP] [-t NTHREADS] [-j JSON_OUTFILE] file.json [file2.json file3.json ...]\n"
	 << "   -t: change number of worker threads (default 1)\n"
	 << "   -P: don't pin threads to cores (default is to pin threads, this should be done on an otherwise idle machine)\n"
	 << "   -j: write json output from thread 0 to specified file (must not already exist)\n";

    if (msg)
	cerr << "Error: " << msg << "\n";

    exit(2);
}


static void usage(const string &msg)
{
    usage(msg.c_str());
}


// -------------------------------------------------------------------------------------------------
//
// global_context


struct global_context {
    bool tflag = false;   // -t: change number of worker threads (default 1)
    bool Pflag = false;   // -P: don't pin threads to cores (default is to pin threads)
    int nthreads = 1;

    run_params rp;

    int ninputs = 0;
    vector<string> input_filenames;  // length ninputs
    vector<Json::Value> input_json;  // length ninputs

    // One-time barrier, used to start pipelines simultaneously in all worker threads.
    std::mutex barrier_lock;
    std::condition_variable barrier_cv;
    int barrier_count = 0;

    // Length nthreads
    vector<Json::Value> output_json;
    string json_outfile;

    global_context(int argc, char **argv);

    shared_ptr<pipeline> make_pipeline() const;

    void wait_at_barrier();
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
	else if (!strcmp(arg, "-j")) {
	    if (iarg >= argc)
		usage("couldn't parse [-j JSON_OUTFILE] argument");
	    if (json_outfile.size() > 0)
		usage("double [-j JSON_OUTFILE] argument specified");
	    
	    this->json_outfile = argv[iarg];
	    iarg++;

	    if (file_exists(json_outfile))
		usage("json outfile '" + json_outfile + "' already exists");
	}	    
	else {
	    if (arglen == 1)
		usage();

	    for (int j = 1; j < arglen; j++) {
		if (arg[j] == 'P')
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

    this->output_json.resize(nthreads);
}


// Calls bind() but not allocate()
shared_ptr<pipeline> global_context::make_pipeline() const
{
    auto ret = make_shared<pipeline> ();
    for (int i = 0; i < ninputs; i++)
	ret->add(pipeline_object::from_json(input_json[i]));

    ret->bind(this->rp);
    return ret;
}


void global_context::wait_at_barrier()
{
    unique_lock<mutex> l(barrier_lock);
    barrier_count++;
    
    if (barrier_count == nthreads) {
	barrier_cv.notify_all();
	return;
    }

    while (barrier_count < nthreads)
        barrier_cv.wait(l);
}


// -------------------------------------------------------------------------------------------------
//
// Worker thread


static void pin_current_thread_to_core(int core_id)
{
#ifdef __APPLE__
    if (core_id == 0)
	cerr << "warning: pinning threads to cores is not implemented in osx\n";
    return;
#else
    int hwcores = std::thread::hardware_concurrency();
    
    if ((core_id < 0) || (core_id >= hwcores))
	throw runtime_error("pin_thread_to_core: core_id=" + to_string(core_id) + " is out of range (hwcores=" + to_string(hwcores) + ")");

    pthread_t thread = pthread_self();

    cpu_set_t cs;
    CPU_ZERO(&cs);
    CPU_SET(core_id, &cs);

    int err = pthread_setaffinity_np(thread, sizeof(cs), &cs);
    if (err)
        throw runtime_error("pthread_setaffinity_np() failed");
#endif
}


static void worker_thread_main(global_context *c, int thread_id)
{
    if (!c->Pflag)
	pin_current_thread_to_core(thread_id);

    auto p = c->make_pipeline();
    p->allocate();
    
    c->wait_at_barrier();
    c->output_json[thread_id] = p->run(c->rp);
}


// -------------------------------------------------------------------------------------------------
//
// print_timing()


using timing_dict_t = std::unordered_map<std::string, double>;


template<typename T>
static vector<Json::Value> _get(const vector<Json::Value> &v, const T &x)
{
    int n = v.size();
    vector<Json::Value> ret(n);

    for (int i = 0; i < n; i++)
	ret[i] = v[i][x];

    return ret;
}


static double _get_cpu_time(const vector<Json::Value> &v)
{
    double tot = 0.0;
    for (size_t i = 0; i < v.size(); i++)
	tot += double_from_json(v[i], "cpu_time");

    return tot / v.size();
}


static void print_timing(timing_dict_t &d, const vector<Json::Value> &v, string name, int indent_level)
{
    int n = v.size();
    rf_assert(n > 0);

    string class_name = string_from_json(v[0], "class_name");
    double cpu_time = _get_cpu_time(v);

    if (name.size() == 0)
	name = string_from_json(v[0], "name");

    if (!has_key(d, class_name))
	d[class_name] = 0.0;

    d[class_name] += cpu_time;

    cout << string(4*indent_level, ' ');
    cout << "[" << cpu_time << " sec] " << name << endl;

    if (class_name == "pipeline") {
	auto v2 = _get(v, "pipeline");
	int m = v2[0].size();
	for (int j = 0; j < m; j++)
	    print_timing(d, _get(v2,j), "", indent_level+1);
    }
    else if (class_name == "wi_sub_pipeline") {
	auto v2 = _get(v, "pipeline");
	rf_assert(v2[0].size() == 3);

	auto v3 = _get(v2, 1);
	string c3 = string_from_json(v3[0], "class_name");

	print_timing(d, _get(v2,0), "", indent_level+1);

	if (c3 != "pipeline")
	    print_timing(d, v3, "", indent_level+1);
	else {
	    v3 = _get(v3, "pipeline");
	    int m = v3[0].size();
	    for (int j = 0; j < m; j++)
		print_timing(d, _get(v3,j), "", indent_level+1);
	}

	print_timing(d, _get(v2,2), "", indent_level+1);
    }
}


static void print_timing(const global_context &c)
{
    timing_dict_t d;
    auto j = _get(c.output_json, "pipeline");

    ssize_t nsamples_fully_processed = ssize_t_from_json(c.output_json[0], "nsamples_fully_processed");
    ssize_t nsamples_partially_processed = ssize_t_from_json(c.output_json[0], "nsamples_partially_processed");

    cout << "Time spent in each transform:\n";
    for (int i = 0; i < c.ninputs; i++)
	print_timing(d, _get(j,i), c.input_filenames[i], 1);

    vector<string> class_names;
    for (const auto &p: d)
	class_names.push_back(p.first);
    
    std::sort(class_names.begin(), class_names.end());

    cout << "\nTotal time spent in each transform type:\n";
    for (const auto &s: class_names) {
	// Skip container classes
	if ((s == "pipeline") || (s == "wi_sub_pipeline"))
	    continue;
	cout << "    " << s << ": " << d[s] << " sec\n";
    }

    cout << "\n"
	 << "Grand total: " << _get_cpu_time(c.output_json) << " sec\n"
	 << "Number of samples fully processed: " << nsamples_fully_processed << "\n"
	 << "Number of samples partially processed: " << nsamples_partially_processed << "\n";
}


// -------------------------------------------------------------------------------------------------


int main(int argc, char **argv)
{
    global_context c(argc, argv);

    // Construct and bind throwaway pipeline, so that some error checking
    // happens before spawning worker threads.
    c.make_pipeline();

    cout << "Command line:";
    for (int i = 0; i < argc; i++)
	cout << " " << argv[i];
    cout << "\n";

    if (c.nthreads > 1) {
	if (!c.Pflag)
	    cout << "Each timing thread will be pinned to its own core (can be disabled with rfp-time -P)\n";
	cout << "Spawning " << c.nthreads << " timing threads\n";
    }

    cout << endl;

    vector<std::thread> threads;

    for (int i = 0; i < c.nthreads; i++)
	threads.push_back(std::thread(worker_thread_main, &c, i));

    for (int i = 0; i < c.nthreads; i++)
	threads[i].join();

    print_timing(c);

    if (c.json_outfile.size() > 0) {
	ofstream f(c.json_outfile);
	if (f.fail()) {
	    cout << "Couldn't open json outfile '" << c.json_outfile << "'\n";
	    exit(1);
	}

	f << c.output_json[0] << endl;
	f.close();
	cout << "Wrote " << c.json_outfile << endl;
    }

    return 0;
}
