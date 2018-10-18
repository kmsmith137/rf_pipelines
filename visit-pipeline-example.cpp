// This toy program shows how to use rf_pipelines::visit_pipeline().
// It prints out a pipeline, by "visiting" each pipeline_object, and printing its name.

#include "rf_pipelines_base_classes.hpp"
#include "rf_pipelines_internals.hpp"

using namespace std;


void print_pipeline(const shared_ptr<rf_pipelines::pipeline_object> &p)
{
}


int main(int argc, char **argv)
{
    if (argc != 2)
	throw runtime_error("usage: test-visitor <chain.json>");

    Json::Value j = rf_pipelines::json_read(argv[1]);
    shared_ptr<rf_pipelines::pipeline_object> p = rf_pipelines::pipeline_object::from_json(j);

    // Define function object which will "visit" each pipeline_object in the pipeline.
    auto visitor = [](const shared_ptr<rf_pipelines::pipeline_object> &p, int depth)
    {
	// Cosmetic: indent to appropriate depth.
	for (int i = 0; i < depth; i++)
	    cout << "    ";
	
	cout << p->name << endl;
    };

    visit_pipeline(visitor, p);
    return 0;
}
