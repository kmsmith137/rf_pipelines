// Right now, all that's here is test_median(), but I may add more tests later!

#include "rf_pipelines_internals.hpp"

using namespace std;
using namespace rf_pipelines;


static float slow_median(const vector<float> &v)
{
    vector<float> w = v;
    std::sort(w.begin(), w.end());

    int n = w.size();
    return (n % 2) ? w[n/2] : (0.5 * (w[n/2-1] + w[n/2]));
}


static void test_median(std::mt19937 &rng)
{
    for (int iouter = 0; iouter < 10000; iouter++) {
	int n = randint(rng, 1, 10);
	vector<float> v = uniform_randvec(rng, n, 0.0, 1.0);

	float mslow = slow_median(v);
	float mfast = median(v);

	rf_assert(fabs(mslow-mfast) < 1.0e-6);
    }

    cout << "test_median: pass\n";
}


int main(int argc, char **argv)
{
    std::random_device rd;
    std::mt19937 rng(rd());

    test_median(rng);
    return 0;
}
