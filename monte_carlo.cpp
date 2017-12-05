#include <random>
#include <iostream>

using namespace std;

int main()
{
   std::random_device rd;
   std::mt19937 rng(rd());
   std::normal_distribution<float> dist;
   dist(rng);
}
