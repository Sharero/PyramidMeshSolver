#include <vector>
#include <math.h>

using namespace std;

class MSG {

   vector<int> gi;
   vector<int> gj;

   vector<double> diL;
   vector<double> ggL;
   vector<double> di;
   vector<double> gg;
   vector<double> rp;
   vector<double> x0;
   vector<double> r;
   vector<double> z;
   vector<double> p;
   vector<double> s;

   int n;

   void MakeLLTDecomposition();
   void CalculateLLT(const vector<double>& f, vector<double>& x);
   void CalculateLT(const vector<double>& f, vector<double>& x);
   void CalculateL(const vector<double>& f, vector<double>& x);
   void Ax(const vector<double>& f, vector<double>& x);

   double Scalar(const vector<double>& a, const vector<double>& b);

public:
   void Calculate(vector<double>& solution);
   void Initialize(const vector<int>& gi_s, const vector<int>& gj_s, const vector<double>& di_s, const vector<double>& gg_s, const vector<double>& rp_s, int n_s);
};
