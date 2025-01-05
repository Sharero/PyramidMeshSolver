#include "msg.h"

void MSG::Initialize(const vector<int>& gi_s, const vector<int>& gj_s, const vector<double>& di_s, const vector<double>& gg_s, const vector<double>& rp_s, int n_s) {

   gi = gi_s;
   gj = gj_s;
   di = di_s;
   gg = gg_s;
   rp = rp_s;
   n = n_s;

   x0.resize(n);
   r.resize(n);
   z.resize(n);
   p.resize(n);
   s.resize(n);

   diL.resize(n);
   ggL.resize(gi[n]);

   for (int i = 0; i < n; i++) {

      diL[i] = di[i];
      x0[i] = 0;
   }

   for (int i = 0; i < gi[n]; i++)
      ggL[i] = gg[i];
}

void MSG::Ax(const vector<double>& f, vector<double>& x) {

   for (int i = 0; i < n; i++) {

      x[i] = di[i] * f[i];

      for (int k = gi[i], k1 = gi[i + 1]; k < k1; k++) {

         int j = gj[k];
         x[i] += gg[k] * f[j];
         x[j] += gg[k] * f[i];
      }
   }
}

void MSG::CalculateLLT(const vector<double>& f, vector<double>& x) {

   CalculateL(f, x);
   CalculateLT(x, x);
}

void MSG::CalculateLT(const vector<double>& f, vector<double>& x) {

   vector<double> temp = f;

   for (int k = n, k1 = n - 1; k > 0; k--, k1--) {

      x[k1] = temp[k1] / diL[k1];

      for (int i = gi[k1]; i < gi[k]; i++)
         temp[gj[i]] -= ggL[i] * x[k1];
   }
}

void MSG::CalculateL(const vector<double>& f, vector<double>& x) {

   for (int k = 1, k1 = 0; k <= n; k++, k1++) {

      double sum = 0;

      for (int i = gi[k1]; i < gi[k]; i++)
         sum += ggL[i] * x[gj[i]];

      x[k1] = (f[k1] - sum) / diL[k1];
   }
}

void MSG::MakeLLTDecomposition() {

   double sumD;
   double sumL;

   for (int k = 0; k < n; k++) {

      sumD = 0;

      int iStartIndex = gi[k];
      int iEndIndex = gi[k + 1];

      for (int i = iStartIndex; i < iEndIndex; i++) {

         sumL = 0;

         int jStartIndex = gi[gj[i]];
         int jEndIndex = gi[gj[i] + 1];

         for (int m = iStartIndex; m < i; m++) {

            for (int j = jStartIndex; j < jEndIndex; j++) {

               if (gj[m] == gj[j]) {

                  sumL += ggL[m] * ggL[j];
                  jStartIndex++;
               }
            }
         }

         ggL[i] = (ggL[i] - sumL) / diL[gj[i]];

         sumD += ggL[i] * ggL[i];
      }

      diL[k] = sqrt(diL[k] - sumD);
   }
}

void MSG::Calculate(vector<double>& solution) {

   int maximumIterations = 1000;
   
   bool end = false;

   double epsilon = 1E-15;
   double discrepansy;
   double rpNorm;
   double scalar1;
   double scalar2;
   double betta;
   double alpha;

   Ax(x0, r);

   for (int i = 0; i < n; i++)
      r[i] = rp[i] - r[i];

   MakeLLTDecomposition();
   CalculateLLT(r, z);

   for (int i = 0; i < n; i++)
      p[i] = z[i];

   rpNorm = sqrt(Scalar(rp, rp));
   scalar1 = Scalar(p, r);

   for (int iter = 0; iter < maximumIterations && !end; iter++) {

      discrepansy = sqrt(Scalar(r, r));

      if (epsilon < discrepansy / rpNorm) {

         Ax(z, s);

         alpha = scalar1 / Scalar(s, z);

         for (int i = 0; i < n; i++) {

            x0[i] += alpha * z[i];
            r[i] -= alpha * s[i];
         }

         CalculateLLT(r, p);
         scalar2 = Scalar(p, r);

         betta = scalar2 / scalar1;

         scalar1 = scalar2;

         for (int i = 0; i < n; i++)
            z[i] = p[i] + betta * z[i];
      }
      else
         end = true;
   }

   for (int i = 0; i < n; i++)
      solution[i] = x0[i];
}

double MSG::Scalar(const vector<double>& a, const vector<double>& b) {

   double scalar = 0;

   for (int i = 0; i < n; i++)
      scalar += a[i] * b[i];

   return scalar;
}