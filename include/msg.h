#ifndef MSG_H
#define MSG_H

#include <vector>

class MSG {
   private:
    std::vector<int> gi, gj;
    std::vector<double> di, gg;

    std::vector<double> di_of_lower_triangle, gg_of_lower_triangle;
    std::vector<double> rp, x0, r, z, p, s;

    std::size_t slae_size{0};

    void makeLLTDecomposition();

    void calculateLLT(const std::vector<double>& f_vector,
                      std::vector<double>& x_vector);

    void calculateLT(const std::vector<double>& f_vector,
                     std::vector<double>& x_vector);

    void calculateL(const std::vector<double>& f_vector,
                    std::vector<double>& x_vector);

    void multiplyMatrixAndVector(const std::vector<double>& f_vector,
                                 std::vector<double>& x_vector);

    [[nodiscard]] double calculateScalarProduct(
        const std::vector<double>& first_vector,
        const std::vector<double>& second_vector) const;

   public:
    MSG(const std::vector<int>& rows_s, const std::vector<double>& di_s,
        const std::vector<int>& columns_s, const std::vector<double>& gg_s,
        std::size_t n_s, const std::vector<double>& rp_s)
        : gi(rows_s),
          gj(columns_s),
          di(di_s),
          gg(gg_s),
          rp(rp_s),
          slae_size(n_s) {
        x0.resize(slae_size, 0);
        r.resize(slae_size, 0);
        z.resize(slae_size, 0);
        p.resize(slae_size, 0);
        s.resize(slae_size, 0);

        di_of_lower_triangle.resize(di_s.size(), 0);
        gg_of_lower_triangle.resize(gg_s.size(), 0);

#pragma unroll 4
        for (int i = 0; i < di_s.size(); i++) {
            di_of_lower_triangle[i] = di[i];
        }

#pragma unroll 4
        for (int i = 0; i < gg_s.size(); i++) {
            gg_of_lower_triangle[i] = gg[i];
        }
    }

    void calculateSLAE(std::vector<double>& solution);
};
#endif
