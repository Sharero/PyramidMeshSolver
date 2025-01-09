#ifndef SLAE_H
#define SLAE_H

#include <set>
#include <tuple>
#include <vector>

#include "msg.h"

class SLAE {
   private:
    std::size_t slae_size{0};

    std::vector<int> ig, jg;
    std::vector<double> di, gg, f, q;

   public:
    void setSlaeSize(std::size_t& new_size) {
        slae_size = new_size;
    }

    [[nodiscard]]
    const std::vector<double>& getResultVector() const {
        return q;
    }

    void setIndexArrays(std::vector<std::set<int>>& connections);

    void allocateMemory(std::size_t size);

    void addLowerTriangleElement(
        const std::tuple<int, int, double>& element_parameters);

    void addDiagonalElement(int diagonal_element_index,
                            double diagonal_element_value);

    void addRightPartElement(int right_part_element_index,
                             double right_part_element_value);

    void applyFirstBoundaryConditions(const std::vector<std::pair<int, double>>&
                                          first_boundary_condition_nodes);

    void solveSLAE();
};
#endif
