#ifndef ELEMENT_H
#define ELEMENT_H

#include <array>

const int NUMBER_OF_VERTICES_OF_PYRAMID = 5;

class Element {
   private:
    std::array<int, NUMBER_OF_VERTICES_OF_PYRAMID> node_indexes;

   public:
    Element() = default;
    Element(int vertice_index_1, int vertice_index_2, int vertice_index_3,
            int vertice_index_4, int vertice_index_5)
        : node_indexes{vertice_index_1, vertice_index_2, vertice_index_3,
                       vertice_index_4, vertice_index_5} {}

    [[nodiscard]]
    const std::array<int, NUMBER_OF_VERTICES_OF_PYRAMID>& getNodeIndexes()
        const {
        return node_indexes;
    }

    void setNodeIndexes(
        const std::array<int, NUMBER_OF_VERTICES_OF_PYRAMID>& new_indexes) {
        node_indexes = new_indexes;
    }
};

#endif
