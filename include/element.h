#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>

class Element {
   private:
    std::vector<int> node_indexes;

   public:
    Element() = default;
    Element(std::initializer_list<int> indexes) : node_indexes(indexes) {}

    [[nodiscard]]
    const std::vector<int>& getNodeIndexes() const {
        return node_indexes;
    }

    void setNodeIndexes(const std::vector<int>& new_indexes) {
        node_indexes = new_indexes;
    }
};

#endif
