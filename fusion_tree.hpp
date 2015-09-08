#ifndef FUSION_TREE_H
#define FUSION_TREE_H

#include <bitset>
#include <cmath>
#include <array>
#include <vector>

/*
	KeyLength == log(u) where u is |U| - universe size


*/
template< size_t KeyLength, class Value >
class fusion_tree_static
{
private:
	template< size_t Size >
	using bitset = std::bitset<Size>;

	template< class T, size_t Size >
	using array = std::array<T, Size>;

	template< class T >
	using vector = std::vector<T>;
public:
	using Key = bitset<KeyLength>;



	iterator at(const Key& key) {
		Node* node = m_root->findBucket(key);

	}

private:
	class Node {
	public:
		static constexpr size_t capacity;

		virtual ~Node() = default;

	protected:
		Key m_key;
	};

	class Node_Regular;
	class Node_Leaf;


	Node_Regular* m_root;
};

template< size_t KeyLength, class Value >
class fusion_tree_static<KeyLength, Value>::Node_Regular : public Node {
public:
	Node_Regular() :
		m_children(),
		m_sketch()
	{
		m_children.fill(nullptr);
	}
	
	virtual ~Node_Regular() override = default;

	virtual bool isLeaf() const override {
		return false; 
	}

private:
	using Sketch = bitset<KeyLength>;
	using Children = array<Node*, capacity>;

	Sketch m_sketch;

	/* Child nodes. 
	 * Keep it sorted.
	*/
	Children m_children;
	size_t m_childrenCount;


	Node* findBucket(const Key& key) {
		if (m_childrenCount == 0) {
			return nullptr;
		}

		Key sketch;
		computeKeySketch(key, sketch);

		for (size_t i = 0; i < m_childrenCount; ++i) {

		}
	}

	/* Computes a node sketch.
	 * Have to be called after every node modification.
	 * Complexity: O(capacity ^ 4)
	*/
	void computeNodeSketch() {
		m_sketch.reset();

		if (m_childrenCount < 2) {
			return;
		}
		
		for (size_t nodeIdx = 0; nodeIdx < capacity; ++nodeIdx) {
			const Node* node = m_children[nodeIdx];
			if (node == nullptr) {
				continue;
			}

			// Compute child's sketch
			Key sketch;
			computeKeySketch(node->m_key, sketch);

			m_sketch |= sketch;
		}
	}

	void computeKeySketch(const Key& key, Key& sketch) {
		Key branchingPoints;
		const size_t branchingPointsCount = determineBranchingPoints(key, branchingPoints);

		vector<size_t> branchingPointsPositions;
		for (size_t i = 0; i < branchingPointsCount; ++i) {
			if (branchingPoints[point] == true) {
				branchingPointsPositions.push_back(i);
			}
		}

		Key sketch = branchingPoints & node->m_key;

		// Create magic number m
		vector<bool> mParts((size_t)std::pow(branchingPointsCount, 3));
		mParts[KeyLength - branchingPointsPositions[0]] = true;
		size_t partsCount = 1;

		for (size_t i = 1; i < branchingPointsCount; ++i) {
			size_t j = 0;
			size_t k = 0;
			size_t l = 0;

			while (j != branchingPointsCount) {
				if (k == branchingPointsCount) { ++j; k = 0; l = 0; continue; }
				if (l == partsCount) { ++k; l = 0; continue; }

				size_t number = branchingPointsPositions[i] + mParts[l] - branchingPointsPositions[j];
				if (mParts[number] == false) {
					mParts[number] = true;
					++partsCount;
				}
					
				++l;
			}
				
			for (auto it = mParts.begin(), iend = mParts.end(); it != iend; ++it) {
				if (*it == false) {
					*it = true;
					break;
				}
			}
		}
			
		unsigned long long m = 0;
			
		size_t minPower = (size_t)(-1u);

		for (size_t i = 0, iend = mParts.size(), id = 0; i != iend; ++i) {
			if (mParts[i] == true) {
				const size_t m_ith = branchingPointsPositions[id] + (size_t)std::pow(mParts.size(), i);
				if (m_ith < minPower) {
					minPower = m_ith;
				}
				m += 1ull << (i + m_ith);
				++id;
			}
		}

		unsigned long long product = static_cast<unsigned long long>(sketch) * m;
		product >>= minPower;

		sketch = product;
		sketch.set(mParts.size());
		sketch <<= nodeIdx;
	}

	size_t determineBranchingPoints(const Key& nodeKey, Key& branchingPoints) {
		const size_t key = static_cast<size_t>(nodeKey);
		size_t pos = capacity / 2;
		size_t count = 0;
		for (size_t bit_idx = 1; bit_idx < KeyLength; ++bit_idx) {
			size_t lb, hb; // Search boundary

			if (key < pos) { // We're in the left subtree, search in right
				lb = pos;
				hb = pos + capacity >> bit_idx;
			} else { // We're in the right subtree, search in the left
				lb = pos - capacity >> bit_idx;
				hb = pos;
			}
			while ((lb < hb) && (m_children[lb] == nullptr)) {
				++lb;
			}

			// Determine if there are a branching point
			branchingPoints.set(bit_idx - 1, lb != hb);
			count += lb != hb;

			if (key < pos) {
				pos /= 2;
			} else {
				pos += pos / 2;
			}
		}
		return count;
	}
}; 

template< size_t KeyLength, class Value >
class fusion_tree_static<KeyLength, Value>::Node_Leaf : public Node {
public:
	virtual ~Node_Leaf() override = default;

	virtual bool isLeaf() const override {
		return true;
	}

private:
	Value m_data;
};



template<size_t KeyLength, class Value>
constexpr size_t fusion_tree_static<KeyLength, Value>::Node::capacity = (size_t)std::pow(KeyLength, 1.0 / 5.0);

#endif // FUSION_TREE_H

