#ifndef FUSION_TREE_H
#define FUSION_TREE_H

#include <bitset>
#include <cmath>
#include <array>
#include <vector>
#include <cstdint>

/* KeyLength == log(u) where u is |U| - universe size
 * See http://courses.csail.mit.edu/6.897/spring03/scribe_notes/L4/lecture4.pdf
 * for conceptual notes.
 * This implementation accepts KeyLength value up to 2^64.
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
	class iterator;
	class const_iterator;
	using reference = Value&;
	using const_reference = const Value&;


	fusion_tree_static() :
		m_root(nullptr)
	{}
	~fusion_tree_static() {
		if (m_root != nullptr) {
			delete m_root;
		}
	}

	const_reference at(const Key& key) const {
		 iterator it(m_root->findBucket(key));
		 return *it;
	}
	reference at(const Key& key) {
		iterator it(m_root->findBucket(key));
		return *it;
	}

	const_iterator cbegin() const {
		return const_iterator(m_root ? findMinNode(m_root) : nullptr);
	}
	const_iterator begin() const {
		return cbegin();
	}
	iterator begin() const {
		return iterator(m_root ? findMinNode(m_root) : nullptr);
	}

	const_iterator cend() const {
		return const_iterator(nullptr);
	}
	const_iterator end() const {
		return cend();
	}
	iterator end() const {
		return iterator(nullptr);
	}

	iterator insert(const Key& key, const Value& value) {
		// TO DO:

	}

	const_iterator find(const Key& key) const {
		if (empty() == true) {
			return cend();
		}
		return const_iterator(m_root->find(key));
	}
	iterator find(const Key& key) {
		if (empty() == true) {
			return cend();
		}
		return iterator(m_root->find(key));
	}

	bool empty() const {
		return m_root == nullptr;
	}

private:
	class Node;
	class Node_Regular;
	class Node_Leaf;
	

	Node_Regular* m_root;
};

template< size_t KeyLength, class Value >
class fusion_tree_static<KeyLength, Value>::Node /*Abstract*/ {
public:
	Node() = default;
	virtual Node(const Node& other) = default;
	virtual Node(Node&& other) = default;
	virtual ~Node() = default;

	virtual Node& operator=(const Node& other) = default;
	virtual Node& operator=(Node&& other) = default;

	virtual static bool isLeaf() const = 0;
};

template< size_t KeyLength, class Value >
class fusion_tree_static<KeyLength, Value>::Node_Regular : public Node {
private:
	static constexpr size_t capacity = (size_t)std::pow(KeyLength, 1.0 / 5.0);
	using Sketch = bitset<KeyLength>;
	struct Children {
		using BinHeap = bitset<capacity + (capacity - 1)>;

		/* Represents a binary search tree as a binary heap.
		* If Tree[x] == true then node x is exists.
		*/
		using Tree = BinHeap;
		using Entry = std::pair<Key, Node*>;
		using Values = vector<Entry>;

		Tree tree;
		Values values;

		Children() :
			values()
		{
			values.reserve(capacity);
		}

		~Children() {
			for (Entry& e : values) {
				if (e.second != nullptr) {
					delete e.second;
				}
			}
		}

		const Entry& operator[](const size_t childIndex) const {
			return values[childIndex];
		}

		Entry& operator[](const size_t childIndex) {
			return values[childIndex];
		}

		size_t count() const {
			return values.size();
		}
	};

public:
	virtual ~Node_Regular() override = default;

	virtual static bool isLeaf() const override {
		return false;
	};

	void insert(const Key& key, const Value& value) {
		Node* node = find(key);
		if (node != nullptr) {
			static_cast<Node_Leaf*>(node)->data() = value;
			return node;
		}


	}

	bool contains(const Key& key) const {
		return find(key) != nullptr;
	}

	Node_Leaf* find(const Key& key) const {
		if (m_children.count() == 0) {
			return nullptr;
		}

		Children::Entry& entry = findBucket(key);
		if (entry.second->isLeaf() == true) {
			if (entry.first == key) {
				return static_cast<Node_Leaf*>(entry.second);
			} else {
				return nullptr;
			}
		} else {
			return entry.second->find(key);
		}
	}

private:
	Sketch m_sketch;
	Children m_children;
	
	enum class Endianness : bool {
		bigEndian = 0,
		littleEndian = 1		
	};
	static Endianness WORD_ORDER = Endianness::bigEndian;
	static Endianness determineEndianness() {
		volatile int num = 1; // disable compiler optimization
		WORD_ORDER = static_cast<Endianness>(*(char *)&num == 1);
		return WORD_ORDER;
	}

	/* Returns position of Most Significant Bit set in number.
	*/
	template< class Number >
	inline static size_t findMSB(const Number& number) {
		/*
			The code loads a 64-bit (IEEE-754 floating-point) double with a 32-bit integer
			(with no paddding bits) by storing the integer in the mantissa while the exponent is set to 2^52.
			From this newly minted double, 2^52 (expressed as a double) is subtracted,
			which sets the resulting exponent to the log base 2 of the input value, number.
			All that is left is shifting the exponent bits into position (20 bits right) and
			subtracting the bias, 0x3FF (which is 1023 decimal).

			This technique only takes 5 operations, but many CPUs are slow at manipulating doubles,
			and the endianess of the architecture must be accommodated.
		*/
		determineEndianness();

		std::uint32_t u[2] = {0, 0}; 
		double d; 
		
		u[WORD_ORDER == Endianness::littleEndian] = 0x43300000;
		u[WORD_ORDER != Endianness::littleEndian] = static_cast<unsigned int>(number);
		memcpy(&d, u, sizeof(double));
		d -= 4503599627370496.0;
		memcpy(u, &d, sizeof(double));
		size_t result = (u[WORD_ORDER == Endianness::littleEndian] >> 20) - 0x3FF;
		return result;
	}

	Children::Entry& findBucket(const Key& key) {
		Key sketch;
		computeKeySketch(key, sketch);

		size_t blockLength = m_children.count() * m_children.count() * m_children.count() * m_children.count() + 1u;

		unsigned long long replicatedSketch = 0;
		for (size_t i = 0; i < m_children.count(); ++i) {
			replicatedSketch += 1ull << (i * blockLength);
		}

		unsigned long long mask = 0;
		for (size_t i = 0; i < m_children.count() - 1u; ++i) {
			mask += 1ull << (i * blockLength + blockLength - 1u);
		}

		replicatedSketch *= static_cast<unsigned long long>(sketch);
		replicatedSketch = static_cast<unsigned long long>(m_sketch) - replicatedSketch;
		replicatedSketch &= mask;

		size_t; childIndex = m_children.count() - 1;
		if (replicatedSketch != 0) {
			childIndex = findMSB(replicatedSketch);
		}

		return m_children[childIndex];
	}

	/* Computes a node sketch.
	 * Have to be called after every node modification.
	 * Complexity: O(capacity ^ 4)
	*/
	void computeNodeSketch() {
		m_sketch.reset();

		if (m_children.count() < 2u) {
			return;
		}		
		
		Key sketch;
		for (size_t i = 0; i < m_children.count(); ++i) {
			// Compute child's sketch
			computeKeySketch(m_children[i].first, sketch);
			m_sketch |= (sketch << i);
		}
	}

	void computeKeySketch(const Key& key, Key& sketch) const {
		Key branchingPoints;
		const size_t branchingPointsCount = determineBranchingPoints(key, branchingPoints);

		vector<size_t> branchingPointsPositions;
		for (size_t i = 0; i < branchingPointsCount; ++i) {
			if (branchingPoints[point] == true) {
				branchingPointsPositions.push_back(i);
			}
		}

		sketch = branchingPoints & key;

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

		sketch = (static_cast<unsigned long long>(sketch) * m) >> minPower;
		sketch.set(mParts.size());
	}

	size_t determineBranchingPoints(const Key& nodeKey, Key& branchingPoints) const {
		size_t keyPos = 0;
		size_t treePos = 0;
		constexpr size_t nonLeafHeight = capacity - 1;
		while (treePos < nonLeafHeight) {			
			const size_t sibling = treePos * 2 + (2 - nodeKey[keyPos]);

			// There is a branching point if sibling node exists
			branchingPoints.set(keyPos, m_children.tree[sibling]);
			
			treePos *= 2;
			treePos += 1 + nodeKey[keyPos];
			++keyPos;
		}
		return count;
	}
}; 

template< size_t KeyLength, class Value >
class fusion_tree_static<KeyLength, Value>::Node_Leaf : public Node {
public:
	Node_Leaf(const Value& value) :
		m_value(value)
	{}
	Node_Leaf(Value&& value) :
		m_value(std::move(value))
	{}

	virtual ~Node_Leaf() override = default;

	virtual static bool isLeaf() const override {
		return true;
	}

	const Value& data() const {
		return m_data;
	}
	Value& data() {
		return m_data;
	}

private:
	Value m_data;
};


#endif // FUSION_TREE_H

