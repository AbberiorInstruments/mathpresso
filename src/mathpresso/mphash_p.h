// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MATHPRESSO_MPHASH_P_H
#define _MATHPRESSO_MPHASH_P_H

// [Dependencies]
#include  <mathpresso/mathpresso_p.h>

namespace mathpresso
{

	// ============================================================================
	// [mathpresso::HashUtils]
	// ============================================================================

		// For hashing and finding primes
	namespace HashUtils
	{
		// \internal
		static uint32_t hashPointer(const void* kPtr)
		{
			uintptr_t p = (uintptr_t)kPtr;
			return static_cast<uint32_t>(
				((p >> 3) ^ (p >> 7) ^ (p >> 12) ^ (p >> 20) ^ (p >> 27)) & 0xFFFFFFFFU);
		}

		// \internal
		static uint32_t hashChar(uint32_t hash, uint32_t c)
		{
			return hash * 65599 + c;
		}

		// \internal
		//
		// Get a hash of the given string `kStr` of `kLen` length. This function doesn't
		// require `kStr` to be NULL terminated.
		uint32_t hashString(const char* kStr, size_t kLen);

		// \internal
		//
		// Get a prime number that is close to `x`, but always greater than or equal to `x`.
		uint32_t closestPrime(uint32_t x);
	};

	// ============================================================================
	// [mathpresso::HashNode]
	// ============================================================================

	// store a hash and a reference to the next hasNode
	struct HashNode
	{
		HashNode(uint32_t hVal = 0) : _next(nullptr), _hVal(hVal)
		{
		}

		//! Next node in the chain, nullptr if last node.
		HashNode* _next;
		//! Hash code.
		uint32_t _hVal;
	};

	// ============================================================================
	// [mathpresso::HashBase]
	// ============================================================================

	struct HashBase
	{
		MATHPRESSO_NO_COPY(HashBase);

		enum
		{
			kExtraFirst = 0,
			kExtraCount = 1
		};

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		HashBase(ZoneHeap* heap) :
			_bucketsCount(1),
			_heap(heap),
			_length(0),
			_bucketsGrow(1)
		{
			_data = _embedded;
			for (uint32_t i = 0; i <= kExtraCount; i++)
				_embedded[i] = nullptr;
		}

		~HashBase()
		{
			if (_data != _embedded)
				_heap->release(_data, static_cast<size_t>(_bucketsCount + kExtraCount) * sizeof(void*));
		}

		// --------------------------------------------------------------------------
		// [Ops]
		// --------------------------------------------------------------------------
		
		// std::unordered_mpa::rehash
		void _rehash(uint32_t newCount);

		// effectively merge this HashBase to other
		// std::unordered_map::merge (cpp17?)
		void _mergeToInvisibleSlot(HashBase& other);

		// std::unordered_map::emplace
		HashNode* _put(HashNode* node);

		// std::unordered_map::rease
		HashNode* _del(HashNode* node);

		// i would call these internal data -> priave or protected...
		uint32_t _bucketsCount;
		HashNode** _data;

	protected:
		ZoneHeap* _heap;

		uint32_t _length;
		uint32_t _bucketsGrow;
	
		HashNode* _embedded[1 + kExtraCount];
	};

	// ============================================================================
	// [mathpresso::Hash<Key, Node>]
	// ============================================================================

	//! \internal
	//!
	//! Low level hash table container used by MathPresso, with some "special"
	//! features.
	//!
	//! Notes:
	//!
	//! 1. This hash table allows duplicates to be inserted (the API is so low
	//!    level that it's up to you if you allow it or not, as you should first
	//!    `get()` the node and then modify it or insert a new node by using `put()`,
	//!    depending on the intention).
	//!
	//! 2. This hash table also contains "invisible" nodes that are not used by
	//!    the basic hash functions, but can be used by functions having "invisible"
	//!    in their name. These are used by the parser to merge symbols from the
	//!    current scope that is being closed into the root local scope.
	//!
	//! Hash is currently used by AST to keep references of global and local
	//! symbols and by AST to IR translator to associate IR specific data with AST.
	template<typename Key, typename Node>
	struct Hash : public HashBase
	{
		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		Hash(ZoneHeap* heap)
			: HashBase(heap)
		{
		}

		// --------------------------------------------------------------------------
		// [Ops]
		// --------------------------------------------------------------------------

		// std::unordere_map::reset? eventually problems with the release handler?
		// ReleaseHandler = AstScopeReleaseHandler (in mpast.cpp, Z. 175)
		// deletes every AstSymbol, via deleteSymbol(), then resets to defaults
		// called by destructor of AstScope.
		template<typename ReleaseHandler>
		void reset(ReleaseHandler& handler)
		{
			HashNode** data = _data;
			uint32_t count = _bucketsCount + kExtraCount;

			for (uint32_t i = 0; i < count; i++)
			{
				HashNode* node = data[i];

				while (node != nullptr)
				{
					HashNode* next = node->_next;
					handler.release(static_cast<Node*>(node));
					node = next;
				}
			}

			if (data != _embedded)
				_heap->release(data, static_cast<size_t>(count + kExtraCount) * sizeof(void*));

			_bucketsCount = 1;
			_bucketsGrow = 1;

			_length = 0;
			_data = _embedded;

			for (uint32_t i = 0; i <= kExtraCount; i++)
				_embedded[i] = nullptr;
		}

		
		Node* get(const Key& key, uint32_t hVal) const
		{
			uint32_t hMod = hVal % _bucketsCount;
			Node* node = static_cast<Node*>(_data[hMod]);
			
			while (node != nullptr)
			{
				if (node->eq(key))
					return node;
				node = static_cast<Node*>(node->_next);
			}

			return nullptr;
		}
		
		// wrappers, not necessary?
		void mergeToInvisibleSlot(Hash<Key, Node>& other) { _mergeToInvisibleSlot(other); }
		Node* put(Node* node) { return static_cast<Node*>(_put(node)); }
		Node* del(Node* node) { return static_cast<Node*>(_del(node)); }
	};

	// ============================================================================
	// [mathpresso::HashIterator<Key, Node>]
	// ============================================================================
	 // just an Iterator over Hash -> use std::iterator instead
	template<typename Key, typename Node>
	struct HashIterator
	{
		HashIterator(const Hash<Key, Node>& hash) { init(hash); }

		bool init(const Hash<Key, Node>& hash)
		{
			Node** buckets = reinterpret_cast<Node**>(hash._data);
			Node* node = buckets[0];

			uint32_t index = 0;
			uint32_t count = hash._bucketsCount;

			while (node == nullptr && ++index < count)
				node = buckets[index];

			_node = node;
			_buckets = buckets;

			_index = index;
			_count = count;

			return node != nullptr;
		}

		bool next()
		{
			// Can't be called after it reached the end.
			MATHPRESSO_ASSERT(has());

			Node* node = static_cast<Node*>(_node->_next);
			if (node == nullptr)
			{
				uint32_t index = _index;
				uint32_t count = _count;

				while (++index < count)
				{
					node = _buckets[index];
					if (node != nullptr) break;
				}
				_index = index;
			}

			_node = node;
			return node != nullptr;
		}

		bool has() const { return _node != nullptr; }
		Node* get() const { return _node; }

		Node* _node;
		Node** _buckets;

		uint32_t _index;
		uint32_t _count;
	};

} // mathpresso namespace

// [Guard]
#endif // _MATHPRESSO_MPHASH_P_H
