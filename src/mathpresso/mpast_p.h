// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MATHPRESSO_MPAST_P_H
#define _MATHPRESSO_MPAST_P_H

// [Dependencies]
#include <mathpresso/mphash_p.h>
#include <mathpresso/mpoperation.h>

#include <iostream>
#include <complex>
#include <memory>

namespace mathpresso
{

	// ============================================================================
	// [Forward Declarations]
	// ============================================================================

	struct AstBuilder;
	struct AstScope;
	struct AstSymbol;

	struct AstNode;
	struct AstProgram;
	struct AstUnary;


	// ============================================================================
	// [mathpresso::AstScopeType]
	// ============================================================================

	enum AstScopeType
	{
		//! Global scope.
		kAstScopeGlobal = 0,

		//! Shadow scope acts like a global scope, however, it's mutable and can be
		//! modified by the optimizer. Shadow scope is never used to store locals.
		kAstScopeShadow = 1,

		//! Local scope.
		kAstScopeLocal = 2, // unused

		//! Nested scope.
		//!
		//! Always statically allocated and merged with the local scope before the
		//! scope is destroyed.
		kAstScopeNested = 3
	};

	// ============================================================================
	// [mathpresso::AstSymbolType]
	// ============================================================================

	//! \internal
	//!
	//! Symbol type.
	enum AstSymbolType
	{
		//! Not used.
		kAstSymbolNone = 0,
		//! Symbol is an intrinsic (converted to operator internally).
		kAstSymbolIntrinsic,
		//! Symbol is a variable.
		kAstSymbolVariable,
		//! Symbol is a function.
		kAstSymbolFunction
	};

	// ============================================================================
	// [mathpresso::AstSymbolFlags]
	// ============================================================================

	enum AstSymbolFlags
	{
		//! The symbol was declared in global scope.
		kAstSymbolIsGlobal = 0x0001,

		//! The symbol was declared and can be used.
		//!
		//! If this flag is not set it means that the parser is parsing its assignment
		//! (code like "int x = ...") and the symbol can't be used at this time. It's
		//! for parser to make sure that the symbol is declared before it's used.
		kAstSymbolIsDeclared = 0x0002,

		//! Used during optimizing phase and to create global constants.
		kAstSymbolIsAssigned = 0x0004,

		//! The symbol (variable) is read-only.
		kAstSymbolIsReadOnly = 0x0008,


		//! The variable has been altered (written), at least once.
		//!
		//! Currently only useful for global variables so the JIT compiler can
		//! perform write operation at the end of the generated function.
		kAstSymbolIsAltered = 0x0010,


		//! The variable is a complex value.
		kAstSymbolIsComplex = 0x0020,


		//! The function returns a complex value
		kAstSymbolRealFunctionReturnsComplex = 0x00040, // unused
		//! The function takes complex arguments
		kAstSymbolComplexFunctionReturnsReal = 0x00080, // unused
		//! See function flags
		kAstSymbolHasState = 0x00100 // unused
	};

	// ============================================================================
	// [mathpresso::AstNodeType]
	// ============================================================================

	//! \internal
	//!
	//! `AstNode` type.
	enum AstNodeType
	{
		//! Not used.
		kAstNodeNone = 0,

		// --------------------------------------------------------------------------
		// [Block]
		// --------------------------------------------------------------------------

		//! Node is `AstProgram`.
		kAstNodeProgram,
		//! Node is `AstBlock`.
		kAstNodeBlock,

		// --------------------------------------------------------------------------
		// [Variable, Immediate]
		// --------------------------------------------------------------------------

		//! Node is `AstVarDecl`.
		kAstNodeVarDecl,
		//! Node is `AstVar` of type double.
		kAstNodeVar,
		//! Node is `AstImm`.
		kAstNodeImm,

		// --------------------------------------------------------------------------
		// [Op]
		// --------------------------------------------------------------------------

		//! Node is `AstUnaryOp`.
		kAstNodeUnaryOp,
		//! Node is `AstBinaryOp`.
		kAstNodeBinaryOp,
		//! Node is `AstTernaryOp`.
		kAstNodeTernaryOp,
		//! Node is `AstCall`.
		kAstNodeCall
	};

	// ============================================================================
	// [mathpresso::AstNodeFlags]
	// ============================================================================

	//! \internal
	//!
	//! `AstNode` flags.
	enum AstNodeFlags
	{
		kAstNone = 0x00,
		kAstNodeHasSideEffect = 0x01,
		kAstTakesComplex = 0x02, // set, if a complex 'parameter' is expected.
		kAstReturnsComplex = 0x04 // set, if the 'return' is complex.
	};

	// ============================================================================
	// [mathpresso::AstBuilder]
	// ============================================================================

	//! \internal
	struct AstBuilder
	{
		MATHPRESSO_NO_COPY(AstBuilder);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstBuilder(ZoneHeap* heap);
		~AstBuilder();

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		ZoneHeap* getHeap() const { return _heap; }

		AstScope* getRootScope() const { return _rootScope; }
		AstProgram* getProgramNode() const { return _programNode; }

		// --------------------------------------------------------------------------
		// [Factory]
		// --------------------------------------------------------------------------

		AstScope* newScope(AstScope* parent, uint32_t scopeType);
		void deleteScope(AstScope* scope);

		AstSymbol* newSymbol(const std::string& key, uint32_t hVal, uint32_t symbolType, uint32_t scopeType);
		AstSymbol* shadowSymbol(const AstSymbol* other);
		void deleteSymbol(AstSymbol* symbol);

#define MATHPRESSO_ALLOC_AST_OBJECT(_Size_) \
  void* obj = _heap->alloc(_Size_); \
  if (MATHPRESSO_UNLIKELY(obj == nullptr)) return nullptr

		template<typename T>
		T* newNode()
		{
			MATHPRESSO_ALLOC_AST_OBJECT(sizeof(T));
			return new(obj) T(this);
		}

		template<typename T, typename P0>
		T* newNode(P0 p0)
		{
			MATHPRESSO_ALLOC_AST_OBJECT(sizeof(T));
			return new(obj) T(this, p0);
		}

		template<typename T, typename P0, typename P1>
		T* newNode(P0 p0, P1 p1)
		{
			MATHPRESSO_ALLOC_AST_OBJECT(sizeof(T));
			return new(obj) T(this, p0, p1);
		}

#undef MATHPRESSO_ALLOC_AST_OBJECT

		void deleteNode(AstNode* node);

		uint32_t newSlotId() { return _numSlots++; }

		// --------------------------------------------------------------------------
		// [Init]
		// --------------------------------------------------------------------------

		Error initProgramScope();

		// --------------------------------------------------------------------------
		// [Dump]
		// --------------------------------------------------------------------------

		Error dump(StringBuilder& sb, const Operations * ops);

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		//! Heap.
		ZoneHeap* _heap;
		//! String builder to build possible output messages.
		StringBuilder _sb;

		//! Root scope.
		AstScope* _rootScope;
		//! Root node.
		AstProgram* _programNode;

		//! Number of variable slots used.
		uint32_t _numSlots;
	};

	// ============================================================================
	// [mathpresso::AstSymbol]
	// ============================================================================

	struct AstSymbol : public HashNode
	{
		MATHPRESSO_NO_COPY(AstSymbol);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstSymbol(const char* name, uint32_t length, uint32_t hVal, uint32_t symbolType, uint32_t scopeType)
			: HashNode(hVal),
			_length(length),
			_name(name),
			_node(nullptr),
			_symbolType(static_cast<uint8_t>(symbolType)),
			_symbolFlags(scopeType == AstScopeType::kAstScopeGlobal ? (int)AstSymbolFlags::kAstSymbolIsGlobal : 0),
			_valueComp(),
			_usedCount(0),
			_writeCount(0)
		{
		}



		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		bool eq(const std::string& s) const
		{
			return eq(s.c_str(), s.length());
		}

		//! Get whether the symbol name is equal to string `s` of `len`.
		bool eq(const char* s, size_t len) const
		{
			return static_cast<size_t>(_length) == len && ::memcmp(_name, s, len) == 0;
		}

		//! Get symbol name length.
		uint32_t getLength() const { return _length; }
		//! Get symbol name.
		const char* getName() const { return _name; }

		//! Check if the symbol has associated node with it.
		bool hasNode() const { return _node != nullptr; }
		//! Get node associated with the symbol (can be `NULL` for built-ins).
		AstNode* getNode() const { return _node; }
		//! Associate node with the symbol (basically the node that declares it).
		void setNode(AstNode* node) { _node = node; }

		//! Get hash value of the symbol name.
		uint32_t getHVal() const { return _hVal; }

		//! Get symbol type, see \ref AstSymbolType.
		uint32_t getSymbolType() const { return _symbolType; }

		//! Get symbol flags, see \ref AstSymbolFlags.
		uint32_t getSymbolFlags() const { return _symbolFlags; }

		bool hasSymbolFlag(uint32_t flag) const { return (_symbolFlags & flag) != 0; }
		void setSymbolFlag(uint32_t flag) { _symbolFlags |= static_cast<uint16_t>(flag); }
		void clearSymbolFlag(uint32_t flag) { _symbolFlags &= ~static_cast<uint16_t>(flag); }

		//! Check if the symbol is global (i.e. it was declared in a global scope).
		bool isGlobal() const { return hasSymbolFlag(AstSymbolFlags::kAstSymbolIsGlobal); }
		//! Check if the symbol was declared.
		bool isDeclared() const { return hasSymbolFlag(AstSymbolFlags::kAstSymbolIsDeclared); }
		//! Set the symbol to be declared (\ref kAstSymbolIsDeclared flag).
		void setDeclared() { setSymbolFlag(AstSymbolFlags::kAstSymbolIsDeclared); }

		uint32_t getVarSlotId() const { return _varSlotId; }
		void setVarSlotId(uint32_t slotId) { _varSlotId = slotId; }

		int32_t getVarOffset() const { return _varOffset; }
		void setVarOffset(int32_t offset) { _varOffset = offset; }

		//! Get whether the variable has assigned a constant value at the moment.
		//!
		//! If true, the `_value` is a valid constant that can be used to replace
		//! the variable node by a constant value. The value can change during AST
		//! traversal in case that the variable is mutable.
		bool isAssigned() const { return hasSymbolFlag(AstSymbolFlags::kAstSymbolIsAssigned); }
		//! Set symbol to be assigned (sets the \ref kAstSymbolIsAssigned flag).
		void setAssigned() { setSymbolFlag(AstSymbolFlags::kAstSymbolIsAssigned); }
		//! Set symbol to not be assigned (clears the \ref kAstSymbolIsAssigned flag).
		void clearAssigned() { clearSymbolFlag(AstSymbolFlags::kAstSymbolIsAssigned); }

		//! Get whether the global symbol has been altered.
		bool isAltered() const { return hasSymbolFlag(AstSymbolFlags::kAstSymbolIsAltered); }
		//! Make a global symbol altered.
		void setAltered() { setSymbolFlag(AstSymbolFlags::kAstSymbolIsAltered); }

		//! Get the constant value, see `isAssigned()`.
		double getValue() const { return _valueComp.real(); }
		std::complex<double> getValueComp() const { return _valueComp; }

		void setValue(double value) { _valueComp.real(value); }
		void setValue(std::complex<double> value) { _valueComp = value; }

		uint32_t getUsedCount() const { return _usedCount; }
		uint32_t getReadCount() const { return _usedCount - _writeCount; }
		uint32_t getWriteCount() const { return _writeCount; }

		void incUsedCount(uint32_t n = 1) { _usedCount += n; }
		void incWriteCount(uint32_t n = 1) { _writeCount += n; _usedCount += n; }

		void decUsedCount(uint32_t n = 1) { _usedCount -= n; }
		void decWriteCount(uint32_t n = 1) { _writeCount -= n; _usedCount -= n; }

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

	private:
		//! Symbol name length.
		uint32_t _length;
		//! Symbol name (key).
		const char* _name;

		//! Node where the symbol is defined.
		AstNode* _node;

		//! Type of the symbol, see \ref AstSymbolType.
		uint8_t _symbolType;
		//! Flags, see \ref AstSymbolFlags.
		uint16_t _symbolFlags;


		//! Number of times the variable is used (both read and write count).
		uint32_t _usedCount;
		//! Number of times the variable is written.
		uint32_t _writeCount;

		//! Variable slot id.
		uint32_t _varSlotId;
		//! Variable offset in data structure (in case the symbol is a global variable).
		int32_t _varOffset;
		//! The current value of the symbol (in case the symbol is an immediate).
		//! if the symbol is real, _valueComp.imag() is set to 0.
		std::complex<double> _valueComp;
	};

	typedef Hash<std::string, AstSymbol> AstSymbolHash;
	typedef HashIterator<std::string, AstSymbol> AstSymbolHashIterator;

	// ============================================================================
	// [mathpresso::AstScope]
	// ============================================================================

	struct AstScope
	{
		MATHPRESSO_NO_COPY(AstScope);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		MATHPRESSO_NOAPI AstScope(AstBuilder* ast, AstScope* parent, uint32_t scopeType);
		MATHPRESSO_NOAPI ~AstScope();

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		//! Get the scope context.
		AstBuilder* getAst() const { return _ast; }
		//! Get the parent scope (or NULL).
		AstScope* getParent() const { return _parent; }
		//! Get symbols hash-table.
		const AstSymbolHash& getSymbols() const { return _symbols; }
		//! Get scope type, see \ref AstScopeType.
		uint32_t getScopeType() const { return _scopeType; }

		//! Get whether the scope type is `AstScopeType::kAstScopeGlobal`.
		bool isGlobal() const { return _scopeType == AstScopeType::kAstScopeGlobal; }

		//! Make this scope a shadow of `ctxScope`.
		void shadowContextScope(AstScope* ctxScope)
		{
			_parent = ctxScope;
			_scopeType = AstScopeType::kAstScopeShadow;
		}

		// --------------------------------------------------------------------------
		// [Ops]
		// --------------------------------------------------------------------------

		//! Get the symbol defined only in this scope.
		AstSymbol* getSymbol(const std::string& name, uint32_t hVal)
		{
			return _symbols.get(name, hVal);
		}

		//! Put a given symbol to this scope.
		//!
		//! NOTE: The function doesn't care about duplicates. The correct flow is
		//! to call `resolveSymbol()` or `getSymbol()` and then `putSymbol()` based
		//! on the result. You should never call `putSymbol()` without checking if
		//! the symbol is already there.
		void putSymbol(AstSymbol* symbol)
		{
			_symbols.put(symbol);
		}

		//! Resolve the symbol by traversing all parent scopes if not found in this
		//! one. An optional `scopeOut` argument can be used to get scope where the
		//! `name` has been found.
		MATHPRESSO_NOAPI AstSymbol* resolveSymbol(const std::string& name, uint32_t hVal, AstScope** scopeOut = nullptr);

		AstSymbol* resolveSymbol(const std::string& name)
		{
			return resolveSymbol(name, HashUtils::hashString(name.c_str(), name.length()));
		}

		MATHPRESSO_NOAPI AstSymbol* removeSymbol(AstSymbol* symbol)
		{
			return _symbols.del(symbol);
		}

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		//! Context.
		AstBuilder* _ast;
		//! Parent scope.
		AstScope* _parent;

		//! Symbols defined in this scope.
		AstSymbolHash _symbols;

		//! Scope type, see \ref AstScopeType.
		uint32_t _scopeType;
	};

	// ============================================================================
	// [mathpresso::AstNode]
	// ============================================================================

#define MATHPRESSO_AST_CHILD(_Index_, _Type_, _Name_, _Memb_) \
  bool has##_Name_() const { return _Memb_ != nullptr; } \
  _Type_* get##_Name_() const { return _Memb_; } \
  \
  _Type_* set##_Name_(_Type_* node) { \
    return static_cast<_Type_*>(replaceAt(_Index_, node)); \
  } \
  \
  _Type_* unlink##_Name_() { \
    _Type_* node = _Memb_; \
    \
    MATHPRESSO_ASSERT(node != nullptr); \
    MATHPRESSO_ASSERT(node->_parent == this); \
    \
    node->_parent = nullptr; \
    _Memb_ = nullptr; \
    \
    return node; \
  } \
  private:\
  _Type_* _Memb_;\
  public:


	struct AstNode
	{
		MATHPRESSO_NO_COPY(AstNode);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstNode(AstBuilder* ast, uint32_t nodeType, AstNode** children = nullptr, uint32_t length = 0)
			: _ast(ast),
			_parent(nullptr),
			_children(children),
			_mpOp(nullptr),
			_opName(),
			_nodeType(static_cast<uint8_t>(nodeType)),
			_nodeFlags(AstNodeFlags::kAstNone),
			_position(~static_cast<uint32_t>(0)),
			_length(length)
		{
		}

		virtual ~AstNode()
		{
		}

		void destroy(AstBuilder* ast)
		{}

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		//! Get the `AstBuilder` instance that created this node.
		AstBuilder* getAst() const { return _ast; }

		//! Check if the node has a parent.
		bool hasParent() const { return _parent != nullptr; }
		//! Get the parent node.
		AstNode* getParent() const { return _parent; }

		//! Get whether the node has children.
		//!
		//! NOTE: Nodes that always have children (even if they are implicitly set
		//! to NULL) always return `true`. This function if useful mostly if the
		//! node is of `AstBlock` type.
		bool hasChildren() const { return _length != 0; }
		//! Get children array.
		AstNode** getChildren() const { return reinterpret_cast<AstNode**>(_children); }
		//! Get length of the children array.
		size_t getLength() const { return _length; }

		//! Get node type.
		uint32_t getNodeType() const { return _nodeType; }
		//! Get whether the node is `AstVar`.
		bool isVar() const { return _nodeType == AstNodeType::kAstNodeVar; }
		//! Get whether the node is `AstImm`.
		bool isImm() const { return _nodeType == AstNodeType::kAstNodeImm; }

		//! Get whether the node has flag `flag`.
		bool hasNodeFlag(uint32_t flag) const { return (static_cast<uint32_t>(_nodeFlags) & flag) != 0; }
		//! Get node flags.
		uint32_t getNodeFlags() const { return _nodeFlags; }
		//! Set node flags.
		void setNodeFlags(uint32_t flags) { _nodeFlags = static_cast<uint8_t>(flags); }
		//! Add node flags.
		void addNodeFlags(uint32_t flags) { _nodeFlags |= static_cast<uint8_t>(flags); }
		//! remove a flag.
		void removeNodeFlags(uint32_t flags) { _nodeFlags &= ~static_cast<uint8_t>(flags); }

		bool takesComplex()   const { return hasNodeFlag(AstNodeFlags::kAstTakesComplex); }
		bool returnsComplex() const { return hasNodeFlag(AstNodeFlags::kAstReturnsComplex); }

		//! Get whether the node has associated position in source code.
		bool hasPosition() const { return _position != ~static_cast<uint32_t>(0); }
		//! Get source code position of the node.
		uint32_t getPosition() const { return _position; }
		//! Set source code position of the node.
		void setPosition(uint32_t position) { _position = position; }
		//! Reset source code position of the node.
		void resetPosition() { _position = ~static_cast<uint32_t>(0); }

		// --------------------------------------------------------------------------
		// [Children]
		// --------------------------------------------------------------------------

		AstNode* getAt(size_t index) const
		{
			MATHPRESSO_ASSERT(index < _length);
			return _children[index];
		}

		//! Replace `refNode` by `node`.
		AstNode* replaceNode(AstNode* refNode, AstNode* node);
		//! Replace node at index `index` by `node`.
		AstNode* replaceAt(uint32_t index, AstNode* node);

		//! Inject `node` between this node and `refNode`.
		AstNode* injectNode(AstNode* refNode, AstUnary* node);
		//! Inject `node` between this node and node at index `index`.
		AstNode* injectAt(uint32_t index, AstUnary* node);

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		//! AST builder.
		AstBuilder* _ast;
		//! Parent node.
		AstNode* _parent;
		//! Child nodes.
		AstNode** _children;

		std::shared_ptr<MpOperation> _mpOp;

		std::string _opName;

	private:
		//! Node type, see `AstNodeType`.
		uint8_t _nodeType;
		//! Node flags, see `AstNodeFlags`.
		uint8_t _nodeFlags;
		//! Position (in characters) to the beginning of the program (default -1).
		uint32_t _position;

	protected:
		//! Count of child-nodes.
		size_t _length;
	};

	// ============================================================================
	// [mathpresso::AstBlock]
	// ============================================================================

	struct AstBlock : public AstNode
	{
		MATHPRESSO_NO_COPY(AstBlock);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstBlock(AstBuilder* ast, uint32_t nodeType = AstNodeType::kAstNodeBlock)
			: AstNode(ast, nodeType),
			_capacity(0)
		{
		}

		// --------------------------------------------------------------------------
		// [Ops]
		// --------------------------------------------------------------------------

		//! Reserve the capacity of the AstBlock so one more node can be added into it.
		//!
		//! NOTE: This has to be called before you use `appendNode()` or `insertAt()`.
		//! The reason is that it's easier to deal with possible allocation failure
		//! here (before the node to be added is created) than after the node is
		//! created, but failed to add into the block.
		Error willAdd();

		//! Append the given `node` to the block.
		//!
		//! NOTE: You have to call `willAdd()` before you use `appendNode()` for every
		//! node you want to add to the block.
		void appendNode(AstNode* node)
		{
			MATHPRESSO_ASSERT(node != nullptr);
			MATHPRESSO_ASSERT(node->getParent() == nullptr);

			// We expect `willAdd()` to be called before `appendNode()`.
			MATHPRESSO_ASSERT(_length < _capacity);

			node->_parent = this;

			_children[_length] = node;
			_length++;
		}

		//! Insert the given `node` to the block at index `i`.
		//!
		//! NOTE: You have to call `willAdd()` before you use `insertAt()` for every
		//! node you want to add to the block.
		void insertAt(size_t i, AstNode* node)
		{
			MATHPRESSO_ASSERT(node != nullptr);
			MATHPRESSO_ASSERT(node->getParent() == nullptr);

			// We expect `willAdd()` to be called before `insertAt()`.
			MATHPRESSO_ASSERT(_length < _capacity);

			AstNode** p = getChildren();
			node->_parent = this;

			size_t j = _length;
			while (i < j)
			{
				p[j] = p[j - 1];
				j--;
			}

			p[j] = node;
			_length++;
		}

		//! Remove the given `node`.
		AstNode* removeNode(AstNode* node);
		//! Remove the node at index `index`.
		AstNode* removeAt(size_t index);

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		uint32_t _capacity;
	};

	// ============================================================================
	// [mathpresso::AstUnary]
	// ============================================================================

	struct AstUnary : public AstNode
	{
		MATHPRESSO_NO_COPY(AstUnary);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstUnary(AstBuilder* ast, uint32_t nodeType)
			: AstNode(ast, nodeType, &_child, 1),
			_child(nullptr)
		{
		}

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		AstNode** getChildren() const { return (AstNode**)&_child; }

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		MATHPRESSO_AST_CHILD(0, AstNode, Child, _child);
	};

	// ============================================================================
	// [mathpresso::AstBinary]
	// ============================================================================

	struct AstBinary : public AstNode
	{
		MATHPRESSO_NO_COPY(AstBinary);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstBinary(AstBuilder* ast, uint32_t nodeType)
			: AstNode(ast, nodeType, &_left, 2),
			_left(nullptr),
			_right(nullptr)
		{
		}

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		AstNode** getChildren() const { return (AstNode**)&_left; }

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		MATHPRESSO_AST_CHILD(0, AstNode, Left, _left);
		MATHPRESSO_AST_CHILD(1, AstNode, Right, _right);
	};


	// ============================================================================
	// [mathpresso::AstTernary]
	// ============================================================================

	struct AstTernary : public AstNode
	{
		MATHPRESSO_NO_COPY(AstTernary);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstTernary(AstBuilder* ast, uint32_t nodeType)
			: AstNode(ast, nodeType, &_condition, 3),
			_condition(nullptr),
			_left(nullptr),
			_right(nullptr)
		{
		}

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		AstNode** getChildren() const { return (AstNode**)&_condition; }

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		MATHPRESSO_AST_CHILD(0, AstNode, Condition, _condition);
		MATHPRESSO_AST_CHILD(1, AstNode, Left, _left);
		MATHPRESSO_AST_CHILD(2, AstNode, Right, _right);
	};

	// ============================================================================
	// [mathpresso::AstProgram]
	// ============================================================================

	struct AstProgram : public AstBlock
	{
		MATHPRESSO_NO_COPY(AstProgram);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstProgram(AstBuilder* ast)
			: AstBlock(ast, AstNodeType::kAstNodeProgram)
		{
		}


	};

	// ============================================================================
	// [mathpresso::AstVarDecl]
	// ============================================================================

	struct AstVarDecl : public AstUnary
	{
		MATHPRESSO_NO_COPY(AstVarDecl);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstVarDecl(AstBuilder* ast)
			: AstUnary(ast, AstNodeType::kAstNodeVarDecl),
			_symbol(nullptr)
		{
		}

		void destroy(AstBuilder* ast)
		{
			AstSymbol* sym = getSymbol();
			if (sym != nullptr)
				sym->decUsedCount();
		}

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		AstSymbol* getSymbol() const { return _symbol; }
		void setSymbol(AstSymbol* symbol) { _symbol = symbol; }

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------
	private:
		AstSymbol* _symbol;

	};

	// ============================================================================
	// [mathpresso::AstVar]
	// ============================================================================

	struct AstVar : public AstNode
	{
		MATHPRESSO_NO_COPY(AstVar);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstVar(AstBuilder* ast)
			: AstNode(ast, AstNodeType::kAstNodeVar),
			_symbol(nullptr)
		{
		}

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		AstSymbol* getSymbol() const { return _symbol; }
		void setSymbol(AstSymbol* symbol) { _symbol = symbol; }

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------
	private:
		AstSymbol* _symbol;
	};


	// ============================================================================
	// [mathpresso::AstImm]
	// ============================================================================

	struct AstImm : public AstNode
	{
		MATHPRESSO_NO_COPY(AstImm);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstImm(AstBuilder* ast, double value = 0.0)
			: AstNode(ast, AstNodeType::kAstNodeImm),
			_value({ value, 0 })
		{
		}

		AstImm(AstBuilder* ast, std::complex<double> value)
			: AstNode(ast, AstNodeType::kAstNodeImm),
			_value(value)
		{
			addNodeFlags(AstNodeFlags::kAstReturnsComplex);
		}

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		template<typename T>
		T getValue() const;
		

		void setValue(double value)
		{
			_value = { value, 0 };
			removeNodeFlags(kAstReturnsComplex);
		}
		void setValue(std::complex<double> value)
		{
			_value = value;
			addNodeFlags(kAstReturnsComplex);
		}

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------
	private:
		std::complex<double> _value;
	};



	// ============================================================================
	// [mathpresso::AstUnaryOp]
	// ============================================================================

	struct AstUnaryOp : public AstUnary
	{
		MATHPRESSO_NO_COPY(AstUnaryOp);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstUnaryOp(AstBuilder* ast)
			: AstUnary(ast, AstNodeType::kAstNodeUnaryOp)
		{
		}

	};

	// ============================================================================
	// [mathpresso::AstBinaryOp]
	// ============================================================================

	struct AstBinaryOp : public AstBinary
	{
		MATHPRESSO_NO_COPY(AstBinaryOp);

		AstBinaryOp(AstBuilder* ast) : AstBinary(ast, AstNodeType::kAstNodeBinaryOp)
		{
		}

		void destroy(AstBuilder* ast)
		{
			if (_mpOp && (_mpOp->flags() & MpOperation::IsAssignment) && hasLeft())
			{
				AstVar* var = static_cast<AstVar*>(getLeft());
				AstSymbol* sym = var->getSymbol();

				if (sym != nullptr)
					sym->decWriteCount();
			}
		}
	};

	// ============================================================================
	// [mathpresso::AstTernaryOp]
	// ============================================================================

	struct AstTernaryOp : public AstTernary
	{
		MATHPRESSO_NO_COPY(AstTernaryOp);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstTernaryOp(AstBuilder* ast) :
			AstTernary(ast, AstNodeType::kAstNodeTernaryOp)
		{
		}

	};

	// ============================================================================
	// [mathpresso::AstCall]
	// ============================================================================

	struct AstCall : public AstBlock
	{
		MATHPRESSO_NO_COPY(AstCall);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstCall(AstBuilder* ast)
			: AstBlock(ast, AstNodeType::kAstNodeCall),
			_symbol(nullptr)
		{
		}

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		AstSymbol* getSymbol() const { return _symbol; }
		void setSymbol(AstSymbol* symbol) { _symbol = symbol; }

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		AstSymbol* _symbol;
	};

	// ============================================================================
	// [mathpresso::AstVisitor]
	// ============================================================================

	struct AstVisitor
	{
		MATHPRESSO_NO_COPY(AstVisitor);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstVisitor(AstBuilder* ast);
		virtual ~AstVisitor();

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		AstBuilder* getAst() const { return _ast; }

		// --------------------------------------------------------------------------
		// [OnNode]
		// --------------------------------------------------------------------------

		virtual Error onNode(AstNode* node);

		virtual Error onProgram(AstProgram* node);
		virtual Error onBlock(AstBlock* node) = 0;
		virtual Error onVarDecl(AstVarDecl* node) = 0;
		virtual Error onVar(AstVar* node) = 0;
		virtual Error onImm(AstImm* node) = 0;
		virtual Error onUnaryOp(AstUnaryOp* node) = 0;
		virtual Error onBinaryOp(AstBinaryOp* node) = 0;
		virtual Error onTernaryOp(AstTernaryOp* node) = 0;
		virtual Error onCall(AstCall* node) = 0;

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		AstBuilder* _ast;
	};

	// ============================================================================
	// [mathpresso::AstDump]
	// ============================================================================

	struct AstDump : public AstVisitor
	{
		MATHPRESSO_NO_COPY(AstDump);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstDump(AstBuilder* ast, StringBuilder& sb, const Operations * ctx);
		virtual ~AstDump();

		// --------------------------------------------------------------------------
		// [OnNode]
		// --------------------------------------------------------------------------

		virtual Error onBlock(AstBlock* node);
		virtual Error onVarDecl(AstVarDecl* node);
		virtual Error onVar(AstVar* node);
		virtual Error onImm(AstImm* node);
		virtual Error onUnaryOp(AstUnaryOp* node);
		virtual Error onBinaryOp(AstBinaryOp* node);
		virtual Error onTernaryOp(AstTernaryOp * node);
		virtual Error onCall(AstCall* node);

		// --------------------------------------------------------------------------
		// [Helpers]
		// --------------------------------------------------------------------------

		Error info(const char* fmt, ...);
		Error nest(const char* fmt, ...);
		Error denest();

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		StringBuilder& _sb;
		uint32_t _level;
		const Operations * _ops;
	};

} // mathpresso namespace

// [Guard]
#endif // _MATHPRESSO_MPAST_P_H
