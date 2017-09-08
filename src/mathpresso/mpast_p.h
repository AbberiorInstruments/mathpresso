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
	enum class AstSymbolType
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
	enum class AstNodeType
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
	struct AstBuilder : public std::enable_shared_from_this<AstBuilder>
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

		std::shared_ptr<AstProgram> getProgramNode() const { return _programNode; }

		// --------------------------------------------------------------------------
		// [Factory]
		// --------------------------------------------------------------------------

		std::shared_ptr<AstSymbol> newSymbol(const std::string& key, AstSymbolType symbolType, AstScopeType scopeType);
		std::shared_ptr<AstSymbol> shadowSymbol(const std::shared_ptr<AstSymbol> other);
		void deleteSymbol(std::shared_ptr<AstSymbol> symbol);

#define MATHPRESSO_ALLOC_AST_OBJECT(_Size_) \
  void* obj = _heap->alloc(_Size_); \
  if (MATHPRESSO_UNLIKELY(obj == nullptr)) return nullptr

		template<typename T>
		std::shared_ptr<T> newNode()
		{
			return std::make_shared<T>(shared_from_this());
		}

		template<typename T, typename P0>
		std::shared_ptr<T> newNode(P0 p0)
		{
			return std::make_shared<T>(shared_from_this(), p0);
		}

		template<typename T, typename P0, typename P1>
		std::shared_ptr<T> newNode(P0 p0, P1 p1)
		{
			return std::make_shared<T>(shared_from_this(), p0, p1);
		}

#undef MATHPRESSO_ALLOC_AST_OBJECT

		void deleteNode(std::shared_ptr<AstNode> node);

		uint32_t newSlotId() { return _numSlots++; }

		// --------------------------------------------------------------------------
		// [Init]
		// --------------------------------------------------------------------------

		Error initProgramScope();

		// --------------------------------------------------------------------------
		// [Dump]
		// --------------------------------------------------------------------------

		Error dump(StringBuilder& sb, const Symbols * ops);

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		//! String builder to build possible output messages.
		StringBuilder _sb;

		//! Root node.
		std::shared_ptr<AstProgram> _programNode;

		//! Number of variable slots used.
		uint32_t _numSlots;
	};

	// ============================================================================
	// [mathpresso::AstSymbol]
	// ============================================================================

	struct AstSymbol : public std::enable_shared_from_this<AstSymbol>
	{
		MATHPRESSO_NO_COPY(AstSymbol);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstSymbol(const std::string & name, AstSymbolType symbolType, uint32_t scopeType)
			: _name(name),
			_node(nullptr),
			_symbolType(symbolType),
			_symbolFlags(scopeType == AstScopeType::kAstScopeGlobal ? (int)AstSymbolFlags::kAstSymbolIsGlobal : 0),
			_valueComp(),
			_usedCount(0),
			_writeCount(0)
		{
		}

		~AstSymbol() {}

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		bool eq(const std::string& s) const
		{
			return s == _name;
		}

		//! Get symbol name length.
		size_t getLength() const { return _name.length(); }
		//! Get symbol name.
		std::string getName() const { return _name; }

		//! Check if the symbol has associated node with it.
		bool hasNode() const { return _node != nullptr; }
		//! Get node associated with the symbol (can be `NULL` for built-ins).
		std::shared_ptr<AstNode> getNode() const { return _node; }
		//! Associate node with the symbol (basically the node that declares it).
		void setNode(std::shared_ptr<AstNode> node) { _node = node; }

		//! Get symbol type, see \ref AstSymbolType.
		AstSymbolType getSymbolType() const { return _symbolType; }

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
		void incWriteCount(uint32_t n = 1) { _writeCount += n; incUsedCount(n); }

		void decUsedCount(uint32_t n = 1)
		{
			_usedCount -= n;
			if (_usedCount == 0)
			{
				this->~AstSymbol();
			}
		}
		void decWriteCount(uint32_t n = 1) { _writeCount -= n; decUsedCount(n); }

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

	private:
		//! Symbol name length.
		size_t _length;
		//! Symbol name (key).
		std::string _name;

		//! Node where the symbol is defined.
		std::shared_ptr<AstNode> _node;

		//! Type of the symbol, see \ref AstSymbolType.
		AstSymbolType _symbolType;
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

	// ============================================================================
	// [mathpresso::AstNode]
	// ============================================================================

#define MATHPRESSO_AST_CHILD(_Index_, _Type_, _Name_, _Memb_) \
  bool has##_Name_() const { return _children[_Index_] != nullptr; } \
  std::shared_ptr<_Type_> get##_Name_() const { return _children[_Index_]; } \
  \
  std::shared_ptr<_Type_> set##_Name_(std::shared_ptr<_Type_> node) { \
    _children[_Index_] = node; \
    return std::static_pointer_cast<_Type_>(replaceAt(_Index_, node)); \
  } \
  \
  std::shared_ptr<_Type_> unlink##_Name_() { \
    std::shared_ptr<_Type_> node = _children[_Index_]; \
    \
    MATHPRESSO_ASSERT(node != nullptr); \
    MATHPRESSO_ASSERT(node->getParent() == shared_from_this()); \
    \
    node->_parent.reset(); \
    _children[_Index_] = nullptr; \
    \
    return node; \
  }


	struct AstNode : public std::enable_shared_from_this<AstNode>
	{
		MATHPRESSO_NO_COPY(AstNode);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstNode(std::shared_ptr<AstBuilder> ast, AstNodeType nodeType, std::vector<std::shared_ptr<AstNode>> children = {}, uint32_t length = 0)
			: _ast(ast),
			_parent(),
			_children(children),
			_mpOp(nullptr),
			_opName(),
			_nodeType(nodeType),
			_nodeFlags(AstNodeFlags::kAstNone),
			_position(~static_cast<uint32_t>(0))
		{
		}

		virtual ~AstNode()
		{
		}

		void destroy(std::shared_ptr<AstBuilder> ast)
		{
		}

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		//! Get the `AstBuilder` instance that created this node.
		std::shared_ptr<AstBuilder> getAst() const { return _ast; }

		//! Check if the node has a parent.
		bool hasParent() const { return _parent.lock() != nullptr; }
		//! Get the parent node.
		std::shared_ptr<AstNode> getParent() const { return _parent.lock(); }

		//! Get whether the node has children.
		//!
		//! NOTE: Nodes that always have children (even if they are implicitly set
		//! to NULL) always return `true`. This function if useful mostly if the
		//! node is of `AstBlock` type.
		bool hasChildren() const { return getLength() != 0; }
		//! Get children array.
		std::vector<std::shared_ptr<AstNode>> getChildren() const { return _children; }
		//! Get length of the children array.
		size_t getLength() const { return _children.size(); }

		//! Get node type.
		AstNodeType getNodeType() const { return _nodeType; }
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

		std::shared_ptr<AstNode> getAt(size_t index) const
		{
			MATHPRESSO_ASSERT(index < _children.size());
			return _children[index];
		}

		//! Replace `refNode` by `node`.
		std::shared_ptr<AstNode> replaceNode(std::shared_ptr<AstNode> refNode, std::shared_ptr<AstNode> node);
		//! Replace node at index `index` by `node`.
		std::shared_ptr<AstNode> replaceAt(uint32_t index, std::shared_ptr<AstNode> node);

		//! Inject `node` between this node and `refNode`.
		std::shared_ptr<AstNode> injectNode(std::shared_ptr<AstNode> refNode, std::shared_ptr<AstUnary> node);
		//! Inject `node` between this node and node at index `index`.
		std::shared_ptr<AstNode> injectAt(uint32_t index, std::shared_ptr<AstUnary> node);

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		//! AST builder.
		std::shared_ptr<AstBuilder> _ast; 
		//! Parent node.
		std::weak_ptr<AstNode> _parent;
		//! Child nodes.
		std::vector<std::shared_ptr<AstNode>> _children;

		std::shared_ptr<MpOperation> _mpOp;

		std::string _opName;

	private:
		//! Node type, see `AstNodeType`.
		AstNodeType _nodeType;
		//! Node flags, see `AstNodeFlags`.
		uint8_t _nodeFlags;
		//! Position (in characters) to the beginning of the program (default -1).
		uint32_t _position;
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

		AstBlock(std::shared_ptr<AstBuilder> ast, AstNodeType nodeType = AstNodeType::kAstNodeBlock)
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
		void appendNode(std::shared_ptr<AstNode> node)
		{
			MATHPRESSO_ASSERT(node != nullptr);
			MATHPRESSO_ASSERT(node->getParent() == nullptr);

			// We expect `willAdd()` to be called before `appendNode()`.
			//MATHPRESSO_ASSERT(getLength() < _capacity);

			node->_parent = shared_from_this();

			_children.push_back(node);
		}

		//! Insert the given `node` to the block at index `i`.
		//!
		//! NOTE: You have to call `willAdd()` before you use `insertAt()` for every
		//! node you want to add to the block.
		void insertAt(size_t i, std::shared_ptr<AstNode> node)
		{
			MATHPRESSO_ASSERT(node != nullptr);
			MATHPRESSO_ASSERT(node->getParent() == nullptr);

			// We expect `willAdd()` to be called before `insertAt()`.
			//MATHPRESSO_ASSERT(getLength() < _capacity);

			std::vector<std::shared_ptr<AstNode>> p = getChildren();
			node->_parent = std::static_pointer_cast<AstBlock>(shared_from_this());

			size_t j = getLength();
			while (i < j)
			{
				p[j] = p[j - 1];
				j--;
			}

			p[j] = node;
		}

		//! Remove the given `node`.
		std::shared_ptr<AstNode> removeNode(std::shared_ptr<AstNode> node);
		//! Remove the node at index `index`.
		std::shared_ptr<AstNode> removeAt(size_t index);

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

		AstUnary(std::shared_ptr<AstBuilder> ast, AstNodeType nodeType)
			: AstNode(ast, nodeType, { nullptr }, 1)
		{
		}


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

		AstBinary(std::shared_ptr<AstBuilder> ast, AstNodeType nodeType)
			: AstNode(ast, nodeType, { nullptr, nullptr }, 2)
		{
		}

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

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

		AstTernary(std::shared_ptr<AstBuilder> ast, AstNodeType nodeType)
			: AstNode(ast, nodeType, { nullptr, nullptr, nullptr }, 3)
		{
		}

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

		AstProgram(std::shared_ptr<AstBuilder> ast)
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

		AstVarDecl(std::shared_ptr<AstBuilder> ast)
			: AstUnary(ast, AstNodeType::kAstNodeVarDecl),
			_symbol(nullptr)
		{
		}

		void destroy(std::shared_ptr<AstBuilder> ast)
		{
			std::shared_ptr<AstSymbol> sym = getSymbol();
			if (sym != nullptr)
			{
				sym->decUsedCount();
			}

		}

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		std::shared_ptr<AstSymbol> getSymbol() const { return _symbol; }
		void setSymbol(std::shared_ptr<AstSymbol> symbol) { _symbol = symbol; }

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------
	private:
		std::shared_ptr<AstSymbol> _symbol;

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

		AstVar(std::shared_ptr<AstBuilder> ast)
			: AstNode(ast, AstNodeType::kAstNodeVar),
			_symbol(nullptr)
		{
		}

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		std::shared_ptr<AstSymbol> getSymbol() const { return _symbol; }
		void setSymbol(std::shared_ptr<AstSymbol> symbol) { _symbol = symbol; }

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------
	private:
		std::shared_ptr<AstSymbol> _symbol;
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

		AstImm(std::shared_ptr<AstBuilder> ast, double value = 0.0)
			: AstNode(ast, AstNodeType::kAstNodeImm),
			_value({ value, 0 })
		{
		}

		AstImm(std::shared_ptr<AstBuilder> ast, std::complex<double> value)
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

		AstUnaryOp(std::shared_ptr<AstBuilder> ast)
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

		AstBinaryOp(std::shared_ptr<AstBuilder> ast) : AstBinary(ast, AstNodeType::kAstNodeBinaryOp)
		{
		}

		void destroy(std::shared_ptr<AstBuilder> ast)
		{
			if (_mpOp && (_mpOp->flags() & MpOperation::IsAssignment) && hasLeft())
			{
				std::shared_ptr<AstVar> var = std::static_pointer_cast<AstVar>(getLeft());
				std::shared_ptr<AstSymbol> sym = var->getSymbol();

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

		AstTernaryOp(std::shared_ptr<AstBuilder> ast) :
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

		AstCall(std::shared_ptr<AstBuilder> ast)
			: AstBlock(ast, AstNodeType::kAstNodeCall),
			_symbol(nullptr)
		{
		}

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		std::shared_ptr<AstSymbol> getSymbol() const { return _symbol; }
		void setSymbol(std::shared_ptr<AstSymbol> symbol) { _symbol = symbol; }

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		std::shared_ptr<AstSymbol> _symbol;
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

		AstVisitor(std::shared_ptr<AstBuilder> ast);
		virtual ~AstVisitor();

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		std::shared_ptr<AstBuilder> getAst() const { return _ast; }

		// --------------------------------------------------------------------------
		// [OnNode]
		// --------------------------------------------------------------------------

		virtual Error onNode(std::shared_ptr<AstNode> node);

		virtual Error onProgram(std::shared_ptr<AstProgram> node) ;
		virtual Error onBlock(std::shared_ptr<AstBlock> node) = 0;
		virtual Error onVarDecl(std::shared_ptr<AstVarDecl> node) = 0;
		virtual Error onVar(std::shared_ptr<AstVar> node) = 0;
		virtual Error onImm(std::shared_ptr<AstImm> node) = 0;
		virtual Error onUnaryOp(std::shared_ptr<AstUnaryOp> node) = 0;
		virtual Error onBinaryOp(std::shared_ptr<AstBinaryOp> node) = 0;
		virtual Error onTernaryOp(std::shared_ptr<AstTernaryOp> node) = 0;
		virtual Error onCall(std::shared_ptr<AstCall> node) = 0;

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		std::shared_ptr<AstBuilder> _ast;
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

		AstDump(std::shared_ptr<AstBuilder> ast, StringBuilder& sb, const Symbols * ctx);
		virtual ~AstDump();

		// --------------------------------------------------------------------------
		// [OnNode]
		// --------------------------------------------------------------------------

		virtual Error onBlock(std::shared_ptr<AstBlock> node) override;
		virtual Error onVarDecl(std::shared_ptr<AstVarDecl> node) override;
		virtual Error onVar(std::shared_ptr<AstVar> node) override;
		virtual Error onImm(std::shared_ptr<AstImm> node) override;
		virtual Error onUnaryOp(std::shared_ptr<AstUnaryOp> node) override;
		virtual Error onBinaryOp(std::shared_ptr<AstBinaryOp> node) override;
		virtual Error onTernaryOp(std::shared_ptr<AstTernaryOp> node) override;
		virtual Error onCall(std::shared_ptr<AstCall> node) override;

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
		const Symbols * _ops;
	};

} // mathpresso namespace

// [Guard]
#endif // _MATHPRESSO_MPAST_P_H
