// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Export]
#define MATHPRESSO_EXPORTS

// [Dependencies]
#include <mathpresso/mpast_p.h>
#include <mathpresso/mathpresso.h>

namespace mathpresso
{

	// ============================================================================
	// [mathpresso::mpAstNodeSize]
	// ============================================================================

	struct AstNodeSize
	{
		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		uint32_t getNodeType() const { return _nodeType; }
		uint32_t getNodeSize() const { return _nodeSize; }

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		uint8_t _nodeType;
		uint8_t _reserved;
		uint16_t _nodeSize;
	};

#define ROW(type, size) { type, 0, static_cast<uint8_t>(size) }
	static const AstNodeSize mpAstNodeSize[] = {
	  ROW(AstNodeType::kAstNodeNone     , 0),
	  ROW(AstNodeType::kAstNodeProgram  , sizeof(AstProgram)),
	  ROW(AstNodeType::kAstNodeBlock    , sizeof(AstBlock)),
	  ROW(AstNodeType::kAstNodeVarDecl  , sizeof(AstVarDecl)),
	  ROW(AstNodeType::kAstNodeVar      , sizeof(AstVar)),
	  ROW(AstNodeType::kAstNodeImm      , sizeof(AstImm)),
	  ROW(AstNodeType::kAstNodeUnaryOp  , sizeof(AstUnaryOp)),
	  ROW(AstNodeType::kAstNodeBinaryOp , sizeof(AstBinaryOp)),
	  ROW(AstNodeType::kAstNodeTernaryOp, sizeof(AstTernaryOp)),
	  ROW(AstNodeType::kAstNodeCall     , sizeof(AstCall))
	};
#undef ROW

	// ============================================================================
	// [mathpresso::AstBuilder - Construction / Destruction]
	// ============================================================================

	AstBuilder::AstBuilder(ZoneHeap* heap)
		: _heap(heap),
		_rootScope(nullptr),
		_programNode(nullptr),
		_numSlots(0)
	{
	}
	AstBuilder::~AstBuilder()
	{
		if (_programNode)
		{
			deleteNode(_programNode);
		}
	}

	// ============================================================================
	// [mathpresso::AstBuilder - Factory]
	// ============================================================================

	AstScope* AstBuilder::newScope(AstScope* parent, uint32_t scopeType)
	{
		void* p = _heap->alloc(sizeof(AstScope));
		if (p == nullptr)
			return nullptr;
		return new(p) AstScope(this, parent, scopeType);
	}

	void AstBuilder::deleteScope(AstScope* scope)
	{
		scope->~AstScope();
		_heap->release(scope, sizeof(AstScope));
	}

	AstSymbol* AstBuilder::newSymbol(const std::string& key, uint32_t hVal, AstSymbolType symbolType, uint32_t scopeType)
	{
		size_t kLen = key.length();
		void* p = _heap->alloc(sizeof(AstSymbol) + kLen + 1);

		if (p == nullptr)
			return nullptr;

		char* kStr = static_cast<char*>(p) + sizeof(AstSymbol);
		::memcpy(kStr, key.c_str(), kLen);

		kStr[kLen] = '\0';
		return new(p) AstSymbol(kStr, static_cast<uint32_t>(kLen), hVal, symbolType, scopeType);
	}

	AstSymbol* AstBuilder::shadowSymbol(const AstSymbol* other)
	{
		std::string name(other->getName(), other->getLength());
		AstSymbol* sym = newSymbol(name, other->getHVal(), other->getSymbolType(), AstScopeType::kAstScopeShadow);

		if (sym == nullptr)
			return nullptr;

		sym->setSymbolFlag(other->getSymbolFlags());

		if (sym->getSymbolType() == AstSymbolType::kAstSymbolVariable)
		{
			sym->setVarSlotId(other->getVarSlotId());
			sym->setVarOffset(other->getVarOffset());
			sym->setValue(other->getValueComp());
		}

		return sym;
	}

	void AstBuilder::deleteSymbol(AstSymbol* symbol)
	{
		size_t kLen = symbol->getLength();
		symbol->~AstSymbol();
		_heap->release(symbol, sizeof(AstSymbol) + kLen + 1);
	}

	void AstBuilder::deleteNode(AstNode* node)
	{
		size_t length = node->getLength();
		AstNode** children = node->getChildren();

		uint32_t nodeType = node->getNodeType();
		MATHPRESSO_ASSERT(mpAstNodeSize[nodeType].getNodeType() == nodeType);

		switch (nodeType)
		{
			case AstNodeType::kAstNodeProgram: static_cast<AstProgram*>(node)->destroy(this); break;
			case AstNodeType::kAstNodeBlock: static_cast<AstBlock*>(node)->destroy(this); break;
			case AstNodeType::kAstNodeVarDecl: static_cast<AstVarDecl*>(node)->destroy(this); break;
			case AstNodeType::kAstNodeVar: static_cast<AstVar*>(node)->destroy(this); break;
			case AstNodeType::kAstNodeImm: static_cast<AstImm*>(node)->destroy(this); break;
			case AstNodeType::kAstNodeUnaryOp: static_cast<AstUnaryOp*>(node)->destroy(this); break;
			case AstNodeType::kAstNodeBinaryOp: static_cast<AstBinaryOp*>(node)->destroy(this); break;
			case AstNodeType::kAstNodeTernaryOp: static_cast<AstTernaryOp*>(node)->destroy(this); break;
			case AstNodeType::kAstNodeCall: static_cast<AstCall*>(node)->destroy(this); break;
		}

		for (uint32_t i = 0; i < length; i++)
		{
			AstNode* child = children[i];
			if (child != nullptr)
			{
				deleteNode(child);
			}
		}

		node ->~AstNode();

		_heap->release(node, mpAstNodeSize[nodeType].getNodeSize());
	}

	// ============================================================================
	// [mathpresso::AstBuilder - Initialization]
	// ============================================================================

	Error AstBuilder::initProgramScope()
	{
		if (_rootScope == nullptr)
		{
			_rootScope = newScope(nullptr, AstScopeType::kAstScopeGlobal);
			MATHPRESSO_NULLCHECK(_rootScope);
		}

		if (_programNode == nullptr)
		{
			_programNode = newNode<AstProgram>();
			MATHPRESSO_NULLCHECK(_programNode);
		}

		return ErrorCode::kErrorOk;
	}

	// ============================================================================
	// [mathpresso::AstBuilder - Dump]
	// ============================================================================

	Error AstBuilder::dump(StringBuilder& sb, std::shared_ptr<SubContext> ops)
	{
		return AstDump(this, sb, ops).onProgram(getProgramNode());
	}

	// ============================================================================
	// [mathpresso::AstScope - Construction / Destruction]
	// ============================================================================

	struct AstScopeReleaseHandler
	{
		AstScopeReleaseHandler(AstBuilder* ast) : _ast(ast) {}
		void release(AstSymbol* node) { _ast->deleteSymbol(node); }

		AstBuilder* _ast;
	};

	AstScope::AstScope(AstBuilder* ast, AstScope* parent, uint32_t scopeType)
		: _ast(ast),
		_parent(parent),
		_symbols(ast->getHeap()),
		_scopeType(static_cast<uint8_t>(scopeType))
	{
	}

	AstScope::~AstScope()
	{
		AstScopeReleaseHandler handler(_ast);
		_symbols.reset(handler);
	}

	// ============================================================================
	// [mathpresso::AstScope - Ops]
	// ============================================================================

	AstSymbol* AstScope::resolveSymbol(const std::string& name, uint32_t hVal, AstScope** scopeOut)
	{
		AstScope* scope = this;
		AstSymbol* symbol;

		do
		{
			symbol = scope->_symbols.get(name, hVal);
		} while (symbol == nullptr && (scope = scope->getParent()) != nullptr);

		if (scopeOut != nullptr)
			*scopeOut = scope;

		return symbol;
	}

	// ============================================================================
	// [mathpresso::AstNode - Ops]
	// ============================================================================

	AstNode* AstNode::replaceNode(AstNode* refNode, AstNode* node)
	{
		MATHPRESSO_ASSERT(refNode != nullptr);
		MATHPRESSO_ASSERT(refNode->getParent() == this);
		MATHPRESSO_ASSERT(node == nullptr || !node->hasParent());

		size_t length = _length;
		AstNode** children = getChildren();

		for (uint32_t i = 0; i < length; i++)
		{
			AstNode* child = children[i];

			if (child != refNode)
				continue;

			children[i] = node;
			refNode->_parent = nullptr;

			if (node != nullptr)
				node->_parent = this;

			return refNode;
		}

		return nullptr;
	}

	AstNode* AstNode::replaceAt(uint32_t index, AstNode* node)
	{
		AstNode* child = getAt(index);
		_children[index] = node;

		if (child != nullptr)
			child->_parent = nullptr;

		if (node != nullptr)
			node->_parent = this;

		return child;
	}

	AstNode* AstNode::injectNode(AstNode* refNode, AstUnary* node)
	{
		MATHPRESSO_ASSERT(refNode != nullptr && refNode->getParent() == this);
		MATHPRESSO_ASSERT(node != nullptr && node->getParent() == nullptr);

		size_t length = _length;
		AstNode** children = getChildren();

		for (uint32_t i = 0; i < length; i++)
		{
			AstNode* child = children[i];

			if (child != refNode)
				continue;

			children[i] = node;
			refNode->_parent = node;

			node->_parent = this;
			node->setChild(refNode);

			return refNode;
		}

		return nullptr;
	}

	AstNode* AstNode::injectAt(uint32_t index, AstUnary* node)
	{
		AstNode* child = getAt(index);

		MATHPRESSO_ASSERT(node != nullptr && node->getParent() == nullptr);
		MATHPRESSO_ASSERT(child != nullptr);

		_children[index] = node;
		child->_parent = node;

		node->_parent = this;
		node->setChild(child);

		return child;
	}

	template<>
	std::complex<double> AstImm::getValue() const
	{
		return _value;
	}

	template<>
	double AstImm::getValue() const
	{
		return _value.real();
	}

	// ============================================================================
	// [mathpresso::AstBlock - Ops]
	// ============================================================================

	static Error mpBlockNodeGrow(AstBlock* self)
	{
		size_t oldCapacity = self->_capacity;
		size_t newCapacity = oldCapacity;

		size_t length = self->getLength();
		MATHPRESSO_ASSERT(oldCapacity == length);

		// Grow, we prefer growing quickly until we reach 128 and then 1024 nodes. We
		// don't expect to reach these limits in the most used expressions; only test
		// cases can exploit this assumption.
		//
		// Growing schema:
		//   0..4..8..16..32..64..128..256..384..512..640..768..896..1024..[+256]
		if (newCapacity == 0)
			newCapacity = 4;
		else if (newCapacity < 128)
			newCapacity *= 2;
		else if (newCapacity < 1024)
			newCapacity += 128;
		else
			newCapacity += 256;

		ZoneHeap* heap = self->getAst()->getHeap();

		AstNode** oldArray = self->getChildren();
		AstNode** newArray = static_cast<AstNode**>(heap->alloc(newCapacity * sizeof(AstNode), newCapacity));

		MATHPRESSO_NULLCHECK(newArray);
		newCapacity /= sizeof(AstNode*);

		self->_children = newArray;
		self->_capacity = static_cast<uint32_t>(newCapacity);

		if (oldCapacity != 0)
		{
			::memcpy(newArray, oldArray, length * sizeof(AstNode*));
			heap->release(oldArray, oldCapacity * sizeof(AstNode*));
		}

		return ErrorCode::kErrorOk;
	}

	// Tell the AST, we want to add a node, so it can allocate memory if necessary
	Error AstBlock::willAdd()
	{
		// Grow if needed.
		if (_length == _capacity)
			MATHPRESSO_PROPAGATE(mpBlockNodeGrow(this));
		return ErrorCode::kErrorOk;
	}

	AstNode* AstBlock::removeNode(AstNode* node)
	{
		MATHPRESSO_ASSERT(node != nullptr);
		MATHPRESSO_ASSERT(node->getParent() == this);

		AstNode** p = getChildren();
		AstNode** pEnd = p + _length;

		while (p != pEnd)
		{
			if (p[0] == node)
				goto _Found;
			p++;
		}

		// If removeNode() has been called we expect the node to be found. Otherwise
		// there is a bug somewhere.
		MATHPRESSO_ASSERT(!"Reached");
		return nullptr;

	_Found:
		_length--;
		::memmove(p, p + 1, static_cast<size_t>(pEnd - p - 1) * sizeof(AstNode*));

		node->_parent = nullptr;
		return node;
	}

	AstNode* AstBlock::removeAt(size_t index)
	{
		MATHPRESSO_ASSERT(index < _length);

		if (index >= _length)
			return nullptr;

		AstNode** p = getChildren() + index;
		AstNode* oldNode = p[0];

		_length--;
		::memmove(p, p + 1, static_cast<size_t>(_length - index) * sizeof(AstNode*));

		oldNode->_parent = nullptr;
		return oldNode;
	}

	// ============================================================================
	// [mathpresso::AstVisitor - Construction / Destruction]
	// ============================================================================

	AstVisitor::AstVisitor(AstBuilder* ast)
		: _ast(ast)
	{
	}
	AstVisitor::~AstVisitor() {}

	// ============================================================================
	// [mathpresso::AstVisitor - OnNode]
	// ============================================================================

	Error AstVisitor::onNode(AstNode* node)
	{
		switch (node->getNodeType())
		{
			case AstNodeType::kAstNodeProgram: return onProgram(static_cast<AstProgram*>(node));
			case AstNodeType::kAstNodeBlock: return onBlock(static_cast<AstBlock*>(node));
			case AstNodeType::kAstNodeVarDecl: return onVarDecl(static_cast<AstVarDecl*>(node));
			case AstNodeType::kAstNodeVar: return onVar(static_cast<AstVar*>(node));
			case AstNodeType::kAstNodeImm: return onImm(static_cast<AstImm*>(node));
			case AstNodeType::kAstNodeUnaryOp: return onUnaryOp(static_cast<AstUnaryOp*>(node));
			case AstNodeType::kAstNodeBinaryOp: return onBinaryOp(static_cast<AstBinaryOp*>(node));
			case AstNodeType::kAstNodeTernaryOp: return onTernaryOp(static_cast<AstTernaryOp*>(node));
			case AstNodeType::kAstNodeCall: return onCall(static_cast<AstCall*>(node));

			default:
				return MATHPRESSO_TRACE_ERROR(ErrorCode::kErrorInvalidState);
		}
	}

	Error AstVisitor::onProgram(AstProgram* node)
	{
		return onBlock(node);
	}

	// ============================================================================
	// [mathpresso::AstDump - Construction / Destruction]
	// ============================================================================

	AstDump::AstDump(AstBuilder* ast, StringBuilder& sb, std::shared_ptr<SubContext> ops)
		: AstVisitor(ast),
		_sb(sb),
		_level(0),
		_context(ops)
	{
	}
	AstDump::~AstDump() {}

	// ============================================================================
	// [mathpresso::AstDump - OnNode]
	// ============================================================================

	Error AstDump::onBlock(AstBlock* node)
	{
		AstNode** children = node->getChildren();
		size_t i, count = node->getLength();

		for (i = 0; i < count; i++)
			onNode(children[i]);

		return ErrorCode::kErrorOk;
	}

	Error AstDump::onVarDecl(AstVarDecl* node)
	{
		AstSymbol* sym = node->getSymbol();

		nest("%s [VarDecl%s]", sym ? sym->getName() : static_cast<const char*>(nullptr), (node->takesComplex() ? ", complex" : ""));
		if (node->hasChild())
			MATHPRESSO_PROPAGATE(onNode(node->getChild()));
		return denest();
	}

	template<class T>
	std::string sym_name(T * node)
	{
		auto sym = node->getSymbol();
		return sym ? sym->getName() : "(null)";
	}

	std::string op_name(AstNode * node, std::shared_ptr<SubContext> ops)
	{
		return ops->functionName(node->_mpOp);
	}

	const char * node_type(AstNode * node)
	{
		return node->returnsComplex() ? "<cplx>" : "<real>";
	}

	const char * parm_type(AstNode * node)
	{
		return node->takesComplex() ? "<cplx>" : "<real>";
	}

	Error AstDump::onVar(AstVar* node)
	{
		return info("%s %s", sym_name(node).c_str(), node_type(node));
	}

	Error AstDump::onImm(AstImm* node)
	{
		auto v = node->getValue<std::complex<double>>();

		if (node->returnsComplex())
			return info("%lf%+lfi, %s", v.real(), v.imag(), node_type(node));
		else
			return info("%lf, %s", v.real(), node_type(node));
	}


	Error AstDump::onUnaryOp(AstUnaryOp* node)
	{
		nest("%s [Unary, %s -> %s]", op_name(node, _context).c_str(), parm_type(node), node_type(node));
		if (node->hasChild())
			MATHPRESSO_PROPAGATE(onNode(node->getChild()));
		return denest();
	}

	Error AstDump::onBinaryOp(AstBinaryOp* node)
	{
		nest("%s [Binary, %s -> %s]", op_name(node, _context).c_str(), parm_type(node), node_type(node));
		if (node->hasLeft())
			MATHPRESSO_PROPAGATE(onNode(node->getLeft()));
		if (node->hasRight())
			MATHPRESSO_PROPAGATE(onNode(node->getRight()));
		return denest();
	}

	Error AstDump::onTernaryOp(AstTernaryOp* node)
	{
		nest("%s [Ternary, %s -> %s]", op_name(node, _context).c_str(), parm_type(node), node_type(node));
		if (node->hasCondition())
			MATHPRESSO_PROPAGATE(onNode(node->getCondition()));
		if (node->hasLeft())
			MATHPRESSO_PROPAGATE(onNode(node->getLeft()));
		if (node->hasRight())
			MATHPRESSO_PROPAGATE(onNode(node->getRight()));
		return denest();
	}

	Error AstDump::onCall(AstCall* node)
	{
		AstSymbol* sym = node->getSymbol();

		nest("%s(), %s -> %s", node->_opName.c_str(), parm_type(node), node_type(node));
		onBlock(node);
		return denest();
	}

	// ============================================================================
	// [mathpresso::AstDump - Helpers]
	// ============================================================================

	Error AstDump::info(const char* fmt, ...)
	{
		va_list ap;
		va_start(ap, fmt);

		_sb.appendChars(' ', static_cast<size_t>(_level) * 2);
		_sb.appendFormatVA(fmt, ap);
		_sb.appendChar('\n');

		va_end(ap);
		return ErrorCode::kErrorOk;
	}

	Error AstDump::nest(const char* fmt, ...)
	{
		va_list ap;
		va_start(ap, fmt);

		_sb.appendChars(' ', static_cast<size_t>(_level) * 2);
		_sb.appendFormatVA(fmt, ap);
		_sb.appendChar('\n');

		va_end(ap);
		_level++;

		return ErrorCode::kErrorOk;
	}

	Error AstDump::denest()
	{
		MATHPRESSO_ASSERT(_level > 0);
		_level--;

		return ErrorCode::kErrorOk;
	}

} // mathpresso namespace
