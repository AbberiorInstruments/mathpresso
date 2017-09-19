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


	static uint32_t getNodeSize(AstNodeType nodeType)
	{
		switch (nodeType)
		{
			case AstNodeType::kAstNodeNone: return 0;
			case AstNodeType::kAstNodeProgram: return sizeof(AstProgram);
			case AstNodeType::kAstNodeBlock: return sizeof(AstBlock);
			case AstNodeType::kAstNodeVarDecl: return sizeof(AstVarDecl);
			case AstNodeType::kAstNodeVar: return sizeof(AstVar);
			case AstNodeType::kAstNodeImm: return sizeof(AstImm);
			case AstNodeType::kAstNodeUnaryOp: return sizeof(AstUnaryOp);
			case AstNodeType::kAstNodeBinaryOp: return sizeof(AstBinaryOp);
			case AstNodeType::kAstNodeTernaryOp: return sizeof(AstTernaryOp);
			case AstNodeType::kAstNodeCall: return sizeof(AstCall);
			default: return 0;
		}
	}

	// ============================================================================
	// [mathpresso::AstBuilder - Construction / Destruction]
	// ============================================================================

	AstBuilder::AstBuilder()
		: _programNode(nullptr),
		_numSlots(0)
	{
	}
	AstBuilder::~AstBuilder()
	{
	}

	// ============================================================================
	// [mathpresso::AstBuilder - Factory]
	// ============================================================================

	std::shared_ptr<AstSymbol> AstBuilder::shadowSymbol(const std::shared_ptr<AstSymbol> other)
	{
		std::string name(other->getName());
		std::shared_ptr<AstSymbol> sym = std::make_shared<AstSymbol>(name, other->getSymbolType(), false);

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

	// ============================================================================
	// [mathpresso::AstBuilder - Initialization]
	// ============================================================================

	Error AstBuilder::initProgramScope()
	{
		if (_programNode == nullptr)
		{
			_programNode = std::make_shared<AstProgram>();
			MATHPRESSO_NULLCHECK(_programNode);
		}

		return ErrorCode::kErrorOk;
	}

	// ============================================================================
	// [mathpresso::AstBuilder - Dump]
	// ============================================================================

	Error AstBuilder::dump(StringBuilder& sb, const std::shared_ptr<Symbols> syms)
	{
		return AstDump(sb, syms).onProgram(programNode());
	}

	// ============================================================================
	// [mathpresso::AstNode - Ops]
	// ============================================================================

	std::shared_ptr<AstNode> AstNode::replaceNode(std::shared_ptr<AstNode> refNode, std::shared_ptr<AstNode> node)
	{
		MATHPRESSO_ASSERT(refNode != nullptr);
		MATHPRESSO_ASSERT(refNode->getParent() == shared_from_this());
		MATHPRESSO_ASSERT(node == nullptr || !node->hasParent());

		size_t length = getLength();
		std::vector<std::shared_ptr<AstNode>> children = getChildren();

		for (uint32_t i = 0; i < length; i++)
		{
			std::shared_ptr<AstNode> child = children[i];

			if (child != refNode)
				continue;

			return replaceAt(i, node);
		}

		return nullptr;
	}

	std::shared_ptr<AstNode> AstNode::replaceAt(uint32_t index, std::shared_ptr<AstNode> node)
	{
		std::shared_ptr<AstNode> child = getAt(index);
		_children[index] = node;

		if (child != nullptr)
			child->_parent.reset();

		if (node != nullptr)
			node->_parent = shared_from_this();

		return child;
	}

	std::shared_ptr<AstNode> AstNode::injectNode(std::shared_ptr<AstNode> refNode, std::shared_ptr<AstUnary> node)
	{
		MATHPRESSO_ASSERT(refNode != nullptr && refNode->getParent() == shared_from_this());
		MATHPRESSO_ASSERT(node != nullptr && node->getParent() == nullptr);

		size_t length = getLength();
		std::vector<std::shared_ptr<AstNode>> children = getChildren();

		for (uint32_t i = 0; i < length; i++)
		{
			std::shared_ptr<AstNode> child = children[i];

			if (child != refNode)
				continue;

			children[i] = node;
			refNode->_parent = node;

			node->_parent = shared_from_this();
			node->setChild(refNode);

			return refNode;
		}

		return nullptr;
	}

	std::shared_ptr<AstNode> AstNode::injectAt(uint32_t index, std::shared_ptr<AstUnary> node)
	{
		std::shared_ptr<AstNode> child = getAt(index);

		MATHPRESSO_ASSERT(node != nullptr && node->getParent() == nullptr);
		MATHPRESSO_ASSERT(child != nullptr);

		_children[index] = node;
		child->_parent = node;

		node->_parent = shared_from_this();
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

	std::shared_ptr<AstNode> AstBlock::removeNode(std::shared_ptr<AstNode> node)
	{
		MATHPRESSO_ASSERT(node != nullptr);
		MATHPRESSO_ASSERT(node->getParent() == shared_from_this());

		auto p = _children.begin();
		auto pEnd = _children.end();

		while (p != pEnd)
		{
			if (p[0] == node)
			{
				_children.erase(p);
				node->_parent.reset();
				return node;
			}
			p++;
		}

		// If removeNode() has been called we expect the node to be found. Otherwise
		// there is a bug somewhere.
		MATHPRESSO_ASSERT(!"Reached");
		return nullptr;
	}

	std::shared_ptr<AstNode> AstBlock::removeAt(size_t index)
	{
		MATHPRESSO_ASSERT(index < getLength());

		if (index >= getLength())
			return nullptr;


		auto oldNode = _children[index];
		_children.erase(_children.begin() + index);

		oldNode->_parent.reset();

		return oldNode;
	}

	// ============================================================================
	// [mathpresso::AstVisitor - Construction / Destruction]
	// ============================================================================

	AstVisitor::AstVisitor()
	{
	}
	AstVisitor::~AstVisitor() {}

	// ============================================================================
	// [mathpresso::AstVisitor - OnNode]
	// ============================================================================

	Error AstVisitor::onNode(std::shared_ptr<AstNode> node)
	{
		switch (node->getNodeType())
		{
			case AstNodeType::kAstNodeProgram: return onProgram(std::static_pointer_cast<AstProgram>(node));
			case AstNodeType::kAstNodeBlock: return onBlock(std::static_pointer_cast<AstBlock>(node));
			case AstNodeType::kAstNodeVarDecl: return onVarDecl(std::static_pointer_cast<AstVarDecl>(node));
			case AstNodeType::kAstNodeVar: return onVar(std::static_pointer_cast<AstVar>(node));
			case AstNodeType::kAstNodeImm: return onImm(std::static_pointer_cast<AstImm>(node));
			case AstNodeType::kAstNodeUnaryOp: return onUnaryOp(std::static_pointer_cast<AstUnaryOp>(node));
			case AstNodeType::kAstNodeBinaryOp: return onBinaryOp(std::static_pointer_cast<AstBinaryOp>(node));
			case AstNodeType::kAstNodeTernaryOp: return onTernaryOp(std::static_pointer_cast<AstTernaryOp>(node));
			case AstNodeType::kAstNodeCall: return onCall(std::static_pointer_cast<AstCall>(node));

			default:
				return MATHPRESSO_TRACE_ERROR(ErrorCode::kErrorInvalidState);
		}
	}

	Error AstVisitor::onProgram(std::shared_ptr<AstProgram> node)
	{
		return onBlock(node);
	}

	// ============================================================================
	// [mathpresso::AstDump - Construction / Destruction]
	// ============================================================================

	AstDump::AstDump(StringBuilder& sb, const std::shared_ptr<const Symbols> syms)
		: AstVisitor(),
		_sb(sb),
		_level(0),
		_symbols(syms)
	{
	}
	AstDump::~AstDump() {}

	// ============================================================================
	// [mathpresso::AstDump - OnNode]
	// ============================================================================

	Error AstDump::onBlock(std::shared_ptr<AstBlock> node)
	{
		auto children = node->getChildren();
		size_t i, count = node->getLength();

		for (i = 0; i < count; i++)
			onNode(children[i]);

		return ErrorCode::kErrorOk;
	}

	Error AstDump::onVarDecl(std::shared_ptr<AstVarDecl> node)
	{
		std::shared_ptr<AstSymbol> sym = node->getSymbol();

		nest("%s [VarDecl%s]", sym ? sym->getName().c_str() : "", node->takesComplex() ? ", complex" : ", real");
		if (node->hasChild())
			MATHPRESSO_PROPAGATE(onNode(node->getChild()));
		return denest();
	}

	template<class T>
	std::string sym_name(std::shared_ptr<T> node)
	{
		auto sym = node->getSymbol();
		return sym ? sym->getName() : "(null)";
	}

	std::string op_name(std::shared_ptr<AstNode> node, const std::shared_ptr<const Symbols> ops)
	{
		return ops->name(node->_mpOp);
	}

	const char * node_type(std::shared_ptr<AstNode> node)
	{
		return node->returnsComplex() ? "<cplx>" : "<real>";
	}

	const char * parm_type(std::shared_ptr<AstNode> node)
	{
		return node->takesComplex() ? "<cplx>" : "<real>";
	}

	Error AstDump::onVar(std::shared_ptr<AstVar> node)
	{
		return info("%s %s", sym_name(node).c_str(), node_type(node));
	}

	Error AstDump::onImm(std::shared_ptr<AstImm> node)
	{
		auto v = node->getValue<std::complex<double>>();

		if (node->returnsComplex())
			return info("%lf%+lfi, %s", v.real(), v.imag(), node_type(node));
		else
			return info("%lf, %s", v.real(), node_type(node));
	}


	Error AstDump::onUnaryOp(std::shared_ptr<AstUnaryOp> node)
	{
		nest("%s [Unary, %s -> %s]", op_name(node, _symbols).c_str(), parm_type(node), node_type(node));
		if (node->hasChild())
			MATHPRESSO_PROPAGATE(onNode(node->getChild()));
		return denest();
	}

	Error AstDump::onBinaryOp(std::shared_ptr<AstBinaryOp> node)
	{
		nest("%s [Binary, %s -> %s]", op_name(node, _symbols).c_str(), parm_type(node), node_type(node));
		if (node->hasLeft())
			MATHPRESSO_PROPAGATE(onNode(node->getLeft()));
		if (node->hasRight())
			MATHPRESSO_PROPAGATE(onNode(node->getRight()));
		return denest();
	}

	Error AstDump::onTernaryOp(std::shared_ptr<AstTernaryOp> node)
	{
		nest("%s [Ternary, %s -> %s]", op_name(node, _symbols).c_str(), parm_type(node), node_type(node));
		if (node->hasCondition())
			MATHPRESSO_PROPAGATE(onNode(node->getCondition()));
		if (node->hasLeft())
			MATHPRESSO_PROPAGATE(onNode(node->getLeft()));
		if (node->hasRight())
			MATHPRESSO_PROPAGATE(onNode(node->getRight()));
		return denest();
	}

	Error AstDump::onCall(std::shared_ptr<AstCall> node)
	{
		std::shared_ptr<AstSymbol> sym = node->getSymbol();

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
