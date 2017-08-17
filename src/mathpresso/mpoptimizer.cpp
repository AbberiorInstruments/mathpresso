// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Export]
#define MATHPRESSO_EXPORTS

// [Dependencies]
#include <mathpresso/mpast_p.h>
#include <mathpresso/mpeval_p.h>
#include <mathpresso/mpoptimizer_p.h>
#include <mathpresso/mpoperation.h>

namespace mathpresso
{

	// ============================================================================
	// [mpsl::AstOptimizer - Construction / Destruction]
	// ============================================================================

	AstOptimizer::AstOptimizer(AstBuilder* ast, ErrorReporter* errorReporter, const Operations* ops)
		: AstVisitor(ast),
		_errorReporter(errorReporter),
		_ops(ops)
	{
	}
	AstOptimizer::~AstOptimizer() {}

	// ============================================================================
	// [mpsl::AstOptimizer - OnNode]
	// ============================================================================

	Error AstOptimizer::onNode(AstNode * node)
	{
		switch (node->getNodeType())
		{
			case AstNodeType::kAstNodeProgram: return onProgram(static_cast<AstProgram*>(node));
			case AstNodeType::kAstNodeBlock: return onBlock(static_cast<AstBlock*>(node));
			case AstNodeType::kAstNodeVar: return onVar(static_cast<AstVar*>(node));
			case AstNodeType::kAstNodeImm: return onImm(static_cast<AstImm*>(node));
			case AstNodeType::kAstNodeVarDecl:
			case AstNodeType::kAstNodeUnaryOp:
			case AstNodeType::kAstNodeBinaryOp:
			case AstNodeType::kAstNodeTernaryOp:
			case AstNodeType::kAstNodeCall:
				return optimize(node);

			default:
				return MATHPRESSO_TRACE_ERROR(ErrorCode::kErrorInvalidState);
		}
	}

	Error AstOptimizer::onBlock(AstBlock* node)
	{
		// Prevent removing nodes that are not stored in pure `AstBlock`. For example
		// function call inherits from `AstBlock`, but it needs each expression passed.
		bool alterable = node->getNodeType() == kAstNodeBlock;

		size_t curCount = node->getLength();
		size_t oldCount;
		bool isComplex = false;

		size_t i = 0;
		while (i < curCount)
		{
			MATHPRESSO_PROPAGATE(onNode(node->getAt(i)));

			oldCount = curCount;
			curCount = node->getLength();

			if (curCount < oldCount)
			{
				if (!alterable)
					return MATHPRESSO_TRACE_ERROR(ErrorCode::kErrorInvalidState);
				continue;
			}

			if (alterable && (node->getAt(i)->isImm()))
			{
				_ast->deleteNode(node->removeAt(i));
				curCount--;
				continue;
			}
			// TODO: check semantics here
			isComplex |= node->getAt(i)->returnsComplex();
			i++;
		}

		if (isComplex)
		{
			node->addNodeFlags(AstNodeFlags::kAstReturnsComplex | AstNodeFlags::kAstTakesComplex);
		}

		return ErrorCode::kErrorOk;
	}

	Error AstOptimizer::onVarDecl(AstVarDecl* node)
	{
		if (node->_mpOp)
			return node->_mpOp->optimize(this, node);
		return _errorReporter->onError(ErrorCode::kErrorInvalidState, node->getPosition(),
									   "No MpOperation.");
	}

	Error AstOptimizer::onVar(AstVar* node)
	{
		AstSymbol* sym = node->getSymbol();
		bool b_complex = node->returnsComplex() || sym->hasSymbolFlag(AstSymbolFlags::kAstSymbolIsComplex);

		if (sym->isAssigned() && !node->hasNodeFlag(AstNodeFlags::kAstNodeHasSideEffect))
		{
			AstImm* imm;

			if (!b_complex)
			{
				imm = _ast->newNode<AstImm>(sym->getValue());
			}
			else
			{
				imm = _ast->newNode<AstImm>(sym->getValueComp());
			}
			_ast->deleteNode(node->getParent()->replaceNode(node, imm));
		}
		return ErrorCode::kErrorOk;
	}

	Error AstOptimizer::onImm(AstImm* node)
	{
		return ErrorCode::kErrorOk;
	}

	Error AstOptimizer::optimize(AstNode * node)
	{

		bool takesComplex = false;

		for (size_t i = 0; i < node->getLength(); i++)
		{
			MATHPRESSO_PROPAGATE(onNode(node->getAt(i)));
			takesComplex |= node->getAt(i)->returnsComplex();
		}

		if (node->getNodeType() == AstNodeType::kAstNodeTernaryOp)
		{
			takesComplex = node->getAt(1)->returnsComplex() || node->getAt(2)->returnsComplex();
		}

		// Find optimal signature
		node->_mpOp = _ops->find(node->_opName, node->getLength(), takesComplex);

		if (node->_mpOp)
		{
			if (node->_mpOp->signature().return_type_ == Signature::type::complex)
			{
				node->addNodeFlags(AstNodeFlags::kAstReturnsComplex);
			}

			if (node->_mpOp->signature().areParams(Signature::type::complex))
			{
				node->addNodeFlags(AstNodeFlags::kAstTakesComplex);
			}

			return node->_mpOp->optimize(this, node);
		}
		return _errorReporter->onError(ErrorCode::kErrorInvalidState, node->getPosition(),
									   "No MpOperation.");
	}

	Error AstOptimizer::onUnaryOp(AstUnaryOp* node)
	{
		return optimize(node);
	}

	Error AstOptimizer::onBinaryOp(AstBinaryOp* node)
	{
		return optimize(node);
	}

	Error AstOptimizer::onTernaryOp(AstTernaryOp* node)
	{
		return optimize(node);
	}

	Error AstOptimizer::onCall(AstCall* node)
	{
		return optimize(node);
	}


} // mathpresso namespace
