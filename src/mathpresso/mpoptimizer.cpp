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

namespace mathpresso {

	// ============================================================================
	// [mpsl::AstOptimizer - Construction / Destruction]
	// ============================================================================

	AstOptimizer::AstOptimizer(AstBuilder* ast, ErrorReporter* errorReporter, const symbolMap* syms)
		: AstVisitor(ast),
		_errorReporter(errorReporter),
		_symbols(syms) {}
	AstOptimizer::~AstOptimizer() {}

	// ============================================================================
	// [mpsl::AstOptimizer - OnNode]
	// ============================================================================

	Error AstOptimizer::onBlock(AstBlock* node) {
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

			if (curCount < oldCount) {
				if (!alterable)
					return MATHPRESSO_TRACE_ERROR(kErrorInvalidState);
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
			node->addNodeFlags(kAstReturnsComplex | kAstTakesComplex);
		}

		return kErrorOk;
	}

	Error AstOptimizer::callMpOperation(AstVarDecl* node) {
		if (node->mpOp_)
			return node->mpOp_->optimize(this, node);
		return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
			"No MpOperation.");
	}

	Error AstOptimizer::onVar(AstVar* node) {
		AstSymbol* sym = node->getSymbol();
		bool b_complex = node->returnsComplex() || sym->hasSymbolFlag(kAstSymbolIsComplex);

		if (sym->isAssigned() && !node->hasNodeFlag(kAstNodeHasSideEffect)) 
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
		return kErrorOk;
	}

	Error AstOptimizer::onImm(AstImm* node) {
		return kErrorOk;
	}

	Error AstOptimizer::onUnaryOp(AstUnaryOp* node)
	{
		if (node->mpOp_)
		{
			return node->mpOp_->optimize(this, node);
		}
		return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
			"No MpOperation.");
	}

	Error AstOptimizer::onBinaryOp(AstBinaryOp* node)
	{
		if (node->mpOp_ != nullptr)
		{
			return node->mpOp_->optimize(this, node);
		}
		return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
			"No MpOperation.");
	}

	Error AstOptimizer::onTernaryOp(AstTernaryOp* node) {

		if (node->mpOp_ != nullptr)
		{
			return node->mpOp_->optimize(this, node);
		}
		return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
			"No MpOperation.");

	}

	Error AstOptimizer::onCall(AstCall* node) 
	{
		if (node->getSymbol()->getOp() != nullptr)
		{
			return node->getSymbol()->getOp()->optimize(this, node);
		}
		return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
			"No MpOperation.");
	}



} // mathpresso namespace
