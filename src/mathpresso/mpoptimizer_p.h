// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MATHPRESSO_MPOPTIMIZER_P_H
#define _MATHPRESSO_MPOPTIMIZER_P_H

// [Dependencies]
#include <mathpresso/mpast_p.h>

namespace mathpresso
{

	// ============================================================================
	// [mpsl::AstOptimizer]
	// ============================================================================

	//! This is a Optimizer for the AST, it uses MpOperation::optimize to optimize
	//! instructions, when an MpOperation is associated with the node and otherwise
	//! optimizes them by itself.
	//! Also it chooses the MpOperation with the correct signature, where possible.
	struct AstOptimizer : public AstVisitor
	{
		MATHPRESSO_NO_COPY(AstOptimizer);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		AstOptimizer(ErrorReporter* errorReporter, std::shared_ptr<Context> ctx) noexcept;
		virtual ~AstOptimizer() noexcept;

		// --------------------------------------------------------------------------
		// [OnNode]
		// --------------------------------------------------------------------------

		Error onNode(std::shared_ptr<AstNode> node) override;

		Error onBlock(std::shared_ptr<AstBlock> node) override;
		Error onVarDecl(std::shared_ptr<AstVarDecl> node) override;
		Error onVar(std::shared_ptr<AstVar> node) override;
		Error onImm(std::shared_ptr<AstImm> node) override;
		Error onUnaryOp(std::shared_ptr<AstUnaryOp> node) override;
		Error onBinaryOp(std::shared_ptr<AstBinaryOp> node) override;
		Error onTernaryOp(std::shared_ptr<AstTernaryOp> node) override;
		Error onCall(std::shared_ptr<AstCall> node) override;

		Error optimize(std::shared_ptr<AstNode> node);

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

	private:

		ErrorReporter* _errorReporter;
		std::shared_ptr<Context> _shadowContext;

	};

} // mathpresso namespace

// [Guard]
#endif // _MATHPRESSO_MPOPTIMIZER_P_H
