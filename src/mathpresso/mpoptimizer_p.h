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

namespace mathpresso {

// ============================================================================
// [mpsl::AstOptimizer]
// ============================================================================

struct AstOptimizer : public AstVisitor {
  MATHPRESSO_NO_COPY(AstOptimizer)

  // --------------------------------------------------------------------------
  // [Construction / Destruction]
  // --------------------------------------------------------------------------

  AstOptimizer(std::shared_ptr<AstBuilder> ast, ErrorReporter* errorReporter, std::shared_ptr<Context> ctx);
  virtual ~AstOptimizer();

  // --------------------------------------------------------------------------
  // [OnNode]
  // --------------------------------------------------------------------------

  virtual Error onNode(std::shared_ptr<AstNode> node) override;

  virtual Error onBlock(std::shared_ptr<AstBlock> node) override;
  virtual Error onVarDecl(std::shared_ptr<AstVarDecl> node) override;
  virtual Error onVar(std::shared_ptr<AstVar> node) override;
  virtual Error onImm(std::shared_ptr<AstImm> node) override;
  virtual Error onUnaryOp(std::shared_ptr<AstUnaryOp> node) override;
  virtual Error onBinaryOp(std::shared_ptr<AstBinaryOp> node) override;
  virtual Error onTernaryOp(std::shared_ptr<AstTernaryOp> node) override;
  virtual Error onCall(std::shared_ptr<AstCall> node) override;

  virtual Error optimize(std::shared_ptr<AstNode> node);

  // --------------------------------------------------------------------------
  // [Members]
  // --------------------------------------------------------------------------

  ErrorReporter* _errorReporter;
  std::shared_ptr<Context> _shadowContext;
};

} // mathpresso namespace

// [Guard]
#endif // _MATHPRESSO_MPOPTIMIZER_P_H
