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

  AstOptimizer(AstBuilder* ast, ErrorReporter* errorReporter, std::shared_ptr<SubContext> ops);
  virtual ~AstOptimizer();

  // --------------------------------------------------------------------------
  // [OnNode]
  // --------------------------------------------------------------------------

  virtual Error onNode(AstNode * node) override;

  virtual Error onBlock(AstBlock* node);
  virtual Error onVarDecl(AstVarDecl* node);
  virtual Error onVar(AstVar* node);
  virtual Error onImm(AstImm* node);
  virtual Error optimize(AstNode * node);
  virtual Error onUnaryOp(AstUnaryOp* node);
  virtual Error onBinaryOp(AstBinaryOp* node);
  virtual Error onTernaryOp(AstTernaryOp * node);
  virtual Error onCall(AstCall* node);

  // --------------------------------------------------------------------------
  // [Members]
  // --------------------------------------------------------------------------

  ErrorReporter* _errorReporter;
  std::shared_ptr<SubContext> _context;
};

} // mathpresso namespace

// [Guard]
#endif // _MATHPRESSO_MPOPTIMIZER_P_H
