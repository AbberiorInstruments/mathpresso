// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MATHPRESSO_MPPARSER_P_H
#define _MATHPRESSO_MPPARSER_P_H

// [Dependencies]
#include <mathpresso/mpast_p.h>
#include <mathpresso/mptokenizer_p.h>

namespace mathpresso
{

	// ============================================================================
	// [Forward Declaration]
	// ============================================================================

	struct AstBuilder;
	struct AstBlock;
	struct AstNode;
	struct AstProgram;
	struct AstVar;

	// ============================================================================
	// [mathpresso::Parser]
	// ============================================================================

	struct Parser
	{
		MATHPRESSO_NO_COPY(Parser);

		enum ParserFlags
		{
			kNoFlags = 0x00,
			kEnableVarDecls = 0x01,
			kEnableNestedBlock = 0x02
		};

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		Parser(AstBuilder* ast, ErrorReporter* errorReporter, const char* body, size_t len, const Operations * ops)
			: _ast(ast),
			_errorReporter(errorReporter),
			_currentScope(ast->getRootScope()),
			_tokenizer(body, len),
			_ops(ops)
		{
		}
		~Parser() {}

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		AstScope* getCurrentScope() const { return _currentScope; }

		// --------------------------------------------------------------------------
		// [Parse]
		// --------------------------------------------------------------------------

		Error parseProgram(AstProgram* block);

		Error parseStatement(AstBlock* block, uint32_t flags);
		Error parseBlockOrStatement(AstBlock* block);

		Error parseVariableDecl(AstBlock* block);
		Error parseExpression(AstNode** pNodeOut, bool isNested);
		Error parseCall(AstNode** pNodeOut);

		Error reparseTernary(AstNode* node);

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		AstBuilder* _ast;
		ErrorReporter* _errorReporter;

		AstScope* _currentScope;
		Tokenizer _tokenizer;
		const Operations * _ops;
	};

} // mathpresso namespace

// [Guard]
#endif // _MATHPRESSO_MPPARSER_P_H
