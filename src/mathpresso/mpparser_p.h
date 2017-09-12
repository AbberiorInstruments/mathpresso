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

		Parser(std::shared_ptr<AstBuilder> ast, ErrorReporter* errorReporter, const std::string & body, const std::shared_ptr<Context> shadowContext)
			: _ast(ast),
			_errorReporter(errorReporter),
			_tokenizer(body.c_str(), body.length()),
			_shadowContext(shadowContext)
		{
		}
		~Parser() {}

		// --------------------------------------------------------------------------
		// [Parse]
		// --------------------------------------------------------------------------

		Error parseProgram(std::shared_ptr<AstProgram> block);

		Error parseStatement(std::shared_ptr<AstBlock> block, uint32_t flags);
		Error parseBlockOrStatement(std::shared_ptr<AstBlock> block);

		Error parseVariableDecl(std::shared_ptr<AstBlock> block);
		Error parseExpression(std::shared_ptr<AstNode>* pNodeOut, bool isNested);
		Error parseCall(std::shared_ptr<AstNode>* pNodeOut);

		Error reparseTernary(std::shared_ptr<AstNode> node);

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		std::shared_ptr<AstBuilder> _ast;
		ErrorReporter* _errorReporter;

		Tokenizer _tokenizer;

		std::shared_ptr<Context> _shadowContext; // the current context.
	};

} // mathpresso namespace

// [Guard]
#endif // _MATHPRESSO_MPPARSER_P_H
