// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MATHPRESSO_MPTOKENIZER_P_H
#define _MATHPRESSO_MPTOKENIZER_P_H

// [Dependencies]
#include <mathpresso/mathpresso_p.h>
#include <mathpresso/mpstrtod_p.h>
#include <complex>

namespace mathpresso
{

	// ============================================================================
	// [mathpresso::Token]
	// ============================================================================


	//! \internal
	//!
	//! Token type.
	enum TokenType
	{
		kTokenInvalid = 0,  // <invalid>

		kTokenSymbol,       // <symbol>
		kTokenNumber,       // <number>
		kTokenComplex,	  // <complex number>, written as a+bi

		kTokenVar,          // 'var' keyword
		kTokenReserved,     // reserved keyword

		kTokenDot = 36,     // .
		kTokenComma,        // ,
		kTokenSemicolon,    // ;

		kTokenLCurl,        // {
		kTokenRCurl,        // }

		// unused:
		kTokenLBracket,     // [
		// unused:
		kTokenRBracket,     // ]

		kTokenLParen,       // (
		kTokenRParen,       // )

		kTokenOperator,

		kTokenEnd           // <end>
	};
	//! \internal
	//!
	//! Character classes used by tokenizer.
	enum TokenChar
	{
		// Digit.
		kTokenChar0x0, kTokenChar0x1, kTokenChar0x2, kTokenChar0x3,
		kTokenChar0x4, kTokenChar0x5, kTokenChar0x6, kTokenChar0x7,
		kTokenChar0x8, kTokenChar0x9,

		// Digit-Hex.
		kTokenChar0xA, kTokenChar0xB, kTokenChar0xC, kTokenChar0xD,
		kTokenChar0xE, kTokenChar0xF,

		// imaginary
		kTokenCharImg,

		// Non-Hex ASCII [A-Z] Letter and Underscore [_].
		kTokenCharSym,

		// Punctuation.
		kTokenCharDot = TokenType::kTokenDot,          // .
		kTokenCharCom = TokenType::kTokenComma,        // ,
		kTokenCharSem = TokenType::kTokenSemicolon,    // ;
		kTokenCharLCu = TokenType::kTokenLCurl,        // {
		kTokenCharRCu = TokenType::kTokenRCurl,        // }
		kTokenCharLBr = TokenType::kTokenLBracket,     // [
		kTokenCharRBr = TokenType::kTokenRBracket,     // ]
		kTokenCharLPa = TokenType::kTokenLParen,       // (
		kTokenCharRPa = TokenType::kTokenRParen,       // )

		// Marks a Operator
		kTokenCharOp = TokenType::kTokenOperator,

		// Space.
		kTokenCharSpc = 63,

		// Extended ASCII character (0x80 and above), acts as non-recognized.
		kTokenCharExt,
		// Invalid (non-recognized) character.
		kTokenCharInv,

		kTokenCharSingleCharTokenEnd = kTokenCharRPa
	};

#define C(_Id_) kTokenChar##_Id_
	static const uint8_t mpCharClass[] = {
		C(Inv), C(Inv), C(Inv), C(Inv), C(Inv), C(Inv), C(Inv), C(Inv), // 000-007 ........ | All invalid.
		C(Inv), C(Spc), C(Spc), C(Spc), C(Spc), C(Spc), C(Inv), C(Inv), // 008-015 .     .. | Spaces: 0x9-0xD.
		C(Inv), C(Inv), C(Inv), C(Inv), C(Inv), C(Inv), C(Inv), C(Inv), // 016-023 ........ | All invalid.
		C(Inv), C(Inv), C(Inv), C(Inv), C(Inv), C(Inv), C(Inv), C(Inv), // 024-031 ........ | All invalid.
		C(Spc), C(Op), C(Inv), C(Inv), C(Inv), C(Op), C(Op), C(Inv), // 032-039  !"#$%&' | Unassigned: "#$'.
		C(LPa), C(RPa), C(Op), C(Op), C(Com), C(Op), C(Dot), C(Op), // 040-047 ()*+,-./ |
		C(0x0), C(0x1), C(0x2), C(0x3), C(0x4), C(0x5), C(0x6), C(0x7), // 048-055 01234567 |
		C(0x8), C(0x9), C(Op), C(Sem), C(Op),  C(Op),  C(Op),  C(Op), // 056-063 89:;<=>? |
		C(Inv), C(0xA), C(0xB), C(0xC), C(0xD), C(0xE), C(0xF), C(Sym), // 064-071 @ABCDEFG | Unassigned: @.
		C(Sym), C(Sym), C(Sym), C(Sym), C(Sym), C(Sym), C(Sym), C(Sym), // 072-079 HIJKLMNO |
		C(Sym), C(Sym), C(Sym), C(Sym), C(Sym), C(Sym), C(Sym), C(Sym), // 080-087 PQRSTUVW |
		C(Sym), C(Sym), C(Sym), C(LBr), C(Inv), C(RBr), C(Op), C(Sym), // 088-095 XYZ[\]^_ | Unassigned: \.
		C(Inv), C(0xA), C(0xB), C(0xC), C(0xD), C(0xE), C(0xF), C(Sym), // 096-103 `abcdefg | Unassigned: `.
		C(Sym), C(Img), C(Sym), C(Sym), C(Sym), C(Sym), C(Sym), C(Sym), // 104-111 hijklmno |
		C(Sym), C(Sym), C(Sym), C(Sym), C(Sym), C(Sym), C(Sym), C(Sym), // 112-119 pqrstuvw |
		C(Sym), C(Sym), C(Sym), C(LCu), C(Op),  C(RCu), C(Op), C(Inv), // 120-127 xyz{|}~  |
		C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), // 128-135 ........ | Extended.
		C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), // 136-143 ........ |
		C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), // 144-151 ........ |
		C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), // 152-159 ........ |
		C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), // 160-167 ........ |
		C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), // 168-175 ........ |
		C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), // 176-183 ........ |
		C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), // 184-191 ........ |
		C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), // 192-199 ........ |
		C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), // 200-207 ........ |
		C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), // 208-215 ........ |
		C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), // 216-223 ........ |
		C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), // 224-231 ........ |
		C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), // 232-239 ........ |
		C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), // 240-247 ........ |
		C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext), C(Ext)  // 248-255 ........ |
	};
#undef C


	// ============================================================================
	// [mathpresso::Token]
	// ============================================================================

	//! \internal
	//!
	//! Token.
	struct Token
	{
		// --------------------------------------------------------------------------
		// [Reset]
		// --------------------------------------------------------------------------

		void reset()
		{
			position = 0;
			length = 0;
			value = 0.0;
			token = TokenType::kTokenInvalid;
		}

		// --------------------------------------------------------------------------
		// [Accessors]
		// --------------------------------------------------------------------------

		uint32_t setData(size_t position, size_t length, uint32_t token)
		{
			this->position = position;
			this->length = length;
			this->token = token;
			return token;
		}

		uint32_t getPosAsUInt() const
		{
			MATHPRESSO_ASSERT(position < ~static_cast<uint32_t>(0));
			return static_cast<uint32_t>(position);
		}

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		//! Token position from the beginning of the input.
		size_t position;
		//! Token string length.
		size_t length;

		//! Token type.
		uint32_t token;

		//! Token value (if the token is a number). If token = TokenType::kTokenComplex this 
		// should be interpreted as the imaginary part of a complex number.
		double value;

	};

	// ============================================================================
	// [mathpresso::Tokenizer]
	// ============================================================================

	struct Tokenizer
	{
		MATHPRESSO_NO_COPY(Tokenizer);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		Tokenizer(const char* s, size_t sLen)
			: _p(s),
			_start(s),
			_end(s + sLen),
			_strtod()
		{
			_token.reset();
		}

		// --------------------------------------------------------------------------
		// [Ops]
		// --------------------------------------------------------------------------

		//! Get the current token.
		uint32_t peek(Token* token);
		//! Get the current token and advance.
		uint32_t next(Token* token);

		//! Set the token that will be returned by `next()` and `peek()` functions.
		void set(Token* token)
		{
			// We have to update also _p in case that multiple tokens were put back.
			_p = _start + token->position + token->length;
			_token = *token;
		}

		//! Consume a token got by using peek().
		void consume()
		{
			_token.token = TokenType::kTokenInvalid;
		}

		//! Consume a token got by using peek() and call `peek()`.
		//!
		//! NOTE: Can be called only immediately after peek().
		uint32_t consumeAndPeek(Token* token)
		{
			consume();
			return peek(token);
		}

		//! Consume a token got by using peek() and call `next()`.
		//!
		//! NOTE: Can be called only immediately after peek().
		uint32_t consumeAndNext(Token* token)
		{
			consume();
			return next(token);
		}

		std::string getTokenName(Token token) const
		{
			return std::string(_start + token.position, token.length);
		}

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------
		const char* _start;
	
	private:
		const char* _p;
		const char* _end;

		StrToD _strtod;
		Token _token;
	};

} // mathpresso namespace

// [Guard]
#endif // _MATHPRESSO_MPTOKENIZER_P_H
