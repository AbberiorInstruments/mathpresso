// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Export]
#define MATHPRESSO_EXPORTS

// [Dependencies]
#include <mathpresso/mpeval_p.h>
#include <mathpresso/mphash_p.h>
#include <mathpresso/mptokenizer_p.h>
#include <iostream>

namespace mathpresso {

	// ============================================================================
	// [mathpresso::Tokenizer]
	// ============================================================================

	bool isIn(char c, std::string arr)
	{
		for (auto p : arr)
		{
			if (p == tolower(c))
			{
				return true;
			}
		}
		return false;
	}

	bool isOperator(char c)
	{
		return isIn(c, "+-*/^<>=~%!&|?:");
	}

	bool isNum(char c)
	{
		return isIn(c, "1234567890.");
	}

	bool isSymbolFirst(char c)
	{
		return isIn(c, "abcdefghijklmnopqrstuvwxyz_");
	}

	bool isSymbol(char c)
	{
		return isSymbolFirst(c) || isNum(c);
	}

	bool isSeparator(char c)
	{
		return isIn(c, ",()[]{};");
	}

	bool isSpace(char c)
	{
		return  (c == ' ' || c == '\t');
	}



	//! \internal
	//!
	//! RAW lowercase conversion.
	//!
	//! This method exploits how ASCII table has been designed. It expects ASCII
	//! character on the input that will be lowercased by setting the 0x20 bit on.
	static uint32_t mpGetLower(uint32_t c) 
	{
		return c | 0x20;
	}

	//! \internal
	static constexpr double mpPow10Table[] =
	{
	  1e+0 , 1e+1 , 1e+2 , 1e+3 , 1e+4 , 1e+5 , 1e+6 , 1e+7 ,
	  1e+8 , 1e+9 , 1e+10, 1e+11, 1e+12, 1e+13, 1e+14, 1e+15
	};

	//! \internal
	enum
	{
		kSafeDigits = 15,
		kPow10TableSize = static_cast<int>(MATHPRESSO_ARRAY_SIZE(mpPow10Table))
	};

#define CHAR4X(C0, C1, C2, C3) \
  ( (static_cast<uint32_t>(C0)      ) + \
    (static_cast<uint32_t>(C1) <<  8) + \
    (static_cast<uint32_t>(C2) << 16) + \
    (static_cast<uint32_t>(C3) << 24) )

	//! \internal
	//!
	//! Converts a given symbol `s` of `sLen` to a keyword token.
	static uint32_t mpGetKeyword(const uint8_t* s, size_t sLen)
	{
		if (sLen == 3 && s[0] == 'v' && s[1] == 'a' && s[2] == 'r')
			return TokenType::kTokenVar;

		return TokenType::kTokenSymbol;
	}

	uint32_t Tokenizer::peek(Token* token)
	{
		uint32_t uToken = _token.token;
		if (uToken != TokenType::kTokenInvalid || (uToken = next(&_token)) != TokenType::kTokenInvalid)
			*token = _token;
		return uToken;
	}

	uint32_t Tokenizer::next(Token* token)
	{
		// Skip parsing if the next token is already done, caused by `peek()`.

		// the type of the currently parsed token
		uint32_t c = _token.token;
		uint32_t hVal;

		if (c != TokenType::kTokenInvalid)
		{
			*token = _token;
			_token.token = TokenType::kTokenInvalid;
			return c;
		}

		// Input string.
		const uint8_t* p = reinterpret_cast<const uint8_t*>(_p);
		const uint8_t* pToken = p;

		const uint8_t* pStart = reinterpret_cast<const uint8_t*>(_start);
		const uint8_t* pEnd = reinterpret_cast<const uint8_t*>(_end);

		// --------------------------------------------------------------------------
		// [Spaces]
		// --------------------------------------------------------------------------

	_Repeat:
		for (;;)
		{
			if (p == pEnd)
				goto _EndOfInput;

			hVal = p[0];
			c = mpCharClass[hVal];

			if (!isSpace(p[0]))
				break;
			p++;
		}

		// Save the first character of the token.
		pToken = p;

		// --------------------------------------------------------------------------
		// [Number | Dot]
		// --------------------------------------------------------------------------

		if (isNum(p[0]))
		{
			// Parsing floating point is not that simple as it looks. To simplify the
			// most common cases we parse floating points up to `kSafeDigits` and then
			// use libc `strtod()` function to parse numbers that are more complicated.
			//
			// http://www.exploringbinary.com/fast-path-decimal-to-floating-point-conversion/
			double val = 0.0;
			size_t digits = 0;

			// Skip leading zeros.
			while (p[0] == '0')
			{
				if (++p == pEnd)
					break;
			}

			// Parse significant.
			size_t scale = 0;
			while (p != pEnd)
			{
				c = static_cast<uint32_t>(p[0]) - static_cast<uint32_t>('0');
				if (c > 9)
					break;
				scale++;

				if (c != 0)
				{
					if (scale < kPow10TableSize)
						val = val * mpPow10Table[scale] + static_cast<double>(static_cast<int>(c));
					digits += scale;
					scale = 0;
				}

				p++;
			}
			size_t significantDigits = digits + scale;

			// Parse fraction.
			if (p != pEnd && p[0] == '.')
			{
				size_t scale = 0;

				while (++p != pEnd)
				{
					c = static_cast<uint32_t>(p[0]) - static_cast<uint32_t>('0');
					if (c > 9)
						break;
					scale++;

					if (c != 0)
					{
						if (scale < kPow10TableSize)
							val = val * mpPow10Table[scale] + static_cast<double>(static_cast<int>(c));
						digits += scale;
						scale = 0;
					}
				}

				// Token is '.'.
				if (size_t(p - pToken) == 1)
				{
					_p = reinterpret_cast<const char*>(p);
					return token->setData(size_t(pToken - pStart), size_t(p - pToken), 0, TokenType::kTokenDot);
				}
			}

			bool safe = digits <= kSafeDigits && significantDigits < 999999;
			int exponent = safe ? static_cast<int>(significantDigits) - static_cast<int>(digits) : 0;

			// Parse an optional exponent.
			if (p != pEnd && mpGetLower(p[0]) == 'e')
			{
				if (++p == pEnd)
					goto _Invalid;

				bool negative = p[0] == '-';
				if (negative || p[0] == '+')
					if (++p == pEnd)
						goto _Invalid;

				uint32_t e = 0;
				size_t eLen = 0;

				do
				{
					c = static_cast<uint32_t>(p[0]) - static_cast<uint32_t>('0');
					if (c > 9)
						break;

					e = e * 10 + c;
					eLen++;
				} while (++p != pEnd);

				// Error if there is no number after the 'e' token.
				if (eLen == 0)
					goto _Invalid;

				// If less than 10 digits it's safe to assume the exponent is zero if
				// `e` is zero. Otherwise it could have overflown the 32-bit integer.
				if (e == 0 && eLen < 10)
					eLen = 0; // No exponent.

				if (eLen <= 6)
					exponent += negative ? -static_cast<int>(e) : static_cast<int>(e);
				else
					safe = false;
			}

			// Error if there is an alpha-numeric character right next to the number that is not 'i'.
			if (p != pEnd && isSymbolFirst(p[0]) && p[0] != 'i')
				goto _Invalid;


			// Limit a range of safe values from Xe-15 to Xe15.
			safe = safe && exponent >= -kPow10TableSize && exponent <= kPow10TableSize;
			size_t len = (size_t)(p - pToken);

			// check whether there is a complex number or not and set the output accordingly.
			uint32_t tokenType;
			if (p[0] == 'i')
			{
				p++;
				len++;
				tokenType = TokenType::kTokenComplex;
			}
			else
			{
				tokenType = TokenType::kTokenNumber;
			}


			if (safe)
			{
				if (exponent != 0)
					val = exponent < 0 ? val / mpPow10Table[-exponent] : val * mpPow10Table[exponent];
			}
			else
			{
				// Using libc's strtod is not optimal, but it's precise for complex cases.
				char tmp[512];
				char* buf = tmp;

				if (len >= MATHPRESSO_ARRAY_SIZE(tmp) && (buf = static_cast<char*>(::malloc(len + 1))) == nullptr)
					return TokenType::kTokenInvalid;

				memcpy(buf, pToken, len);
				buf[len] = '\0';

				val = _strtod.conv(buf, nullptr);

				if (buf != tmp)
					::free(buf);
			}

			token->value = val;
			token->setData((size_t)(pToken - pStart), len, 0, tokenType);

			_p = reinterpret_cast<const char*>(p);
			return tokenType;
		}

		// --------------------------------------------------------------------------
		// [Symbol | Keyword]
		// --------------------------------------------------------------------------

		else if (isSymbolFirst(p[0]))
		{
			// We always generate the hVal during tokenization to improve performance.
			while (++p != pEnd)
			{
				uint32_t ord = p[0];
				if (!isSymbol(p[0]))
					break;
				hVal = HashUtils::hashChar(hVal, ord);
			}

			size_t len = (size_t)(p - pToken);
			_p = reinterpret_cast<const char*>(p);
			return token->setData((size_t)(pToken - pStart), len, hVal, mpGetKeyword(pToken, len));
		}

		// --------------------------------------------------------------------------
		// [Single-Char]
		// --------------------------------------------------------------------------

		else if (isSeparator(p[0]))
		{
			_p = reinterpret_cast<const char*>(++p);
			return token->setData(size_t(pToken - pStart), size_t(p - pToken), 0, c);
		}

		// --------------------------------------------------------------------------
		// [Single-Char | Multi-Char]
		// --------------------------------------------------------------------------

		else if (isOperator(p[0]))
		{
			size_t length(1);

			if (++p != pEnd && isOperator(p[0]))
			{
				if (++p != pEnd && isOperator(p[0]))
					length = 3;
				else
					length = 2;
			}

			if (std::string(reinterpret_cast<const char*>(pToken), length) == "//")
			{
				for (;;)
				{
					if (p == pEnd)
						goto _EndOfInput;

					if (p++[0] == '\n')
						goto _Repeat;
				}
			}

			_p = reinterpret_cast<const char*>(p);
			return token->setData(size_t(pToken - pStart), length, 0, TokenType::kTokenOperator);

		}

		// --------------------------------------------------------------------------
		// [Invalid]
		// --------------------------------------------------------------------------

	_Invalid:
		_p = reinterpret_cast<const char*>(pToken);
		return token->setData(size_t(pToken - pStart), size_t(p - pToken), 0, TokenType::kTokenInvalid);

		// --------------------------------------------------------------------------
		// [EOI]
		// --------------------------------------------------------------------------

	_EndOfInput:
		_p = _end;
		return token->setData(size_t(pToken - pStart), size_t(p - pToken), 0, TokenType::kTokenEnd);
	}

} // mathpresso namespace
