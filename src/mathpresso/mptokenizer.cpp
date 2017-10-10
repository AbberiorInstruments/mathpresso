// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Export]
#define MATHPRESSO_EXPORTS

// [De_endencies]
#include <mathpresso/mpeval_p.h>
#include <mathpresso/mptokenizer_p.h>
#include <iostream>

namespace mathpresso
{

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
		return  (c == ' ' || c == '\t' || c == '\n');
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

	namespace InternalConsts
	{
		constexpr int kSafeDigits = 15;
	}

#define CHAR4X(C0, C1, C2, C3) \
  ( (static_cast<uint32_t>(C0)      ) + \
    (static_cast<uint32_t>(C1) <<  8) + \
    (static_cast<uint32_t>(C2) << 16) + \
    (static_cast<uint32_t>(C3) << 24) )

	//! \internal
	//!
	//! Converts a given symbol `s` of `sLen` to a keyword token.
	static uint32_t mpGetKeyword(const char* s, size_t sLen)
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
		const char* p = _p;
		const char* pToken = p;

		// --------------------------------------------------------------------------
		// [Spaces]
		// --------------------------------------------------------------------------

		while (true)
		{
			for (;;)
			{
				if (p == _end)
					return endOfInput(token, pToken, p);

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
					if (++p == _end)
						break;
				}

				// Parse significant.
				size_t scale = 0;
				while (p != _end)
				{
					c = static_cast<uint32_t>(p[0] - '0');
					if (c > 9)
						break;
					scale++;

					if (c != 0)
					{
						if (scale < InternalConsts::kSafeDigits)
							val = val * std::pow(10, scale) + static_cast<double>(c);
						digits += scale;
						scale = 0;
					}

					p++;
				}
				size_t significantDigits = digits + scale;

				// Parse fraction.
				if (p != _end && p[0] == '.')
				{
					size_t scale_fraction = 0;

					while (++p != _end)
					{
						c = static_cast<uint32_t>(p[0] - '0');
						if (c > 9)
							break;
						scale_fraction++;

						if (c != 0)
						{
							if (scale_fraction < InternalConsts::kSafeDigits)
								val = val * std::pow(10, scale_fraction) + static_cast<double>(c);
							digits += scale_fraction;
							scale_fraction = 0;
						}
					}

					// Token is '.'.
					if (p - pToken == 1)
					{
						_p = p;
						return token->setData(pToken - _start, p - pToken, TokenType::kTokenDot);
					}
				}

				bool safe = digits <= InternalConsts::kSafeDigits && significantDigits < 999999;
				int exponent = safe ? static_cast<int>(significantDigits) - static_cast<int>(digits) : 0;

				// Parse an optional exponent.
				if (p != _end && mpGetLower(p[0]) == 'e')
				{
					if (++p == _end)
						return invalid(token, pToken, p);

					bool negative = p[0] == '-';
					if (negative || p[0] == '+')
						if (++p == _end)
							return invalid(token, pToken, p);

					uint32_t e = 0;
					size_t eLen = 0;

					do
					{
						c = static_cast<uint32_t>(p[0] - '0');
						if (c > 9)
							break;

						e = e * 10 + c;
						eLen++;
					} while (++p != _end);

					// Error if there is no number after the 'e' token.
					if (eLen == 0)
						return invalid(token, pToken, p);

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
				if (p != _end && isSymbolFirst(p[0]) && p[0] != 'i')
					return invalid(token, pToken, p);


				// Limit a range of safe values from Xe-15 to Xe15.
				safe = safe && exponent >= -InternalConsts::kSafeDigits && exponent <= InternalConsts::kSafeDigits;
				size_t len = p - pToken;

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
						val = exponent < 0 ? val / std::pow(10, -exponent) : val * std::pow(10, exponent);
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
				token->setData(pToken - _start, len, tokenType);

				_p = p;
				return tokenType;
			}

			// --------------------------------------------------------------------------
			// [Symbol | Keyword]
			// --------------------------------------------------------------------------

			else if (isSymbolFirst(p[0]))
			{
				// We always generate the hVal during tokenization to improve performance.
				while (++p != _end)
				{
					if (!isSymbol(p[0]))
						break;
				}

				size_t len = p - pToken;
				_p = p;
				return token->setData(pToken - _start, len, mpGetKeyword(pToken, len));
			}

			// --------------------------------------------------------------------------
			// [Single-Char]
			// --------------------------------------------------------------------------

			else if (isSeparator(p[0]))
			{
				_p = ++p;
				return token->setData(pToken - _start, p - pToken, c);
			}

			// --------------------------------------------------------------------------
			// [Single-Char | Multi-Char]
			// --------------------------------------------------------------------------

			else if (isOperator(p[0]))
			{
				size_t length(1);

				if (++p != _end && isOperator(p[0]))
				{
					if (++p != _end && isOperator(p[0]))
						length = 3;
					else
						length = 2;
				}

				if (std::string(pToken, length) == "//")
				{
					for (;;)
					{
						if (p == _end)
							return endOfInput(token, pToken, p);

						if (p++[0] == '\n')
							continue;
					}
				}

				_p = p;
				return token->setData(pToken - _start, length, TokenType::kTokenOperator);

			}
		}

		return invalid(token, pToken, p);
	}

	uint32_t Tokenizer::invalid(mathpresso::Token * token, const char * pToken, const char * p)
	{
		_p = pToken;
		return token->setData(pToken - _start, p - pToken, TokenType::kTokenInvalid);
	}

	uint32_t Tokenizer::endOfInput(mathpresso::Token * token, const char * pToken, const char * p)
	{
		_p = _end;
		return token->setData(pToken - _start, p - pToken, TokenType::kTokenEnd);
	}

} // mathpresso namespace
