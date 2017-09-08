// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MATHPRESSO_MPHASH_P_H
#define _MATHPRESSO_MPHASH_P_H

#include  <mathpresso/mathpresso_p.h>

namespace mathpresso
{

	// ============================================================================
	// [mathpresso::HashUtils]
	// ============================================================================

	namespace HashUtils
	{
		// \internal
		static uint32_t hashChar(uint32_t hash, uint32_t c)
		{
			return hash * 65599 + c;
		}
	};

} // mathpresso namespace

// [Guard]
#endif // _MATHPRESSO_MPHASH_P_H
