// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MATHPRESSO_MPCOMPILER_P_H
#define _MATHPRESSO_MPCOMPILER_P_H

// [Dependencies]
#include "./mpast_p.h"

namespace mathpresso
{

	MATHPRESSO_NOAPI CompiledFunc mpCompileFunction(AstBuilder* ast, uint32_t options, OutputLog* log, bool b_complex = false);
	MATHPRESSO_NOAPI void mpFreeFunction(void* fn);

	// ============================================================================
	// [mathpresso::JitVar]
	// ============================================================================

	struct JitVar 
	{
		enum FLAGS 
		{
			FLAG_NONE = 0,
			FLAG_RO = 1
		};

		MATHPRESSO_INLINE JitVar() : op(), flags(FLAG_NONE) {}
		MATHPRESSO_INLINE JitVar(asmjit::Operand op, uint32_t flags) : op(op), flags(flags) {}
		MATHPRESSO_INLINE JitVar(const JitVar& other) : op(other.op), flags(other.flags) {}
		MATHPRESSO_INLINE ~JitVar() {}

		// Reset
		MATHPRESSO_INLINE void reset() 
		{
			op.reset();
			flags = FLAG_NONE;
		}

		// Operator Overload.
		MATHPRESSO_INLINE const JitVar& operator=(const JitVar& other)
		{
			op = other.op;
			flags = other.flags;
			return *this;
		}

		// Swap.
		MATHPRESSO_INLINE void swapWith(JitVar& other)
		{
			JitVar t(*this);
			*this = other;
			other = t;
		}

		// Operand management.
		MATHPRESSO_INLINE const asmjit::Operand& getOperand() const { return op; }
		MATHPRESSO_INLINE const asmjit::X86Mem& getMem() const { return *static_cast<const asmjit::X86Mem*>(&op); }
		MATHPRESSO_INLINE const asmjit::X86Xmm& getXmm() const { return *static_cast<const asmjit::X86Xmm*>(&op); }

		MATHPRESSO_INLINE bool isNone() const { return op.isNone(); }
		MATHPRESSO_INLINE bool isMem() const { return op.isMem(); }
		MATHPRESSO_INLINE bool isXmm() const { return op.isReg(asmjit::X86Reg::kRegXmm); }

		// Flags.
		MATHPRESSO_INLINE bool isRO() const { return (flags & FLAG_RO) != 0; }
		MATHPRESSO_INLINE void setRO() { flags |= FLAG_RO; }
		MATHPRESSO_INLINE void clearRO() { flags &= ~FLAG_RO; }

		// Members.
		asmjit::Operand op;
		uint32_t flags;
	};

	// ============================================================================
	// [mathpresso::JitCompiler]
	// ============================================================================

	struct MATHPRESSO_NOAPI JitCompiler 
	{
		JitCompiler(asmjit::ZoneHeap* heap, asmjit::X86Compiler* cc)
			: heap(heap),
			cc(cc),
			varSlots(nullptr),
			functionBody(nullptr),
			constPool(&cc->_cbDataZone) 
		{

			enableSSE4_1 = asmjit::CpuInfo::getHost().hasFeature(asmjit::CpuInfo::kX86FeatureSSE4_1);
		}

		~JitCompiler() {}

		// Function Generator.
		void beginFunction();
		void endFunction();

		// Variable Management.
		JitVar copyVar(const JitVar& other, uint32_t flags);
		JitVar writableVar(const JitVar& other);
		JitVar registerVar(const JitVar& other);

		JitVar copyVarComplex(const JitVar & other, uint32_t flags);
		JitVar writableVarComplex(const JitVar & other);
		JitVar registerVarComplex(const JitVar & other, bool otherIsComplex = false);
		JitVar registerVarAsComplex(const JitVar & other);

		// Compiler.
		void compile(AstBlock* node, AstScope* rootScope, uint32_t numSlots, bool b_complex);

		JitVar onNode(AstNode* node);
		JitVar onBlock(AstBlock* node);
		JitVar onVarDecl(AstVarDecl* node);
		JitVar onVar(AstVar* node);
		JitVar onImm(AstImm* node);
		JitVar onUnaryOp(AstUnaryOp* node);
		JitVar onBinaryOp(AstBinaryOp* node);
		JitVar onTernaryOp(AstTernaryOp * node);
		JitVar onCall(AstCall* node);

		// Helpers.
		void inlineRound(const asmjit::X86Xmm& dst, const asmjit::X86Xmm& src, uint32_t op, bool takesComplex, bool returnsComplex);
		void inlineCallAbstract(const asmjit::X86Xmm& dst, const asmjit::X86Xmm* args, uint32_t count, uint32_t op, bool takesComplex, bool returnsComplex);
		void inlineCallAbstract(const asmjit::X86Xmm & dst, const asmjit::X86Xmm * args, uint32_t count, void * fp, bool takesComplex, bool returnsComplex);
		void inlineCall(const asmjit::X86Xmm& dst, const asmjit::X86Xmm* args, uint32_t count, void* fn);
		void inlineCallDRetC(const asmjit::X86Xmm & dst, const asmjit::X86Xmm * args, uint32_t count, void * fn);
		void inlineCallCRetD(const asmjit::X86Xmm & dst, const asmjit::X86Xmm * args, const uint32_t count, void * fn);
		void inlineCallComplex(const asmjit::X86Xmm & dst, const asmjit::X86Xmm * args, uint32_t count, void * fn);

		// Constants.
		void prepareConstPool();
		JitVar getConstantU64(uint64_t value);
		JitVar getConstantU64(uint64_t real, uint64_t imag);
		JitVar getConstantU64AsPD(uint64_t value);
		JitVar getConstantD64(double value);
		JitVar getConstantD64(std::complex<double> value);
		JitVar getConstantD64AsPD(double value);

		// Members.
		asmjit::X86Compiler* cc;
		asmjit::ZoneHeap* heap;

		asmjit::X86Gp resultAddress;
		asmjit::X86Gp variablesAddress;

		JitVar* varSlots;
		asmjit::CBNode* functionBody;

		asmjit::Label constLabel;
		asmjit::X86Gp constPtr;
		asmjit::ConstPool constPool;

		bool enableSSE4_1;
	};

} // mathpresso namespace

// [Guard]
#endif // _MATHPRESSO_MPCOMPILER_P_H
