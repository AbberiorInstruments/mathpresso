// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MATHPRESSO_MPCOMPILER_P_H
#define _MATHPRESSO_MPCOMPILER_P_H

// [Dependencies]
#include <mathpresso/mpast_p.h>
#include <memory>

namespace mathpresso
{

	MATHPRESSO_NOAPI CompiledFunc mpCompileFunction(AstBuilder* ast, uint32_t options, OutputLog* log, const Symbols * ops, bool b_complex = false);
	MATHPRESSO_NOAPI void mpFreeFunction(void* fn);

	// ============================================================================
	// [mathpresso::JitVar]
	// ============================================================================

	struct JitVar
	{

		JitVar() : op(), isReadOnly(false)
		{
		}

		JitVar(asmjit::Operand op, bool isROnly) : op(op), isReadOnly(isROnly)
		{
		}

		JitVar(const JitVar& other) : op(other.op), isReadOnly(other.isRO())
		{
		}

		~JitVar()
		{
		}

		// Reset
		void reset()
		{
			op.reset();
			isReadOnly = false;
		}

		// Operator Overload.
		const JitVar& operator=(const JitVar& other)
		{
			op = other.op;
			this->isReadOnly = other.isRO();
			return *this;
		}

		// Operator Overload.
		const bool operator==(const JitVar& other)
		{
			return this->isRO() == other.isRO() && this->op.isEqual(other.op);
		}

		const bool operator!=(const JitVar& other)
		{
			return !(this->operator==(other));
		}

		// Swap.
		void swapWith(JitVar& other)
		{
			JitVar t(*this);
			*this = other;
			other = t;
		}

		// Operand management.
		const asmjit::Operand& getOperand() const
		{
			return op;
		}

		const asmjit::X86Mem& getMem() const
		{
			return *static_cast<const asmjit::X86Mem*>(&op);
		}

		const asmjit::X86Xmm& getXmm() const
		{
			return *static_cast<const asmjit::X86Xmm*>(&op);
		}

		bool isNone() const
		{
			return op.isNone();
		}

		bool isMem() const
		{
			return op.isMem();
		}

		bool isXmm() const
		{
			return op.isReg(asmjit::X86Reg::kRegXmm);
		}

		bool isRO() const
		{
			return isReadOnly;
		}

		void setRO()
		{
			isReadOnly = true;
		}

		void clearRO()
		{
			isReadOnly = false;
		}

		// Members.
	private:
		asmjit::Operand op;
		bool isReadOnly;
	};

	// ============================================================================
	// [mathpresso::JitCompiler]
	// ============================================================================

	struct MATHPRESSO_NOAPI JitCompiler
	{
		JitCompiler(asmjit::ZoneHeap* heap, asmjit::X86Compiler* cc, const Symbols * ops)
			: heap(heap),
			cc(cc),
			varSlots({}),
			functionBody(nullptr),
			constPool(&cc->_cbDataZone),
			_ops(ops)
		{
			enableSSE4_1 = asmjit::CpuInfo::getHost().hasFeature(asmjit::CpuInfo::kX86FeatureSSE4_1);
		}

		~JitCompiler() {}

		// Function Generator.
		void beginFunction();
		void endFunction();

		// Variable Management.
		JitVar copyVar(const JitVar& other, bool isRO);
		JitVar writableVar(const JitVar& other);
		JitVar registerVar(const JitVar& other);

		JitVar copyVarComplex(const JitVar & other, bool isRO);
		JitVar writableVarComplex(const JitVar & other);
		JitVar registerVarComplex(const JitVar & other, bool otherIsNonComplex = false);
		JitVar registerVarAsComplex(const JitVar & other);

		// Compiler.
		void compile(std::shared_ptr<AstBlock> node, AstScope* rootScope, uint32_t numSlots, bool b_complex);

		JitVar onNode(std::shared_ptr<AstNode> node);
		JitVar onBlock(std::shared_ptr<AstBlock> node);
		JitVar onVar(std::shared_ptr<AstVar> node);
		JitVar onImm(std::shared_ptr<AstImm> node);

		// Helpers.
		void inlineCallDRetD(const asmjit::X86Xmm& dst, const asmjit::X86Xmm* args, size_t count, void* fn);
		void inlineCallDRetC(const asmjit::X86Xmm & dst, const asmjit::X86Xmm * args, size_t count, void * fn);
		void inlineCallCRetD(const asmjit::X86Xmm & dst, const asmjit::X86Xmm * args, size_t count, void * fn);
		void inlineCallCRetC(const asmjit::X86Xmm & dst, const asmjit::X86Xmm * args, size_t count, void * fn);

		template<typename RET, typename PARAMS>
		void inlineCall(const asmjit::X86Xmm & dst, const asmjit::X86Xmm * args, size_t count, void * fn);

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

		std::vector<JitVar> varSlots;
		asmjit::CBNode* functionBody;

		asmjit::Label constLabel;
		asmjit::X86Gp constPtr;
		asmjit::ConstPool constPool;

		bool enableSSE4_1;

		const Symbols * _ops;

	};

} // mathpresso namespace

// [Guard]
#endif // _MATHPRESSO_MPCOMPILER_P_H
