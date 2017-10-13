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
	//! Creates a JitCompiler, initializes it and compiles the AST provided in ast.
	//! it then returns the resulting Function. If there was an error while compiling,
	//! execution of the returned function results in undefined behavior.
	CompiledFunc mpCompileFunction(std::shared_ptr<AstBuilder> ast, uint32_t options, OutputLog* log, std::shared_ptr<Context> ctx, bool b_complex = false);

	//! Destructor of a Compiled function.
	void mpFreeFunction(void* fn);

	// ============================================================================
	// [mathpresso::JitVar]
	// ============================================================================

	//! Represents a memory-location or a register, so it can be used as a parameter
	//! of an Operation or to save its return-value.
	struct JitVar
	{

		JitVar() noexcept
			: op(),
			isReadOnly(false)
		{
		}

		JitVar(asmjit::Operand op, bool isROnly) noexcept
			: op(op),
			isReadOnly(isROnly)
		{
		}

		JitVar(const JitVar& other) noexcept
			: op(other.op),
			isReadOnly(other.isRO())
		{
		}

		~JitVar() noexcept
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

		void setRO(bool isRO)
		{
			isReadOnly = isRO;
		}

		// Members.
	private:
		asmjit::Operand op;
		bool isReadOnly;
	};

	// ============================================================================
	// [mathpresso::JitCompiler]
	// ============================================================================

	//! This takes a AST and compiles it into machine-code, that can be executed by
	//! the processor.
	struct JitCompiler
	{
		JitCompiler(asmjit::X86Compiler* cc, std::shared_ptr<Context> ctx) noexcept
			: cc(cc),
			enableSSE4_1(asmjit::CpuInfo::getHost().hasFeature(asmjit::CpuInfo::kX86FeatureSSE4_1)),
			varSlots({}),
			functionBody(nullptr),
			constPool(&cc->_cbDataZone),
			_shadowContext(ctx)
		{
		}

		~JitCompiler() noexcept {}

		//TODO: check the following functions, whether they are all necessary
		// Variable Management.
		JitVar copyVar(const JitVar& other, bool isRO);
		JitVar writableVar(const JitVar& other);
		JitVar registerVar(const JitVar& other);

		JitVar copyVarComplex(const JitVar & other, bool isRO);
		JitVar writableVarComplex(const JitVar & other);
		JitVar registerVarComplex(const JitVar & other, bool otherIsNonComplex = false);
		JitVar registerVarAsComplex(const JitVar & other);

		// Compiler.
		void compile(std::shared_ptr<AstBlock> node, std::shared_ptr<Context> rootContext, uint32_t numSlots, bool b_complex);

		JitVar onNode(std::shared_ptr<AstNode> node);

		template<typename RET, typename ARGS>
		void inlineCall(const asmjit::X86Xmm & dst, const asmjit::X86Xmm * args, size_t count, void * fn);

		// Constants.
		JitVar getConstantU64(uint64_t value);
		JitVar getConstantU64(uint64_t real, uint64_t imag);
		JitVar getConstantU64AsPD(uint64_t value);
		JitVar getConstantD64(double value);
		JitVar getConstantD64(std::complex<double> value);
		JitVar getConstantD64AsPD(double value);

		// Members.
		asmjit::X86Compiler* cc;
		bool enableSSE4_1;
		std::vector<JitVar> varSlots;

	private:

		// Function Generator.
		void beginFunction();
		void endFunction();

		JitVar onBlock(std::shared_ptr<AstBlock> node);
		JitVar onVar(std::shared_ptr<AstVar> node);
		JitVar onImm(std::shared_ptr<AstImm> node);

		void prepareConstPool();

		asmjit::X86Gp resultAddress;
		asmjit::X86Gp variablesAddress;

		asmjit::CBNode* functionBody;

		asmjit::Label constLabel;
		asmjit::X86Gp constPtr;
		asmjit::ConstPool constPool;

		std::shared_ptr<Context> _shadowContext;

	};

} // mathpresso namespace

// [Guard]
#endif // _MATHPRESSO_MPCOMPILER_P_H
