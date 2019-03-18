#pragma once
// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MP_OPERATION_P_H
#define _MP_OPERATION_P_H

#include <mathpresso/mathpresso.h>
#include <complex>
#include <limits>
#include <type_traits>

namespace mathpresso
{

	// Forward Declarations:
	struct JitCompiler;
	struct AstNode;
	struct AstOptimizer;
	struct JitVar;
	struct Context;

	typedef JitVar(*mpAsmFunc)(JitCompiler*, JitVar*);

	//! Description of the signature of a Operation: the number and type of parameters
	//! and type of the return. The type can be real (double) or a complex value (std::complex<double>).
	struct MATHPRESSO_API Signature
	{
		// Types allowed
		enum class type
		{
			real = 0,
			complex = 1
		};

		// Currently all params need to be the same type but we already cater for the next step...
		struct param
		{
			type type_;
			std::string name_;
		};

		template<typename T>
		struct TypeId;

		Signature() noexcept
		{
		}

		Signature(size_t nargs) noexcept;
		Signature(size_t nargs, type cmnType) noexcept;
		Signature(size_t nargs, type argType, type retType) noexcept;

		bool areParams(type _type) const;
		std::string to_string();

		type return_type_;
		std::vector<param> parameters_;
	private:
		std::string typeToString(type _type);
		void init(type retType, std::vector<param> params);
	};

	template<>
	struct Signature::TypeId<std::complex<double>>
	{
		static const type id_ = type::complex;
	};

	template<>
	struct Signature::TypeId<double>
	{
		static const type id_ = type::real;
	};

	//! An interface for functions or operations with a signature.
	//! Implementations of this hold information of how to optimize 
	//! a instruction in the AST and how to compile it to machine-code.
	class MATHPRESSO_API MpOperation : public MpObject
	{
	public:
		enum Flags
		{
			None = 0,
			FlagHasState = 0x00000002,
			RighttoLeft = 0x00000008,
			IsAssignment = 0x00001000,
			IsMethod = 0x10000000
		};

		MpOperation(const Signature &s, uint32_t flags, uint32_t precedence = 0)  noexcept
			: signature_(s),
			flags_(flags),
			precedence_(precedence)
		{
		}

		virtual ~MpOperation() noexcept
		{
		}

		// Add ASM code to compiler stack 
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode>  node) const = 0;
		// Optimize AST 
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const = 0;

		size_t nargs() const
		{
			return signature_.parameters_.size();
		}
		bool isRightToLeft() const
		{
			return hasFlag(RighttoLeft);
		}
		uint32_t precedence() const
		{
			return precedence_;
		}
		Signature signature() const
		{
			return signature_;
		}
		bool hasFlag(uint32_t flag) const
		{
			return flag & flags_;
		}
	protected:
		Signature signature_;
		uint32_t flags_;
		uint32_t precedence_;
	};

	template<typename RET, typename ARGS>
	class MATHPRESSO_API MpOperationEval : public MpOperation
	{
	public:
		MpOperationEval(size_t numargs, uint32_t flags = MpOperation::None, uint32_t priority = 0) noexcept;
		virtual ~MpOperationEval() noexcept
		{
		}
		virtual RET evaluate(const ARGS * args) const = 0;
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const override;
	};

	//! The implementation of a function, that can be called by the expression.
	template<typename RET, typename ARGS>
	class MATHPRESSO_API MpOperationFunc : public MpOperationEval<RET, ARGS>
	{
	public:
		MpOperationFunc(void * fnPtr, size_t numargs, uint32_t flags = MpOperation::None, uint32_t priority = 0) noexcept :
			MpOperationEval<RET, ARGS>(numargs, flags, priority),
			fnPtr_(fnPtr)
		{
		}

		virtual ~MpOperationFunc() noexcept
		{
		}

		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
		virtual RET evaluate(const ARGS * args) const;
	protected:
		// Function-pointer:
		void * fnPtr_;
	};

	//! The implementation of a binary operation like '+' or '-'.
	template<typename T>
	class MATHPRESSO_API MpOperationBinary : public MpOperationEval<T, T>
	{
	public:
		enum Flags
		{
			NopIfLZero = 0x10000000,
			NopIfRZero = 0x20000000,
			NopIfLOne = 0x40000000,
			NopIfROne = 0x80000000,
			NopIfZero = NopIfLZero | NopIfRZero,
			NopIfOne = NopIfLOne | NopIfROne,
			IsCommutativ = 0x00000010
		};

		MpOperationBinary(uint32_t flags = MpOperation::None, uint32_t priority = 0) noexcept
			: MpOperationEval<T, T>(2, flags, priority)
		{
		}

		// calls generatAsm() after setting up.
		virtual JitVar compile(JitCompiler* jc, std::shared_ptr<AstNode> node) const override;
		// uses calculate() to calculate immediate values.
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const override;
		virtual T evaluate(T vL, T vR) const = 0;

		virtual T evaluate(const T * v) const override
		{
			throw std::runtime_error("Coding error!");
		}
	protected:
		// These are called by compile() and should only contain the asm-statements. vl will always
		// be in a register, vr can be in Register or in Memory.
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const = 0;
	};
}

#include <complex>

namespace fobj
{
	template<typename R, typename A, size_t N, uint32_t F = 0>
	struct Caller__
	{
		using ret_t = R;
		using arg_t = A;
		using nargs_t = std::integral_constant<int, N>;
		using flags_t = std::integral_constant<int, F>;
	};

	template<typename R, typename A, size_t N, typename FPTR_T, FPTR_T FPTR>
	class Caller_;

	template<typename R, typename A, typename FPTR_T, FPTR_T FPTR>
	class Caller_<R, A, 0, FPTR_T, FPTR> : public Caller__<R, A, 0>
	{
	public:
		static R call(const A * args)
		{
			return FPTR();
		}
	};

	template<typename R, typename A, typename FPTR_T, FPTR_T FPTR>
	class Caller_<R, A, 1, FPTR_T, FPTR> : public Caller__<R, A, 1>
	{
	public:
		static R call(const A * args)
		{
			return FPTR(args[0]);
		}
	};

	template<typename R, typename A, typename FPTR_T, FPTR_T FPTR>
	class Caller_<R, A, 2, FPTR_T, FPTR> : public Caller__<R, A, 2>
	{
	public:
		static R call(const A * args)
		{
			return FPTR(args[0], args[1]);
		}
	};

	template<typename R, typename A, typename FPTR_T, FPTR_T FPTR>
	class Caller_<R, A, 3, FPTR_T, FPTR> : public Caller__<R, A, 3>
	{
	public:
		static R call(const A * args)
		{
			return FPTR(args[0], args[1], args[2]);
		}
	};

	template<typename FPTR_T, FPTR_T FPTR>
	class Caller;

	template<typename R, typename ...ARGS, R(*FPTR)(ARGS...)>
	class Caller<R(*)(ARGS...), FPTR> : public Caller_<R, std::remove_reference_t<std::common_type_t<ARGS...>>, sizeof...(ARGS), R(*)(ARGS...), FPTR>
	{
	};

	template<typename R, R(*FPTR)()>
	class Caller<R(*)(), FPTR> : public Caller_<R, double, 0, R(*)(), FPTR>
	{
	};

	template<typename CALLER>
	std::shared_ptr<mathpresso::MpOperation> _mpObject(const CALLER &c, uint32_t flags = mathpresso::MpOperation::None, uint32_t priority = 0)
	{
		return 
			std::make_shared<mathpresso::MpOperationFunc<typename CALLER::ret_t, typename CALLER::arg_t>>(
				reinterpret_cast<void*>(CALLER::call),
				CALLER::nargs_t::value,
				flags | CALLER::flags_t::value,
				priority);
	}
}

using cplx_t = std::complex<double>;

#define _OBJ(expr) fobj::_mpObject(fobj::Caller<decltype(expr), expr>())
#define VPTR(function) reinterpret_cast<void*>(function)

#endif //_MP_OPERATION_P_H
