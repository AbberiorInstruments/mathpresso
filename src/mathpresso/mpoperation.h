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

namespace mathpresso
{

	// Forward Declarations:
	struct JitCompiler;
	struct AstNode;
	struct AstOptimizer;
	struct JitVar;
	struct Context;

	typedef JitVar(*mpAsmFunc)(JitCompiler*, JitVar*);

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

		template<>
		struct TypeId<std::complex<double>>
		{
			static const type id_ = type::complex;
		};

		template<>
		struct TypeId<double>
		{
			static const type id_ = type::real;
		};

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
	private :
		std::string typeToString(type _type);
		void init(type retType, std::vector<param> params);
	};

	template<typename S>
	class TypedSignature;

	template<typename RET, typename ...ARGS>
	class TypedSignature<RET(ARGS...)>
	{
		TypedSignature() : Signature(sizeof(ARGS), std::common_type_t<ARGS...>)
		{
		}
	};

	class MATHPRESSO_API MpOperation : public MpObject
	{
	public:
		enum Flags
		{
			None = 0,
			FlagHasState = 0x00000002,
			RighttoLeft = 0x00000008,
			IsAssignment = 0x0001000,
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
			return (flags_ & RighttoLeft) != 0;
		}
		uint32_t precedence() const
		{
			return precedence_;
		}
		uint32_t flags() const
		{
			return flags_;
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

	template<typename RET, typename ARGS>
	class MATHPRESSO_API MpOperationFunc : public MpOperationEval<RET, ARGS>
	{
	public:
		MpOperationFunc(void * fnPtr, size_t numargs, uint32_t flags = MpOperation::None, uint32_t priority = 0) noexcept :
			MpOperationEval(numargs, flags, priority),
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
			: MpOperationEval(2, flags, priority)
		{
		}

		// calls generatAsm() after setting up.
		virtual JitVar compile(JitCompiler* jc, std::shared_ptr<AstNode> node) const override;
		// uses calculate() to calculate immediate values.
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const override;
		virtual T evaluate(T vL, T vR) const = 0;

		virtual T evaluate(const T * v) const override
		{
			throw std::exception("Coding error!");
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
	using obj_ptr = std::shared_ptr<mathpresso::MpOperation>;

	template<typename Signature, typename R, typename A, size_t N>
	struct CallerBase
	{
		static const size_t NUM_ARGS = N;
		using ARG_TYPE = A;
		using RET_TYPE = R;
		using SIGNATURE = Signature;
	};

	// Create a function with argument array pointer that calls a multi-parameter function
	template<typename Signature, typename R, typename A, void * FPTR, size_t N>
	struct Caller;

	template<typename Signature, typename R, typename A, void * FPTR>
	struct Caller<Signature, R, A, FPTR, 1> : CallerBase<Signature, R, A, 1>
	{
		static R call(A * args)
		{
			return SIGNATURE(FPTR)(args[0]);
		}
	};

	template<typename Signature, typename R, typename A, void * FPTR>
	struct Caller<Signature, R, A, FPTR, 2> : CallerBase<Signature, R, A, 2>
	{
		static R call(A * args)
		{
			return SIGNATURE(FPTR)(args[0], args[1]);
		}
	};

	template<typename Signature, typename R, typename A, void * FPTR>
	struct Caller<Signature, R, A, FPTR, 3> : CallerBase<Signature, R, A, 3>
	{
		static R call(A * args)
		{
			return SIGNATURE(FPTR)(args[0], args[1], args[2]);
		}
	};

	template<typename FPTR_T, void * FPTR>
	struct Caller_;

	template<void * FPTR, typename R, typename ...ARGS>
	struct Caller_<R(*)(ARGS...), FPTR> : Caller<R(*)(ARGS...), R, std::common_type_t<ARGS...>, (void *)FPTR, sizeof...(ARGS)>
	{
	};

	template<typename CALLER>
	obj_ptr _mpObject(const CALLER &c, uint32_t flags = mathpresso::MpOperation::None, uint32_t priority = 0)
	{
		return std::make_shared<mathpresso::MpOperationFunc<typename CALLER::RET_TYPE, typename CALLER::ARG_TYPE>>((void *)CALLER::call, CALLER::NUM_ARGS, flags, priority);
	}
}

using cplx_t = std::complex<double>;

#define _OBJ(expr) fobj::_mpObject(fobj::Caller_<decltype(expr), expr>())
#define VPTR(function) reinterpret_cast<void*>(function)

#endif //_MP_OPERATION_P_H