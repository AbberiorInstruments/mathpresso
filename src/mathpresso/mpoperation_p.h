#pragma once
#include <mathpresso/mpoperation.h>

namespace mathpresso
{
	// isfinite
	template<typename T>
	class MpOperationIsFinite : public MpOperationEval<T, T>
	{
	public:
		MpOperationIsFinite() : MpOperationEval<T, T>(1)
		{
		}
		virtual T evaluate(const T *) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	// isinf
	template<typename T>
	class MpOperationIsInfinite : public MpOperationFunc<T, T>
	{
	public:
		MpOperationIsInfinite() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	// isnan
	template<typename T>
	class MpOperationIsNan : public MpOperationFunc<T, T>
	{
	public:
		MpOperationIsNan() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	// real
	class MpOperationGetReal : public MpOperationFunc<double, std::complex<double>>
	{
	public:
		MpOperationGetReal() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	// imag
	class MpOperationGetImag : public MpOperationFunc<double, std::complex<double>>
	{
	public:
		MpOperationGetImag() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	// Square root
	class MpOperationSqrt : public MpOperationFunc<double, double>
	{
	public:
		MpOperationSqrt() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	// Negation
	template<typename T>
	class MpOperationNeg : public MpOperationFunc<T, T>
	{
	public:
		MpOperationNeg() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	protected:
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const override;
	};

	// Not
	template<typename T>
	class MpOperationNot : public MpOperationFunc<T, T>
	{
	public:
		MpOperationNot() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	// conjugate
	class MpOperationConjug :public MpOperationFunc<std::complex<double>, std::complex<double>>
	{
	public:
		MpOperationConjug() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	private:
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const override;
	};

	// Reciprocate
	template<typename T>
	class MpOperationRecip :public MpOperationFunc<T,T>
	{
	public:
		MpOperationRecip() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	// Sign bit
	class MpOperationSignBit : public MpOperationFunc<double, double>
	{
	public:
		MpOperationSignBit() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	// copy sign
	class MpOperationCopySign : public MpOperationFunc<double, double>
	{
	public:
		MpOperationCopySign() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	// Average
	template<typename T>
	class MpOperationAvg : public MpOperationFunc<T,T>
	{
	public:
		MpOperationAvg() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	// Absolute
	class MpOperationAbs : public MpOperationFunc<double, double>
	{
	public:
		MpOperationAbs() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	// round
	class MpOperationRound : public MpOperationFunc<double, double>
	{
	public:
		MpOperationRound() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	protected:
		std::string description_ = "Rounds *.5 towards infinity.";
	};

	// roundeven
	class MpOperationRoundEven : public MpOperationFunc<double, double>
	{
	public:
		MpOperationRoundEven() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	protected:
		std::string description_ = "Rounds *.5 towards nearest even integer.";
	};

	// trunc
	class MpOperationTrunc : public MpOperationFunc<double, double>
	{
	public:
		MpOperationTrunc() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	// frac
	class MpOperationFrac : public MpOperationFunc<double, double>
	{
	public:
		MpOperationFrac() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	// floor
	class MpOperationFloor : public MpOperationFunc<double, double>
	{
	public:
		MpOperationFloor() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	// ceil
	class MpOperationcCeil : public MpOperationFunc<double, double>
	{
	public:
		MpOperationcCeil() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};


	// Addition
	template<typename T>
	class MpOperationAdd : public MpOperationBinary<T>
	{
	public:
		MpOperationAdd() : MpOperationBinary(MpOperationBinary::NopIfZero | MpOperationBinary::IsCommutativ, 6)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T evaluate(T vL, T vR) const override
		{
			return vL + vR;
		}
	};

	// Subtraction
	template<typename T>
	class MpOperationSub : public MpOperationBinary<T>
	{
	public:
		MpOperationSub() : MpOperationBinary(MpOperationBinary::NopIfZero, 6)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T evaluate(T vL, T vR) const override
		{
			return vL - vR;
		}
	};

	// Multiplication
	template<typename T>
	class MpOperationMul : public MpOperationBinary<T>
	{
	public:
		MpOperationMul() : MpOperationBinary(MpOperationBinary::NopIfZero | MpOperationBinary::IsCommutativ, 5)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T evaluate(T vl, T vr) const override
		{
			return vl*vr;
		}
	};

	// Division
	template<typename T>
	class MpOperationDiv : public MpOperationBinary<T>
	{
	public:
		MpOperationDiv() : MpOperationBinary(MpOperationBinary::NopIfLOne, 5)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T evaluate(T vl, T vr) const override
		{
			return vl / vr;
		}
	};

	// Minimum
	class MpOperationMin : public MpOperationBinary<double>
	{
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double evaluate(double vl, double vr) const override
		{
			return std::min(vl, vr);
		}
	};

	// Maximum
	class MpOperationMax : public  MpOperationBinary<double>
	{
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double evaluate(double vl, double vr) const override
		{
			return std::max(vl, vr);
		}
	};

	// Equality
	template<typename T>
	class MpOperationEq : public MpOperationBinary<T>
	{
	public:
		MpOperationEq() : MpOperationBinary(MpOperationBinary::IsCommutativ, 9)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T evaluate(T vl, T vr) const override
		{
			return vl == vr ? 1.0 : 0.0;
		}
	};

	// Inequality
	template<typename T>
	class MpOperationNe : public MpOperationBinary<T>
	{
	public:
		MpOperationNe() : MpOperationBinary(MpOperationBinary::IsCommutativ, 9)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T evaluate(T vl, T vr) const override
		{
			return vl != vr ? 1.0 : 0.0;
		}
	};

	// Lesser than
	class MpOperationLt : public  MpOperationBinary<double>
	{
	public:
		MpOperationLt()	: MpOperationBinary(MpOperation::None, 8)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double evaluate(double vl, double vr) const override
		{
			return vl < vr ? 1. : 0.;
		}
	};

	// Lesser Equal
	class MpOperationLe : public  MpOperationBinary<double>
	{
	public:
		MpOperationLe() : MpOperationBinary(MpOperation::None, 8)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double evaluate(double vl, double vr) const override
		{
			return vl <= vr ? 1. : 0.;
		}
	};

	// Greater than
	class MpOperationGt : public  MpOperationBinary<double>
	{
	public:
		MpOperationGt() : MpOperationBinary(MpOperation::None, 8)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double evaluate(double vl, double vr) const override
		{
			return vl > vr ? 1. : 0.;
		}
	};

	// Greater equal
	class MpOperationGe : public  MpOperationBinary<double>
	{
	public:
		MpOperationGe() : MpOperationBinary(MpOperation::None, 8)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double evaluate(double vl, double vr) const override
		{
			return vl >= vr ? 1. : 0.;
		}
	};

	// Modulo
	class MpOperationModulo : public  MpOperationBinary<double>
	{
	public:
		MpOperationModulo() : MpOperationBinary<double>(MpOperation::None, 5)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double evaluate(double vl, double vr) const override
		{
			return std::fmod(vl, vr);
		}
	};

	// Ternary Operation
	template<typename T>
	class MpOperationTernary : public MpOperation
	{
	public:
		MpOperationTernary() noexcept;

		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const override;
	protected:
		bool isColon_;
	};

	// Variable declaration
	template<typename T>
	class MpOperationVarDeclaration : public MpOperation
	{
	public:
		MpOperationVarDeclaration() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const override;
	};

	// Assignment
	template<typename T>
	class MpOperationAssignment : public MpOperation
	{
	public:
		MpOperationAssignment() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const override;
	};
}