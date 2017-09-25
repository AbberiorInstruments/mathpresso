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
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
		virtual T evaluate(const T *) const override;
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
		MpOperationAdd() noexcept;

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T calculate(T vl, T vr) const override;
	};

	// Subtraction
	template<typename T>
	class MpOperationSub : public MpOperationBinary<T>
	{
	public:
		MpOperationSub() noexcept;

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T calculate(T vl, T vr) const override;
	};

	// Multiplication
	template<typename T>
	class MpOperationMul : public MpOperationBinary<T>
	{
	public:
		MpOperationMul() noexcept;

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T calculate(T vl, T vr) const override;
	};

	// Division
	template<typename T>
	class MpOperationDiv : public MpOperationBinary<T>
	{
	public:
		MpOperationDiv() noexcept;

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T calculate(T vl, T vr) const override;
	};

	// Minimum
	class MpOperationMin : public MpOperationBinary<double>
	{
	public:
		MpOperationMin() noexcept
			: MpOperationBinary<double>(Signature(2), MpOperation::None, 0)
		{
		}

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculate(double vl, double vr) const override;
	};

	// Maximum
	class MpOperationMax : public  MpOperationBinary<double>
	{
	public:
		MpOperationMax() noexcept
			: MpOperationBinary<double>(Signature(2), MpOperation::None, 0)
		{
		}

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculate(double vl, double vr) const override;
	};

	// Equality
	template<typename T>
	class MpOperationEq : public MpOperationBinary<T>
	{
	public:
		MpOperationEq() noexcept;

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T calculate(T vl, T vr) const override;
	};

	// Inequality
	template<typename T>
	class MpOperationNe : public MpOperationBinary<T>
	{
	public:
		MpOperationNe() noexcept;

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T calculate(T vl, T vr) const override;
	};

	// Lesser than
	class MpOperationLt : public  MpOperationBinary<double>
	{
	public:
		MpOperationLt() noexcept
			: MpOperationBinary<double>(Signature(2, Signature::type::real), MpOperation::None, 8)
		{
		}

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculate(double vl, double vr) const override;
	};

	// Lesser Equal
	class MpOperationLe : public  MpOperationBinary<double>
	{
	public:
		MpOperationLe() noexcept
			: MpOperationBinary<double>(Signature(2, Signature::type::real), MpOperation::None, 8)
		{
		}

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculate(double vl, double vr) const override;
	};

	// Greater than
	class MpOperationGt : public  MpOperationBinary<double>
	{
	public:
		MpOperationGt() noexcept
			: MpOperationBinary<double>(Signature(2, Signature::type::real), MpOperation::None, 8)
		{
		}

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculate(double vl, double vr) const override;
	};

	// Greater equal
	class MpOperationGe : public  MpOperationBinary<double>
	{
	public:
		MpOperationGe()  noexcept
			: MpOperationBinary<double>(Signature(2, Signature::type::real), MpOperation::None, 8)
		{
		}

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculate(double vl, double vr) const override;
	};

	// Modulo
	class MpOperationModulo : public  MpOperationBinary<double>
	{
	public:
		MpOperationModulo() noexcept
			: MpOperationBinary<double>(Signature(2, Signature::type::real), MpOperation::None, 5)
		{
		}

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculate(double vl, double vr) const override;
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