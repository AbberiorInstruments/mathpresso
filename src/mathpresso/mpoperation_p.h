#pragma once
#include <mathpresso/mpoperation.h>

namespace mathpresso
{
	// isfinite
	template<typename T>
	class MpOperationIsFinite : public MpOperationFunc<T, T>
	{
	public:
		MpOperationIsFinite();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// isinf
	template<typename T>
	class MpOperationIsInfinite : public MpOperationFunc<T, T>
	{
	public:
		MpOperationIsInfinite();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// isnan
	template<typename T>
	class MpOperationIsNan : public MpOperationFunc<T, T>
	{
	public:
		MpOperationIsNan();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// real
	class MpOperationGetReal : public MpOperationFunc<double, std::complex<double>>
	{
	public:
		MpOperationGetReal();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// imag
	class MpOperationGetImag : public MpOperationFunc<double, std::complex<double>>
	{
	public:
		MpOperationGetImag();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// Square root
	class MpOperationSqrt : public MpOperationFunc<double, double>
	{
	public:
		MpOperationSqrt();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// Negation
	template<typename T>
	class MpOperationNeg : public MpOperationFunc<T, T>
	{
	public:
		MpOperationNeg();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	protected:
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) const override;
	};

	// Not
	template<typename T>
	class MpOperationNot : public MpOperationFunc<T, T>
	{
	public:
		MpOperationNot();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// conjugate
	class MpOperationConjug :public MpOperationFunc<std::complex<double>, std::complex<double>>
	{
	public:
		MpOperationConjug();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	private:
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) const override;
	};

	// Reciprocate
	template<typename T>
	class MpOperationRecip :public MpOperationFunc<T,T>
	{
	public:
		MpOperationRecip();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// Sign bit
	class MpOperationSignBit : public MpOperationFunc<double, double>
	{
	public:
		MpOperationSignBit();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// copy sign
	class MpOperationCopySign : public MpOperationFunc<double, double>
	{
	public:
		MpOperationCopySign();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// Average
	template<typename T>
	class MpOperationAvg : public MpOperationFunc<T,T>
	{
	public:
		MpOperationAvg();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// Absolute
	class MpOperationAbs : public MpOperationFunc<double, double>
	{
	public:
		MpOperationAbs();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// round
	class MpOperationRound : public MpOperationFunc<double, double>
	{
	public:
		MpOperationRound();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	protected:
		std::string description_ = "Rounds *.5 towards infinity.";
	};

	// roundeven
	class MpOperationRoundEven : public MpOperationFunc<double, double>
	{
	public:
		MpOperationRoundEven();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	protected:
		std::string description_ = "Rounds *.5 towards nearest even integer.";
	};

	// trunc
	class MpOperationTrunc : public MpOperationFunc<double, double>
	{
	public:
		MpOperationTrunc();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// frac
	class MpOperationFrac : public MpOperationFunc<double, double>
	{
	public:
		MpOperationFrac();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// floor
	class MpOperationFloor : public MpOperationFunc<double, double>
	{
	public:
		MpOperationFloor();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// ceil
	class MpOperationcCeil : public MpOperationFunc<double, double>
	{
	public:
		MpOperationcCeil();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// Addition
	template<typename T>
	class MpOperationAdd : public MpOperationBinary<T>
	{
	public:
		MpOperationAdd();

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T calculate(T vl, T vr) const override;
	};

	// Subtraction
	template<typename T>
	class MpOperationSub : public MpOperationBinary<T>
	{
	public:
		MpOperationSub();

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T calculate(T vl, T vr) const override;
	};

	// Multiplication
	template<typename T>
	class MpOperationMul : public MpOperationBinary<T>
	{
	public:
		MpOperationMul();

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T calculate(T vl, T vr) const override;
	};

	// Division
	template<typename T>
	class MpOperationDiv : public MpOperationBinary<T>
	{
	public:
		MpOperationDiv();

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T calculate(T vl, T vr) const override;
	};

	// Minimum
	class MpOperationMin : public MpOperationBinary<double>
	{
	public:
		MpOperationMin() :
			MpOperationBinary<double>(Signature(2), 0)
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
		MpOperationMax() :
			MpOperationBinary<double>(Signature(2), 0)
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
		MpOperationEq();

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T calculate(T vl, T vr) const override;
	};

	// Inequality
	template<typename T>
	class MpOperationNe : public MpOperationBinary<T>
	{
	public:
		MpOperationNe();

	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T calculate(T vl, T vr) const override;
	};

	// Lesser than
	class MpOperationLt : public  MpOperationBinary<double>
	{
	public:
		MpOperationLt() :
			MpOperationBinary<double>(Signature(2, Signature::type::real, MpOperationFlags::OpFlagIsOperator), 8)
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
		MpOperationLe() :
			MpOperationBinary<double>(Signature(2, Signature::type::real, MpOperationFlags::OpFlagIsOperator), 8)
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
		MpOperationGt() :
			MpOperationBinary<double>(Signature(2, Signature::type::real, MpOperationFlags::OpFlagIsOperator), 8)
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
		MpOperationGe() :
			MpOperationBinary<double>(Signature(2, Signature::type::real, MpOperationFlags::OpFlagIsOperator), 8)
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
		MpOperationModulo() :
			MpOperationBinary<double>(Signature(2, Signature::type::real,  MpOperationFlags::OpFlagIsOperator), 5)
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
		MpOperationTernary();

		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) const override;
	protected:
		bool isColon_;
	};

	// Assignment
	template<typename T>
	class MpOperationAssignment : public MpOperation
	{
	public:
		MpOperationAssignment();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) const override;
	};
}