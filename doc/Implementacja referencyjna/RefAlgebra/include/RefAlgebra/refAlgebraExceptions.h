#pragma once
#include <string>

#include "types.h"


namespace RefAlgebra
{
	class RefAlgebraException
	{
	public:
		RefAlgebraException(const std::string& className, const std::string& operation);

		virtual std::string what() const;

		const std::string& className() const;
		const std::string& operation() const;

	protected:
		std::string clssNme;
		std::string oper;

		static std::string toString(Shape shape);
	};

	class DimensionMismatchException : public RefAlgebraException
	{
	public:
		DimensionMismatchException(const std::string& className, const std::string& operation, Shape lhsShape, Shape rhsShape);

		virtual std::string what() const override;

		Shape lhsShape() const;
		Shape rhsShape() const;

	protected:
		Shape lhsShp;
		Shape rhsShp;
	};

	class AdditionException : public DimensionMismatchException
	{
	public:
		AdditionException(const std::string& className, Shape lhsShape, Shape rhsShape);
	};

	class SubstractionException : public DimensionMismatchException
	{
	public:
		SubstractionException(const std::string& className, Shape lhsShape, Shape rhsShape);
	};

	class MultiplicationException : public DimensionMismatchException
	{
	public:
		MultiplicationException(const std::string& className, Shape lhsShape, Shape rhsShape);

		virtual std::string what() const override;
	};

	class ExponentationException : public RefAlgebraException
	{
	public:
		ExponentationException(Shape shape);

		virtual std::string what() const override;

		Shape shape() const;

	protected:
		Shape shp;
	};
}

