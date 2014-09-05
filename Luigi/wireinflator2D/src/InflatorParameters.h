#ifndef CELLPARAMETERS_H
#define CELLPARAMETERS_H

#include <utility>
#include <map>

#include <vcg/space/point2.h>

struct TessellationParameters
{
	double max_area;
	double min_angle;

	TessellationParameters(void)
	{
		max_area = 0.0002;
		min_angle = 30;
	}
};

struct ParameterChange
{
	typedef int NodeID;

	typedef enum
	{
		Translation,
		Radius
	} OperationType;

	OperationType                type;
	std::map<NodeID, vcg::Point2d> node_ops;
};

template <size_t NumParameters>
struct CellParameters
{
	enum { NumberOfParameters = NumParameters};

	virtual ~CellParameters(void)
	{
		;
	}

	/*!
	 * \brief numberOfParameters
	 * \return the number of pattern parameters
	 */
	virtual size_t numberOfParameters(void) const
	{
		return NumParameters;
	}

	/*!
	 * \brief parameter
	 * \param i the parameter index.
	 * \return the i-th parameter
	 */
	virtual double & parameter(int i) = 0;

	virtual const double & cParameter(int i) const = 0;

	/*!
	 * \brief parameterRange return the admissible range for the specified parameter
	 * \param i the parameter index.
	 * \return the pair <min value, max value> for the i-th parameter.
	 */
	virtual std::pair<double,double> parameterRange(int i) const = 0;

	/*!
	 * \brief isValid checks wheter the current set of parameters is valid.
	 * \return true if the parameters are valide, false otherwise.
	 */
	virtual bool isValid(void) const
	{
		for (size_t i=0; i<numberOfParameters(); i++)
		{
			const double & par = cParameter(i);
			std::pair<double, double> pr = this->parameterRange(i);
			if (par < pr.first || par > pr.second)
				return false;
		}
		return true;
	}
};

#endif // CELLPARAMETERS_H
