/*
 * vector_dist_operators.hpp
 *
 *  Created on: Jun 11, 2016
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_VECTOR_DIST_OPERATORS_HPP_
#define SRC_VECTOR_VECTOR_DIST_OPERATORS_HPP_

#define SUM 1
#define SUB 2
#define MUL 3
#define DIV 4

template <typename exp1, typename exp2, unsigned int op>
class vector_dist_expression_op
{

};

template <typename exp1, typename exp2>
class vector_dist_expression_op<exp1,exp2,SUM>
{
	const exp1 o1;
	const exp2 o2;

public:

	typedef decltype(o1.value(vect_dist_key_dx(0))) r_type;

	vector_dist_expression_op(const exp1 & o1, const exp2 & o2)
	:o1(o1),o2(o2)
	{}

	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx(0)))>::type > inline r_type value(const vect_dist_key_dx & key) const
	{
		return o1.value(key) + o2.value(key);
	}
};

template <typename exp1, typename exp2>
class vector_dist_expression_op<exp1,exp2,SUB>
{
	const exp1 o1;
	const exp2 o2;

public:

	vector_dist_expression_op(const exp1 & o1, const exp2 & o2)
	:o1(o1),o2(o2)
	{}

	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx(0)))>::type > inline r_type value(const vect_dist_key_dx & key) const
	{
		return o1.value(key) - o2.value(key);
	}

};

template <typename exp1, typename exp2>
class vector_dist_expression_op<exp1,exp2,MUL>
{
	const exp1 o1;
	const exp2 o2;

public:

	vector_dist_expression_op(const exp1 & o1, const exp2 & o2)
	:o1(o1),o2(o2)
	{}

	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx(0)))>::type > inline r_type value(const vect_dist_key_dx & key) const
	{
		return o1.value(key) * o2.value(key);
	}

};

template <typename exp1, typename exp2>
class vector_dist_expression_op<exp1,exp2,DIV>
{
	const exp1 o1;
	const exp2 o2;

public:

	vector_dist_expression_op(const exp1 & o1, const exp2 & o2)
	:o1(o1),o2(o2)
	{}

	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx(0)))>::type > inline r_type value(const vect_dist_key_dx & key) const
	{
		return o1.value(key) / o2.value(key);
	}

};

/*! \brief Main class that encapsulate a vector operation expression
 *
 * \param prp property involved
 * \param exp expression type
 *
 */
template<unsigned int prp, typename vector>
class vector_dist_expression
{
	vector & v;

public:

	vector_dist_expression(vector & v)
	:v(v)
	{}

	inline auto value(const vect_dist_key_dx & k) const -> decltype(v.template getProp<prp>(k))
	{
		return v.template getProp<prp>(k);
	}

	/*! \brief evaluate the expression and store it in the vector
	 *
	 * \param v_exp expression to evaluate
	 *
	 */
	template<typename exp1, typename exp2, unsigned int op> vector & operator=(const vector_dist_expression_op<exp1,exp2,op> & v_exp)
	{
		auto it = v.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			v.template getProp<prp>(key) = v_exp.value(key);

			++it;
		}

		return v;
	}

	/*! \brief evaluate the expression and store it in the vector
	 *
	 * \param v_exp expression to evaluate
	 *
	 */
	vector & operator=(double d)
	{
		auto it = v.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			v.template getProp<prp>(key) = d;

			++it;
		}

		return v;
	}
};

/*! \brief Main class that encapsulate a vector operation expression
 *
 * \param prp property involved
 * \param exp expression type
 *
 */
template<unsigned int prp>
class vector_dist_expression<prp,double>
{
	double d;

public:

	inline vector_dist_expression(double & d)
	:d(d)
	{}

	inline double value(const vect_dist_key_dx & k) const
	{
		return d;
	}
};


/* \brief sum two distributed vectors
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<unsigned int p1, unsigned int p2, typename v1, typename v2>
inline vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression<p2,v2>,SUM>
operator+(const vector_dist_expression<p1,v1> & va, const vector_dist_expression<p2,v2> & vb)
{
	vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression<p2,v2>,SUM> exp_sum(va,vb);

	return exp_sum;
}

/* \brief sum two distributed vectors
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2, unsigned int op1, unsigned int prp1, typename v1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<prp1,v1>,SUM>
operator+(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression<prp1,v1> & vb)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<prp1,v1>,SUM> exp_sum(va,vb);

	return exp_sum;
}

/* \brief sum two distributed vectors
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2, unsigned int op1, unsigned int prp1, typename v1>
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression_op<exp1,exp2,op1>,SUM>
operator+(const vector_dist_expression<prp1,v1> & va, const vector_dist_expression_op<exp1,exp2,op1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression_op<exp1,exp2,op1>,SUM> exp_sum(va,vb);

	return exp_sum;
}

/* \brief sum two distributed vectors
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2, unsigned int op1, typename exp3 , typename exp4, unsigned int op2>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,SUM>
operator+(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression_op<exp3,exp4,op2> & vb)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,SUM> exp_sum(va,vb);

	return exp_sum;
}

/* \brief sum two distributed vectors
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<unsigned int prp1 , typename v1>
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,double>,SUM>
operator+(const vector_dist_expression<prp1,v1> & va, double d)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,double>,SUM> exp_sum(va,vector_dist_expression<0,double>(d));

	return exp_sum;
}

/* \brief sum two distributed vectors
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<unsigned int prp1 , typename v1>
inline vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<prp1,v1>,SUM>
operator+(double d, const vector_dist_expression<prp1,v1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<prp1,v1>,SUM> exp_sum(vector_dist_expression<0,double>(d),vb);

	return exp_sum;
}

/* \brief sum two distributed vectors
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,SUM>
operator+(const vector_dist_expression_op<exp1,exp2,op1> & va, double d)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,SUM> exp_sum(va,vector_dist_expression<0,double>(d));

	return exp_sum;
}


/* \brief sum two distributed vectors
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<unsigned int p1, unsigned int p2, typename v1, typename v2>
inline vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression<p2,v2>,SUB>
operator-(const vector_dist_expression<p1,v1> & va, const vector_dist_expression<p2,v2> & vb)
{
	vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression<p2,v2>,SUB> exp_sum(va,vb);

	return exp_sum;
}


/* \brief sum two distributed vectors
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1, unsigned int p2, typename v2>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<p2,v2>,SUB>
operator-(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression<p2,v2> & vb)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<p2,v2>,SUB> exp_sum(va,vb);

	return exp_sum;
}

/* \brief sum two distributed vectors
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1, unsigned int p2, typename v2>
inline vector_dist_expression_op<vector_dist_expression<p2,v2>,vector_dist_expression_op<exp1,exp2,op1>,SUB>
operator-(const vector_dist_expression<p2,v2> & va, const vector_dist_expression_op<exp1,exp2,op1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<p2,v2>, vector_dist_expression_op<exp1,exp2,op1>,SUB> exp_sum(va,vb);

	return exp_sum;
}

/* \brief sum two distributed vectors
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1, typename exp3, typename exp4, unsigned int op2>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,SUB>
operator-(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression_op<exp3,exp4,op2> & vb)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,SUB> exp_sum(va,vb);

	return exp_sum;
}

/* \brief sum two distributed vectors
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<unsigned int prp1, typename v1>
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,double>,SUB>
operator-(const vector_dist_expression<prp1,v1> & va, double d)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,double>,SUB> exp_sum(va,vector_dist_expression<0,double>(d));

	return exp_sum;
}

/* \brief sum two distributed vectors
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<unsigned int prp1, typename v1>
inline vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<prp1,v1>,SUB>
operator-(double d, const vector_dist_expression<prp1,v1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<prp1,v1>,SUB> exp_sum(vector_dist_expression<0,double>(d),vb);

	return exp_sum;
}

/* \brief multiplication double with a vector expression
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<unsigned int p2, typename v2>
inline vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<p2,v2>,MUL>
operator*(double d, const vector_dist_expression<p2,v2> & vb)
{
	vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<p2,v2>,MUL> exp_sum(vector_dist_expression<0,double>(d),vb);

	return exp_sum;
}

/* \brief multiplication double with a vector expression
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<unsigned int p2, typename v2>
inline vector_dist_expression_op<vector_dist_expression<p2,v2>,vector_dist_expression<0,double>,MUL>
operator*(const vector_dist_expression<p2,v2> & va, double d)
{
	vector_dist_expression_op<vector_dist_expression<p2,v2>,vector_dist_expression<0,double>,MUL> exp_sum(va,vector_dist_expression<0,double>(d));

	return exp_sum;
}

/* \brief multiplication double with a vector expression
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<unsigned int p1, typename v1,unsigned int p2, typename v2>
inline vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression<p2,v2>,MUL>
operator*(const vector_dist_expression<p1,v1> & va, const vector_dist_expression<p2,v2> & vb)
{
	vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression<p2,v2>,MUL> exp_sum(va,vb);

	return exp_sum;
}

/* \brief multiplication double with a vector expression
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<unsigned int p1, typename v1, typename exp1, typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression_op<exp1,exp2,op1>,MUL>
operator*(const vector_dist_expression<p1,v1> & va, const vector_dist_expression_op<exp1,exp2,op1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression_op<exp1,exp2,op1>,MUL> exp_sum(va,vb);

	return exp_sum;
}

/* \brief multiplication double with a vector expression
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<unsigned int p1, typename v1, typename exp1, typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<p1,v1>,MUL>
operator*(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression<p1,v1> & vb)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<p1,v1>,MUL> exp_sum(va,vb);

	return exp_sum;
}

/* \brief multiplication double with a vector expression
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1, typename exp3 , typename exp4, unsigned int op2>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,MUL>
operator*(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression_op<exp3,exp4,op2> & vb)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,MUL> exp_sum(va,vb);

	return exp_sum;
}

/* \brief multiplication double with a vector expression
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,DIV>
operator/(const vector_dist_expression_op<exp1,exp2,op1> & va, double d)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,DIV> exp_sum(va,vector_dist_expression<0,double>(d));

	return exp_sum;
}


/* \brief multiplication double with a vector expression
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,DIV>
operator/(double d, const vector_dist_expression_op<exp1,exp2,op1> & va)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,DIV> exp_sum(vector_dist_expression<0,double>(d),va);

	return exp_sum;
}

/* \brief multiplication double with a vector expression
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<unsigned int prp1, typename v1>
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,double>,DIV>
operator/(const vector_dist_expression<prp1,v1> & va, double d)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,double>,DIV> exp_sum(va,vector_dist_expression<0,double>(d));

	return exp_sum;
}

/* \brief multiplication double with a vector expression
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<unsigned int prp1, typename v1>
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,double>,DIV>
operator/(double d, const vector_dist_expression<prp1,v1> & va)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,double>,DIV> exp_sum(va,vector_dist_expression<0,double>(d));

	return exp_sum;
}

/* \brief multiplication double with a vector expression
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<unsigned int prp1, typename v1, unsigned int prp2, typename v2>
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<prp2,v2>,DIV>
operator/(const vector_dist_expression<prp1,v1> & va, const vector_dist_expression<prp2,v2> & vb)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<prp2,v2>,DIV> exp_sum(va,vb);

	return exp_sum;
}

/* \brief multiplication double with a vector expression
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<unsigned int prp1, typename v1, typename exp1,typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression_op<exp1,exp2,op1>,DIV>
operator/(const vector_dist_expression<prp1,v1> & va, const vector_dist_expression_op<exp1,exp2,op1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression_op<exp1,exp2,op1>,DIV> exp_sum(va,vb);

	return exp_sum;
}

/* \brief multiplication double with a vector expression
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<unsigned int prp1, typename v1, typename exp1,typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<prp1,v1>,DIV>
operator/(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression<prp1,v1> & vb)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<prp1,v1>,DIV> exp_sum(va,vb);

	return exp_sum;
}

/* \brief multiplication double with a vector expression
 *
 * \param v1 vector1
 * \param v2 vector2
 *
 * \return a grid_key_dx_expression that encapsulate the expression
 *
 */
template<typename exp1,typename exp2, unsigned int op1, typename exp3, typename exp4, unsigned int op2>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,DIV>
operator/(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression_op<exp3,exp4,op2> & vb)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,DIV> exp_sum(va,vb);

	return exp_sum;
}

#endif /* SRC_VECTOR_VECTOR_DIST_OPERATORS_HPP_ */
