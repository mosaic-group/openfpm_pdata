//
// Created by landfried on 24.01.22.
//

#ifndef OPENFPM_PDATA_OPERATIONPROXY_HPP
#define OPENFPM_PDATA_OPERATIONPROXY_HPP


#include <string>
#include <Vector/vector_dist.hpp>




template <typename T>
class OperationProxy {
    T data;
public:
    OperationProxy(T data_in) : data(data_in) {}
    //    OperationProxy(Point<n1, T> &point_in) : data(point_in.data) {}


/*    OperationProxy<T> operator+= (OperationProxy<T> rhs) {
        data += rhs.data;
        return *this;
    }

    OperationProxy<T> operator+= (T rhs) {
        data += rhs;
        return *this;
    }

    T operator* (T rhs) {
        rhs *= data;
        return rhs;
    }*/
};

/*template <typename T, unsigned int n1>
class OperationProxy<Point<n1, T>> {
    Point<n1, T> &data;

public:

    OperationProxy(Point<n1, T> &data_in) : data(data_in) {}



};*/


template <typename T, unsigned int n1>
class OperationProxy<T[n1]> {

    T (&data)[n1];


public:

    OperationProxy(T (&data_in)[n1]) : data(data_in) {}

    T& operator[] (int c) {
        return data[c];
    }

    // Assignment operators

    // OperationProxy

    OperationProxy<T[n1]> operator= (OperationProxy<T[n1]> rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] = rhs.data[i];
        }
        return *this;
    }

    OperationProxy<T[n1]> operator+= (OperationProxy<T[n1]> rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] += rhs.data[i];
        }
        return *this;
    }

    OperationProxy<T[n1]> operator-= (OperationProxy<T[n1]> rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] -= rhs.data[i];
        }
        return *this;
    }

    OperationProxy<T[n1]> operator*= (OperationProxy<T[n1]> rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] *= rhs.data[i];
        }
        return *this;
    }

    OperationProxy<T[n1]> operator/= (OperationProxy<T[n1]> rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] /= rhs.data[i];
        }
        return *this;
    }


    // Point

    OperationProxy<T[n1]> operator= (Point<n1, T> rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] = rhs[i];
        }
        return *this;
    }

    OperationProxy<T[n1]> operator+= (Point<n1, T> rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] += rhs[i];
        }
        return *this;
    }

    OperationProxy<T[n1]> operator-= (Point<n1, T> rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] -= rhs[i];
        }
        return *this;
    }

    OperationProxy<T[n1]> operator*= (Point<n1, T> rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] *= rhs[i];
        }
        return *this;
    }

    OperationProxy<T[n1]> operator/= (Point<n1, T> rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] /= rhs[i];
        }
        return *this;
    }




    // Primitive types

    OperationProxy<T[n1]> operator= (T rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] = rhs;
        }
        return *this;
    }

    OperationProxy<T[n1]> operator+= (T rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] += rhs;
        }
        return *this;
    }

    OperationProxy<T[n1]> operator-= (T rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] -= rhs;
        }
        return *this;
    }

    OperationProxy<T[n1]> operator*= (T rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] *= rhs;
        }
        return *this;
    }

    OperationProxy<T[n1]> operator/= (T rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] /= rhs;
        }
        return *this;
    }



    // Mathematical operations: +, - *, /


    // OperationProxy

    Point<n1, T> operator+ (OperationProxy<T[n1]> rhs) {
        Point<n1, T> result(data);
        for (int i = 0; i < n1; ++i) {
            result[i] += rhs.data[i];
        }
        return result;
    }

    Point<n1, T> operator- (OperationProxy<T[n1]> rhs) {
        Point<n1, T> result(data);
        for (int i = 0; i < n1; ++i) {
            result[i] -= rhs.data[i];
        }
        return result;
    }

    Point<n1, T> operator* (OperationProxy<T[n1]> rhs) {
        Point<n1, T> result(data);
        for (int i = 0; i < n1; ++i) {
            result[i] *= rhs.data[i];
        }
        return result;
    }

    Point<n1, T> operator/ (OperationProxy<T[n1]> rhs) {
        Point<n1, T> result(data);
        for (int i = 0; i < n1; ++i) {
            result[i] /= rhs.data[i];
        }
        return result;
    }


    // Point

    Point<n1, T> operator+ (Point<n1, T> rhs) {
        for (int i = 0; i < n1; ++i) {
            rhs[i] += data[i];
        }
        return rhs;
    }

    Point<n1, T> operator- (Point<n1, T> rhs) {
        for (int i = 0; i < n1; ++i) {
            rhs[i] -= data[i];
        }
        return rhs;
    }

    Point<n1, T> operator* (Point<n1, T> rhs) {
        for (int i = 0; i < n1; ++i) {
            rhs[i] *= data[i];
        }
        return rhs;
    }

    Point<n1, T> operator/ (Point<n1, T> rhs) {
        for (int i = 0; i < n1; ++i) {
            rhs[i] /= data[i];
        }
        return rhs;
    }


    // Primitive types

    Point<n1, T> operator+ (T rhs) {
        Point<n1, T> result(data);
        for (int i = 0; i < n1; ++i) {
            result[i] += rhs;
        }
        return result;
    }

    Point<n1, T> operator- (T rhs) {
        Point<n1, T> result(data);
        for (int i = 0; i < n1; ++i) {
            result[i] -= rhs;
        }
        return result;
    }

    Point<n1, T> operator* (T rhs) {
        Point<n1, T> result(data);
        for (int i = 0; i < n1; ++i) {
            result[i] *= rhs;
        }
        return result;
    }

    Point<n1, T> operator/ (T rhs) {
        Point<n1, T> result(data);
        for (int i = 0; i < n1; ++i) {
            result[i] /= rhs;
        }
        return result;
    }





    std::string toString() const {
        std::string output = "{ ";
        for (int i = 0; i < n1; ++i) {
            output.append(std::to_string(data[i]));
            output.append(" ");
        }
        output.append("}");
        return output;
    }
};

// Extern overloads


/*
// Operator overloads for Type T

template <typename T, unsigned int n1>
Point<n1, T> operator* (T lhs, OperationProxy<T[n1]> rhs) {
    return rhs * lhs;
}
*/



// Operator overloads for Point class


template <typename T, unsigned int n1>
Point<n1, T> operator* (Point<n1, T> lhs, OperationProxy<T[n1]> rhs) {
    for (int i = 0; i < n1; ++i) {
        lhs[i] *= rhs[i];
    }
    return lhs;
}

template <typename T, unsigned int n1>
Point<n1, T> operator/ (Point<n1, T> lhs, OperationProxy<T[n1]> rhs) {
    for (int i = 0; i < n1; ++i) {
        lhs[i] /= rhs[i];
    }
    return lhs;
}



#endif //OPENFPM_PDATA_OPERATIONPROXY_HPP
