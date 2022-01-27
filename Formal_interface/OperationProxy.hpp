//
// Created by landfried on 24.01.22.
//

#ifndef OPENFPM_PDATA_OPERATIONPROXY_HPP
#define OPENFPM_PDATA_OPERATIONPROXY_HPP



template <typename T>
class OperationProxy {
    T data;
public:
    OperationProxy(T data) : data(data) {}
};


template <typename T, unsigned int n1>
class OperationProxy<T[n1]> {

    T (&data)[n1];


public:
    OperationProxy(T (&data)[n1]) : data(data) {}

    OperationProxy<T[n1]> operator+= (OperationProxy<T[n1]> rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] += rhs.data[i];
        }
        return *this;
    }

    OperationProxy<T[n1]> operator+= (Point<n1, T> rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] += rhs[i];
        }
        return *this;
    }

    OperationProxy<T[n1]> operator= (T rhs) {
        for (int i = 0; i < n1; ++i) {
            data[i] = rhs;
        }
        return *this;
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

template <typename T>
class OperationProxyValue {
    T data;
public:
    OperationProxyValue(T data) : data(data)  {}
};

template <typename T, unsigned int n1>
class OperationProxyValue<T[n1]> {

    T data[n1];

public:
    OperationProxyValue(T data[n1]) {
        for (int i = 0; i < n1; ++i) {
            this->data[i] = data[i];
        }
    }


    OperationProxyValue<T[n1]> operator* (T rhs) {
        OperationProxyValue<T[n1]> result(data);
        for (int i = 0; i < n1; ++i) {
            result.data[i] *= rhs;
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




#endif //OPENFPM_PDATA_OPERATIONPROXY_HPP
