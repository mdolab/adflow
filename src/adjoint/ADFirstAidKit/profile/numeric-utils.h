#ifndef NUMERIC_UTILS_H
#define NUMERIC_UTILS_H 1

#include <limits>

template <typename T> bool is_nan(T number) {
        using namespace std;
        if (numeric_limits<T>::is_integer)
                return false;
        return ! (number == number);
}

template <int N, typename T> struct _power {
        T value(T a) {
                return a * _power<N-1,T>(a);
        }
};
template <typename T> struct _power<0,T> {
        T value(T a) {
                return 1;
        }
};

template <int N, typename T> T power(T a) {
        return _power<N,T>().value(a);
}
        
template <typename T> bool is_finite(T number) {
        using namespace std;
        if (numeric_limits<T>::is_integer)
                return true;
        const T inf = numeric_limits<T>::infinity();
        if (number == inf || is_nan(number)) // NaN detection
                        return false;
        return true;
};

#endif
