#include <iterator>
#include <iostream>
#include <cstddef>
#include <vector>


template <typename T>
class Polynomial {
 private:
    std::vector<T> C;

 public:
    void normalize() {
        while (!C.empty() && C.back() == T()) {
            C.pop_back();
        }
    }

    std::vector<T> getVector() const {
        return C;
    }

    Polynomial(const std::vector<T>& coefficients) {
        C = coefficients;
        normalize();
    }

    Polynomial() {
        C.push_back(T());
    }

    Polynomial(T free_coefficient) {
        C.push_back(free_coefficient);
        normalize();
    }

    Polynomial(T old_coefficient, int pos, int) {
        for (int j = 1; j <= pos; ++j) {
            C.push_back(T(0));
        }
        C.push_back(old_coefficient);
    }

    template <typename Iter>
    Polynomial(Iter start, Iter finish) {
        while (start != finish) {
            C.push_back(*start);
            ++start;
        }
        normalize();
    }

    T getCoefficient(size_t idx) const {
        if (idx >= 0 && idx < C.size()) {
            return C[idx];
        }
        return T();
    }

    T operator[](size_t deg) const {
        if (deg >= C.size()) {
            return T();
        }
        return C[deg];
    }

    int Degree() const {
        for (int j = C.size() - 1; j >= 0; j--) {
            if (C[j] != T()) {
                return j;
            }
        }
        return -1;
    }

    size_t getSize() const {
        return C.size();
    }

    Polynomial<T>& operator += (const Polynomial<T>& rhs) {
        size_t rhsSz = rhs.getSize();
        size_t max_size = std:: max(rhsSz, getSize());
        C.resize(max_size);
        for (size_t j = 0; j < rhsSz; ++j) {
            C[j] += rhs.getCoefficient(j);
        }
        normalize();
        return *this;
    }

    Polynomial<T> operator + (const Polynomial<T>& rhs) const {
        size_t rhsSz = rhs.getSize();
        std::vector<T> new_C;
        size_t max_size = std:: max(rhsSz, getSize());
        new_C.resize(max_size);
        for (size_t j = 0; j < max_size; ++j) {
            new_C[j] = getCoefficient(j) + rhs.getCoefficient(j);
        }
        return Polynomial<T>(new_C);
    }

    Polynomial<T>& operator -= (const Polynomial<T>& rhs) {
        size_t rhsSz = rhs.getSize();
        C.resize(std::max(rhsSz, getSize()));
        for (size_t j = 0; j < rhsSz; ++j) {
            C[j] -= rhs.getCoefficient(j);
        }
        normalize();
        return *this;
    }

    Polynomial<T> operator -(const Polynomial<T>& rhs) const {
        size_t rhsSz = rhs.getSize();
        std::vector<T> new_C;
        size_t max_size = std:: max(rhsSz, getSize());
        new_C.resize(max_size);
        for (size_t j = 0; j < max_size; ++j) {
            new_C[j] = getCoefficient(j) - rhs.getCoefficient(j);
        }
        return Polynomial<T>(new_C);
    }

    Polynomial<T>& operator *= (const Polynomial<T>& rhs) {
        size_t lhsSz = getSize(), rhsSz = rhs.getSize();
        std::vector<T> new_C(lhsSz + rhsSz);
        for (size_t i = 0; i < lhsSz; ++i) {
            for (size_t j = 0; j < rhsSz; ++j) {
                new_C[i+j] += getCoefficient(i) * rhs.getCoefficient(j);
            }
        }
        C = new_C;
        normalize();
        return *this;
    }

    Polynomial<T> operator *(const Polynomial<T>& rhs) const {
        size_t lhsSz = getSize(), rhsSz = rhs.getSize();
        std::vector<T> new_C(lhsSz + rhsSz);
        for (size_t i = 0; i < lhsSz; ++i) {
            for (size_t j = 0; j < rhsSz; ++j) {
                new_C[i+j] += getCoefficient(i) * rhs.getCoefficient(j);
            }
        }
        return Polynomial<T>(new_C);
    }

    typename std::vector<T>::const_iterator begin() const {
        return C.begin();
    }

    typename std::vector<T>::const_iterator end() const {
        return C.end();
    }

    T operator()(T x) const {
        T res = T();
        for (int j = C.size() - 1; j >= 0; --j) {
            res = res * x + C[j];
        }
        return res;
    }

    bool operator ==(const Polynomial<T>& rhs) const {
        int i = Degree();
        int j = rhs.Degree();
        if (i != j) {
            return false;
        }
        while (i >= 0) {
            if (getCoefficient(i) != rhs.getCoefficient(j)) {
                return false;
            }
            --i;
            --j;
        }
        return true;
    }



    bool operator !=(const Polynomial& rhs) const {
        return !((*this) == rhs);
    }

    Polynomial<T> operator&(const Polynomial<T>& rhs) const {
        Polynomial<T> res(T(0));
        for (size_t i = 0; i < C.size(); ++i) {
            if (C[i] != T()) {
                Polynomial<T> tmp(T(1));
                for (size_t j = 1; j <= i; ++j) {
                    tmp *= rhs;
                }
                tmp *= C[i];
                res += tmp;
            }
        }
        return res;
    }

    Polynomial<T> operator/(const Polynomial<T>&rhs) const {
        Polynomial<T> num = *this, den = rhs;
        std::vector <T> c_res(num.getSize());
        int maxDegreeOfRhs = den.Degree();
        while (1) {
            int curMaxDegree = num.Degree();
            if (curMaxDegree < maxDegreeOfRhs) {
                break;
            }
            T quot = num[curMaxDegree] / den[maxDegreeOfRhs];
            num -= (Polynomial(quot, curMaxDegree - maxDegreeOfRhs, 0) * den);
            c_res[curMaxDegree - maxDegreeOfRhs] += quot;
        }
        return Polynomial<T>(c_res);
    }

    Polynomial<T> operator%(const Polynomial<T>&rhs) const {
        return (*this) - (((*this)/rhs) * rhs);
    }

    Polynomial<T>& operator /=(const Polynomial<T>&rhs) {
        Polynomial<T> X = (*this) / rhs;
        (*this) = X;
        return (*this);
    }
    Polynomial<T>& operator %=(const Polynomial<T>&rhs) {
        Polynomial<T> X = (*this) % rhs;
        (*this) = X;
        return (*this);
    }

    Polynomial<T> operator, (const Polynomial<T>&rhs) const {
        Polynomial<T> A = (*this), B = rhs;
        if (A.Degree() > B.Degree()) {
            std:: swap(A, B);
        }
        while (A.Degree() > 0) {
            Polynomial<T> tmp = A;
            B %= A;
            A = B;
            B = tmp;
        }
        if (A != T(0)) {
            return Polynomial(T(1));
        }
        B /= B[B.Degree()];
        return B;
    }
};

template <typename T>
std::ostream& operator <<(std::ostream& out, const Polynomial<T>& A) {
    int sz = A.getSize();
    bool first = true, somethingPrinted = false;
    for (int i = sz - 1; i >= 0; --i) {
        if (A.getCoefficient(i) == T()) {
            continue;
        }
        somethingPrinted = true;
        T currentCoef = A.getCoefficient(i);
        if (first) {
            if (currentCoef != T(1) && currentCoef != T(-1)) {
                if (i > 1) {
                    out << currentCoef << "*x^" << i;
                } else if (i == 1) {
                    out << currentCoef << "*x";
                } else {
                    out << currentCoef;
                }
            } else {
                if (i == 0) {
                    out << currentCoef;
                } else {
                    if (currentCoef == T(-1)) {
                        out << "-";
                    }
                    if (i == 1) {
                        out << "x";
                    } else {
                        out << "x^" << i;
                    }
                }
            }
            first = false;
        } else {
            if (currentCoef > T()) {
                out << "+";
                if (currentCoef == T(1)) {
                    if (i == 0) {
                        out << "1";
                    } else if (i == 1) {
                        out << "x";
                    } else {
                        out << "x^" << i;
                    }
                } else {
                    if (i == 0) {
                        out << currentCoef;
                    } else if (i == 1) {
                        out << currentCoef << "*x";
                    } else {
                        out << currentCoef << "*x^" << i;
                    }
                }
            } else {
                if (currentCoef == T(-1)) {
                    if (i == 0) {
                        out << "-1";
                    } else if (i == 1) {
                        out << "-x";
                    } else {
                        out << "-x^" << i;
                    }
                } else {
                    if (i == 0) {
                        out << currentCoef;
                    } else if (i == 1) {
                        out << currentCoef << "*" << "x";
                    } else {
                        out << currentCoef << "*x^" << i;
                    }
                }
            }
        }
    }
    if (!somethingPrinted) {
        out << "0";
    }
    return out;
}
