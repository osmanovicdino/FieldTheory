#ifndef LOGSPACE_H
#define LOGSPACE_H

class LogSpace {
private:
    double num; //our number is given by e^num
    bool iszero;
    short sign;
    
public:

    LogSpace() {
        num = 0.0;
        iszero = false;
        sign = 1;
    }

    template <class T>
    LogSpace(T a) {
        if(a<0) {
            num = log(-a);
            iszero = false;
            sign = -1;
        }
        else if(a==T(0)) {
            num =1.0;
            iszero = true;
            sign = 1;
        }
        else{
        num = log(a);
        iszero = false;
        sign = 1;
        }
    }

    LogSpace(const LogSpace &a) {
        num = a.num;
        iszero = a.iszero;
        sign = a.sign;
    }

    void setdouble(double numm) {
        num = numm;
    }

    double getdouble() {
        if(iszero) return 0.0;
        return sign*exp(this->num);
    }

    inline LogSpace& operator*=(const LogSpace &a) {
        num += (a.num);
        sign *= (a.sign);
        iszero = iszero || a.iszero;
        return *this;
    }
    inline LogSpace& operator/=(const LogSpace &a)
    {
        if(a.iszero) error("division by zero");
        num -= (a.num);
        sign *= a.sign;
        return *this;
    }
    inline LogSpace& operator+=(const LogSpace &a)
    {
        if(a.iszero) {
            return *this;
        }
        else if(this->iszero) {
            num = a.num;
            sign =a.sign;
            iszero = a.iszero;
            return *this;
        }
        else{
            if(sign ==a.sign) {
                num += log(1 + exp(a.num - num));
                iszero = false;
            }
            else if(sign ==1 && a.sign ==-1){
                if(num > a.num) {
                num += log(1 - exp(a.num - num));
                iszero = false;
                }
                else{
                    num = a.num + log(1 - exp(num - a.num));
                    sign =-1;
                    iszero = false;
                }
            }
            else {
                if (num > a.num)
                {
                    num += log(1 - exp(a.num - num));
                    iszero = false;
                }
                else{
                    num = a.num + log(1 - exp(num - a.num));
                    sign = 1;
                    iszero = false;
                }
            }
            return *this;

        }
    }

    inline LogSpace &operator-=(const LogSpace &a)
    {
        *this += (-a);
    }


    LogSpace& operator=(const LogSpace &v)
    { //assigment operator

        num = v.num;
        iszero = v.iszero;
        sign = v.sign;
        return *this;
    }

    template <class T>
    LogSpace& operator<<(const T &a) {
            if (a < 0)
            {
                num = log(-a);
                iszero = false;
                sign = -1;
            }
            else if (a == 0)
            {
                num = 1.0;
                iszero = true;
                sign = 1;
            }
            else
            {
                num = log(a);
                iszero = false;
                sign = 1;
            }
            return *this;
    }

    // template <class T>
    // LogSpace &operator=(const T &a)
    // { //assigment operator
    //     if (a < 0)
    //     {
    //         num = log(-a);
    //         iszero = false;
    //         sign = -1;
    //     }
    //     else if (a == 0)
    //     {
    //         num = 1.0;
    //         iszero = true;
    //         sign = 1;
    //     }
    //     else
    //     {
    //         num = log(a);
    //         iszero = false;
    //         sign = 1;
    //     }
    // }

    friend LogSpace operator*(const LogSpace &a,const LogSpace &b) {
    LogSpace res(a);
    res *= b;
    return res;
    }

    friend LogSpace operator/(const LogSpace &a, const LogSpace &b)
    {
        LogSpace res(a);
        res /= b;
        return res;
    }
    friend LogSpace operator+(const LogSpace &a, const LogSpace &b)
    {
        LogSpace res(a);
        res += b;
        return res;
    }
    friend LogSpace operator-(const LogSpace &a, const LogSpace &b)
    {
        LogSpace res(a);
        res -= b;
        return res;
    }

    friend LogSpace pow(LogSpace &a, LogSpace &b) {
        LogSpace res;
        res.num = a.num*exp(b.num);
    }

    friend ostream& operator<<(ostream &s, const LogSpace &a){
        if(a.iszero) s << "LogSpace: " << 0;
        else s << a.sign << "*e^"<< a.num;
        return s;
    }

    friend LogSpace operator-(const LogSpace &a) {
        LogSpace b(a);
        b.sign = -a.sign;
        return b;
    }

    friend bool operator==(const LogSpace&a, const LogSpace &b) {
        if(a.num == b.num && a.sign == b.sign && a.iszero == b.iszero ) return true;
        else return false;
    }

    friend double log(const LogSpace &a) {
        double b = a.num;
        return b;
    }

};

#endif /* LOGSPACE_H */
