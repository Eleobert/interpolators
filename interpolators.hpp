#pragma once

#include <algorithm>
#include <vector>
#include <cassert>
#include <functional>
#include <numeric>
#include <cmath>


template <typename Container>
auto linspace(typename Container::value_type a, typename Container::value_type b, size_t n)
{
    assert(b > a);
    assert(n > 1);

    Container res(n);
    const auto step = (b - a) / (n - 1);
    auto val = a;
    for(auto& e: res)
    {
        e = val;
        val += step;
    }
    return res;
}


template<typename RandomAcessContainer>
class lagrange
{
    std::vector<double> dem;
    RandomAcessContainer x;
    RandomAcessContainer y;
    const bool initialized = false;

    // function for polinomial i
    double p(double xi, size_t i);

public:
    lagrange(const RandomAcessContainer& x, const RandomAcessContainer& y);
    lagrange() = default;

    double operator()(double x);
};

template<typename RandomAcessContainer>
lagrange<RandomAcessContainer>::lagrange(const RandomAcessContainer& x, const RandomAcessContainer& y): x(x), y(y), 
                                    initialized(true)
{
    assert(x.size() == y.size());
    dem.resize(x.size());

    // precompute the polinomials denominators
    for(int i = 0; i < x.size(); i++)
    {
        dem[i] = 1;
        for(int j = 0; j < x.size(); j++)
        {
            if(i != j)
            {
                dem[i] *= x[i] - x[j];
            }
        }
    }
}


template<typename RandomAcessContainer>
double lagrange<RandomAcessContainer>::p(double xi, size_t i)
{
    auto num = 1.0;
    for(size_t j = 0; j < x.size(); j++)
    {
        if(i != j)
        {
            num *= xi - x[j];
        }
    }
    return num * y[i] / dem[i];   
}

template<typename RandomAcessContainer>
double lagrange<RandomAcessContainer>::operator()(double xi)
{
    assert(initialized == true);

    auto res = 0.0;
    for(int i = 0; i < x.size(); i++)
    {
        res += p(xi, i);
    }
    return res;
}


template<typename RandomAcessContainer>
class newton
{
    std::vector<std::vector<double>> divided_differences;
    double x0, y0, step;
    const bool initialized = false;

public:

    newton(double x0, double step, const RandomAcessContainer& ys);
    newton() = default;

    double operator()(double x);
};


template<typename RandomAcessContainer>
newton<RandomAcessContainer>::newton(double x0, double step, const RandomAcessContainer& ys): 
                                            step(step), x0(x0), y0(ys[0]), divided_differences(ys.size()),
                                            initialized(true)
{
    divided_differences[0].resize(ys.size());
    std::copy(ys.begin(), ys.end(), divided_differences[0].begin());
    
    for(size_t i = 1; i < ys.size(); i++)
    {
        auto& pdd = divided_differences[i - 1]; // previous divided differences
        auto& cdd = divided_differences[i];     // current  divided differences
        
        cdd.resize(pdd.size() - 1);

        std::transform(std::next(pdd.begin()), pdd.end(), pdd.begin(), cdd.begin(), std::minus{});
    }
}


template<typename RandomAcessContainer>
double newton<RandomAcessContainer>::operator()(double x)
{
    assert(initialized == true);

    auto res  = y0;
    auto q    = (x - x0) / step;
    auto poly = q;

    for(int i = 1; i < divided_differences.size(); i++)
    {
        res  += poly * divided_differences[i][0];
        poly *= (q - i) / (i + 1);
    }
    return res;
}