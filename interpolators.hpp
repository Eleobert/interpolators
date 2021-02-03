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
    // make the last value match b exactly
    res[n - 1] = b;
    return res;
}


template<typename Container>
class poly
{
    using type = typename Container::value_type;
    std::vector<type> xs;
    std::vector<type> data;

    auto& qs(int i, int j){return data[i * xs.size() + j];}
public:

    poly(const Container& xs, const Container& ys): xs(xs), data(xs.size() * xs.size())
    {
        assert(xs.size() == ys.size());
        std::copy(std::begin(xs), std::end(xs), this->xs.begin());
        std::copy(std::begin(ys), std::end(ys), data.begin());
    }
    double operator()(double x);
};


template<typename Container>
double poly<Container>::operator()(double x)
{
    for(size_t i = 1; i < xs.size(); i++)
    {
        for(size_t j = i; j < xs.size(); j++)
        {
            qs(i, j) = ((x - xs[j-i]) * qs(i-1, j) - (x - xs[j]) * qs(i-1, j-1)) / (xs[j] - xs[j-i]);
        }
    }
    return qs(xs.size() - 1, xs.size() - 1);
}


template<typename RandomAcessContainer>
class linear
{
    using type = typename RandomAcessContainer::value_type;
    std::vector<type> xs;
    std::vector<type> ys;
    const bool initialized = false;

public:
    linear(const RandomAcessContainer& xs, const RandomAcessContainer& ys);
    linear() = default;

    double operator()(double x);
};


template<typename Container>
linear<Container>::linear(const Container& xs, const Container& ys):
                                    xs(xs.size()), ys(ys.size()), initialized(true)
{
    assert(xs.size() == ys.size());
    std::copy(std::begin(xs), std::end(xs), this->xs.begin());
    std::copy(std::begin(ys), std::end(ys), this->ys.begin());
}


template<typename Container>
double linear<Container>::operator()(double x)
{
    assert(initialized == true);

    auto x_high = std::upper_bound(xs.begin(), xs.end(), x);
    assert(x_high != xs.end());
    auto i = x_high - xs.begin() - 1u; //index of x_low
    return ys[i] + (x - xs[i]) * (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]);
}