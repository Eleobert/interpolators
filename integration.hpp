#include <cassert>

enum riennam_sum_type {left, right};

template<riennam_sum_type sum_type = left, typename Callable>
auto riemann_sum(double a, double b, int n, Callable& f)
{
    assert(n > 1);
    const auto dx  = (b - a) / (n - 1);
    const auto start = (sum_type == left)? a: a + dx;
    
    auto res = 0.0;
    for(auto x = start; x <= b; x += dx)
    {
        res += f(x);
    }

    return res * dx;
}


template<typename Callable>
auto trapz(double a, double b, int n, Callable& f)
{
    assert(n > 1);
    const auto dx  = (b - a) / (n - 1);
    
    auto res = 0.0;
    for(auto x = a; x <= b; x += dx)
    {
        res += (f(x) + f(x + dx));
    }

    return res * dx / 2;
}


template<typename Callable>
auto simpson(double a, double b, int n, Callable& f)
{
    assert(n > 1);
    const auto dx  = (b - a) / (n - 1);
    
    auto res = 0.0;
    for(auto x = a; x <= b; x += dx)
    {
        const auto l = x;      // lower
        const auto u = x + dx; // upper

        res += ((u - l) / 6.0) * (f(l) + 4 * f((l + u) / 2.0) + f(u)); 
    }

    return res;
}


