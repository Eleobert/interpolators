#include <algorithm>
#include <vector>
#include <cassert>

#include <numeric>



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
    const RandomAcessContainer x;
    const RandomAcessContainer y;

    // function for polinomial i
    double p(double xi, double i);

public:
    lagrange(const RandomAcessContainer& x, const RandomAcessContainer& y): x(x), y(y)
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

    double operator()(double x);
};

template<typename RandomAcessContainer>
double lagrange<RandomAcessContainer>::p(double xi, double i)
{
    auto num = 1.0;
    for(int j = 0; j < x.size(); j++)
    {
        if( i != j)
        {
            num *= xi - x[j];
        }
    }
    return num * y[i] / dem[i];   
}

template<typename RandomAcessContainer>
double lagrange<RandomAcessContainer>::operator()(double xi)
{
    auto res = 0.0;
    for(int i = 0; i < x.size(); i++)
    {
        res += p(xi, i);
    }
    return res;
}