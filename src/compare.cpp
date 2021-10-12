#include <cmath>
#include <algorithm>
bool eq(double a, double b) {
    double diff = std::abs(a - b);

    if (diff <= 1e-12)
        return true;
 
    return (diff <= (std::max(std::abs(a), std::abs(b)) * 1e-8));
}

bool le(double a, double b) {
    if (a < b) return true;
    return eq(a,b);
}

bool ge(double a, double b) {
    if (a > b) return true;
    return eq(a,b);
}

