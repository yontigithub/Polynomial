#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "SparsePoly.h"
using namespace std;

int main() {
    SparsePoly p({{4,21.11},{10,50.3}});
    SparsePoly q({1,2,3.5});
    cout << "p(x) = " << p << '\n';
    cout << "q(x) = " << q << '\n';
    auto div = p/q;
    cout << "p(x) = q(x)*m(x) + r(x)" << '\n';
    cout << "m(x) = " << div.first << '\n';
    cout << "r(x) = " << div.second << '\n';

    cout << '\n';
    cout << q*div.first + div.second << "\n\n";

    ld x = 0.95;
    cout << "p(" << x << ") = " << p.eval(x) << '\n';
    cout << "q(" << x << ")*m(" << x << ") + r(" << x << ") = " <<
        q.eval(x) * div.first.eval(x) + div.second.eval(x) << '\n';

    return 0;
}
