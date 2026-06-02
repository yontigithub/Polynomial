//
// Created by Yonatan Rappoport on 02/06/2026.
//

#ifndef SPARSE_POLY_UTILS_H
#define SPARSE_POLY_UTILS_H

#include <iostream>

typedef long double ld;
typedef long long ll;
typedef std::pair<ll, ld> pif;

static ld Pow(ld x, ll i) {
    if (i == 0) return 1;
    ld tmp = Pow(x,i/2);
    tmp *= tmp;

    if(i%2) tmp *= x;
    return tmp;
}

const ld zb = 0.000000001;

static bool is_zero(const ld &c, const ld zero = zb) {
    return abs(c) < zero;
}
static bool is_equal(const pif& a, const pif& b) {
    return is_zero(a.second-b.second) && (a.first == b.first);
}

#endif //SPARSE_POLY_UTILS_H
