#include "poly.h"
#include <iostream>
#include <map>
#include <stdexcept>
#include <algorithm>

static void clean(std::map<power, coeff, std::greater<power>> &terms)
{
    for (auto it = terms.begin(); it != terms.end(); )
    {
        if (it->second == 0)
            it = terms.erase(it);
        else
            ++it;
    }

    if (terms.empty())
        terms[0] = 0;  // zero polynomial
}

polynomial::polynomial()
{
    terms[0] = 0;
}

polynomial::polynomial(const polynomial &other)
{
    terms = other.terms;
}

polynomial &polynomial::operator=(const polynomial &other)
{
    if (this != &other)
        terms = other.terms;

    return *this;
}

// Templated iterator constructor
template <typename Iter>
polynomial::polynomial(Iter begin, Iter end)
{
    for (auto it = begin; it != end; ++it)
        terms[it->first] += it->second;

    clean(terms);
}

// Explicit instantiation so the linker can find it
template polynomial::polynomial(std::vector<std::pair<power, coeff>>::iterator,
                                std::vector<std::pair<power, coeff>>::iterator);

// ======================================================
// Print (debugging only)
// ======================================================
void polynomial::print() const
{
    for (auto &t : terms)
        std::cout << t.second << "x^" << t.first << " ";
    std::cout << std::endl;
}

// ======================================================
// Addition
// ======================================================
polynomial polynomial::operator+(const polynomial &other) const
{
    polynomial result(*this);

    for (auto &t : other.terms)
        result.terms[t.first] += t.second;

    clean(result.terms);
    return result;
}

polynomial polynomial::operator+(int x) const
{
    polynomial result(*this);

    result.terms[0] += x;
    clean(result.terms);

    return result;
}

// int + polynomial
polynomial operator+(int x, const polynomial &p)
{
    return p + x;
}

// ======================================================
// Multiplication
// ======================================================
polynomial polynomial::operator*(const polynomial &other) const
{
    polynomial result;
    result.terms.clear();

    for (auto &a : terms)
    {
        for (auto &b : other.terms)
        {
            power new_pow = a.first + b.first;
            coeff new_coef = a.second * b.second;

            result.terms[new_pow] += new_coef;
        }
    }

    clean(result.terms);
    return result;
}

polynomial polynomial::operator*(int x) const
{
    polynomial result(*this);

    for (auto &t : result.terms)
        t.second *= x;

    clean(result.terms);
    return result;
}

polynomial operator*(int x, const polynomial &p)
{
    return p * x;
}

// ======================================================
// Modulo (long division remainder)
// ======================================================
polynomial polynomial::operator%(const polynomial &divisor) const
{
    // divisor == 0?
    if (divisor.terms.size() == 1 && divisor.terms.begin()->second == 0)
        throw std::runtime_error("Divide by zero polynomial");

    polynomial remainder(*this);
    polynomial d(divisor);

    clean(remainder.terms);
    clean(d.terms);

    while (remainder.find_degree_of() >= d.find_degree_of() &&
           !(remainder.terms.size() == 1 &&
             remainder.terms.begin()->second == 0))
    {
        size_t deg_r = remainder.find_degree_of();
        size_t deg_d = d.find_degree_of();

        coeff coef_r = remainder.terms[deg_r];
        coeff coef_d = d.terms[deg_d];

        size_t pow_shift = deg_r - deg_d;
        coeff coef_shift = coef_r / coef_d;

        // temp = coef_shift * x^pow_shift
        polynomial temp;
        temp.terms.clear();
        temp.terms[pow_shift] = coef_shift;

        // subtract = temp * d
        polynomial subtract = temp * d;

        // remainder = remainder - subtract
        for (auto &t : subtract.terms)
            remainder.terms[t.first] -= t.second;

        clean(remainder.terms);
    }

    return remainder;
}

size_t polynomial::find_degree_of()
{
    return terms.begin()->first;
}

std::vector<std::pair<power, coeff>> polynomial::canonical_form() const
{
    std::vector<std::pair<power, coeff>> out;

    for (auto &t : terms)
        if (t.second != 0)
            out.push_back(t);

    if (out.empty())
        return {{0, 0}};

    return out;
}
