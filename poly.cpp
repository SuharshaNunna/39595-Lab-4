#include "poly.h"
#include <iostream>
#include <map>
#include <stdexcept>
#include <algorithm>
#include <pthread.h>
#include <unordered_map>

static void clean(std::map<power, coeff, std::greater<power>> &terms)
{
    for (auto it = terms.begin(); it != terms.end();)
    {
        if (it->second == 0)
        {
            it = terms.erase(it);
        }
        else
        {
            it++;
        }
    }

    if (terms.empty())
    {
        terms[0] = 0;  // zero polynomial
    }
}

// parallel multiplication helpers

struct multiplication
{
    const std::vector<std::pair<power, coeff>> *a;
    const std::vector<std::pair<power, coeff>> *b;

    size_t start;
    size_t end;

    std::unordered_map<power, coeff> *partial;
};

static void *multiply(void *arg)
{
    multiplication *task = static_cast<multiplication*>(arg);
    const auto &a = *(task->a);
    const auto &b = *(task->b);
    auto &x = *(task->partial);

    for (size_t i = task->start; i < task->end; i++)
    {
        const auto &at = a[i];
        for (const auto &bt : b)
        {
            power p = at.first + bt.first;
            coeff c = at.second * bt.second;
            x[p] += c;
        }
    }

    return nullptr;
}

// polynomial member functions

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
    {
        terms = other.terms;
    }
    return *this;
}

template <typename Iter>
polynomial::polynomial(Iter begin, Iter end)
{
    for (auto it = begin; it != end; it++)
    {
        terms[it->first] += it->second;
    }
    clean(terms);
}

template polynomial::polynomial(std::vector<std::pair<power, coeff>>::iterator, std::vector<std::pair<power, coeff>>::iterator);

polynomial polynomial::operator+(const polynomial &other) const
{
    polynomial result(*this);

    for (auto &t : other.terms)
    {
        result.terms[t.first] += t.second;
    }
    
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

polynomial operator+(int x, const polynomial &p)
{
    return p + x;
}

// parallel operator* implementation using unordered_map

polynomial polynomial::operator*(const polynomial &other) const
{
    // zero checks
    if ((terms.size() == 1 && terms.begin()->second == 0) || (other.terms.size() == 1 && other.terms.begin()->second == 0))
    {
        return polynomial(); 
    }

    std::vector<std::pair<power, coeff>> a(terms.begin(), terms.end());
    std::vector<std::pair<power, coeff>> b(other.terms.begin(), other.terms.end());

    if (a.empty() || b.empty())
    {
        return polynomial();
    }

    const int MAX_THREADS = 8;
    int num = static_cast<int>(std::min<size_t>(MAX_THREADS, a.size()));

    // small polynomials
    if (num <= 1)
    {
        polynomial result;
        result.terms.clear();

        for (const auto &at : a)
        {
            for (const auto &bt : b)
            {
                power p = at.first + bt.first;
                coeff c = at.second * bt.second;
                result.terms[p] += c;
            }
        }

        clean(result.terms);
        return result;
    }

    std::vector<pthread_t> threads(num);
    std::vector<multiplication> tasks(num);
    std::vector<std::unordered_map<power, coeff>> partials(num);

    size_t n = a.size();
    size_t x = (n + num - 1) / num; 

    // launch threads
    for (int t = 0; t < num; ++t)
    {
        size_t start = t * x;
        size_t end = std::min(n, start + x);

        if (start >= end)
        {
            threads[t] = 0;
            continue;
        }

        tasks[t].a = &a;
        tasks[t].b = &b;
        tasks[t].start = start;
        tasks[t].end = end;
        tasks[t].partial = &partials[t];

        pthread_create(&threads[t], nullptr, multiply, &tasks[t]);
    }

    // join threads
    for (int t = 0; t < num; ++t)
    {
        if (threads[t] != 0)
        {
            pthread_join(threads[t], nullptr);
        }
    }

    // combine partial results 
    polynomial result;
    result.terms.clear();

    for (const auto &local : partials)
    {
        for (const auto &kv : local)
        {
            result.terms[kv.first] += kv.second;
        }
    }

    clean(result.terms);
    return result;
}

polynomial polynomial::operator*(int x) const
{
    polynomial result(*this);

    for (auto &t : result.terms)
    {
        t.second *= x;
    }

    clean(result.terms);
    return result;
}

polynomial operator*(int x, const polynomial &p)
{
    return p * x;
}

polynomial polynomial::operator%(const polynomial &mod) const
{

    if (mod.terms.size() == 1 && mod.terms.begin()->second == 0)
    {
        throw std::runtime_error("error");
    }

    polynomial remainder(*this);
    polynomial d(mod);

    clean(remainder.terms);
    clean(d.terms);

    while (remainder.find_degree_of() >= d.find_degree_of() && !(remainder.terms.size() == 1 && remainder.terms.begin()->second == 0))
    {
        size_t deg_r = remainder.find_degree_of();
        size_t deg_d = d.find_degree_of();

        coeff coef_r = remainder.terms[deg_r];
        coeff coef_d = d.terms[deg_d];

        size_t pow = deg_r - deg_d;
        coeff coef = coef_r / coef_d;
        polynomial temp;
        
        temp.terms.clear();
        temp.terms[pow] = coef;

        polynomial subtract = temp * d;

        for (auto &t : subtract.terms)
        {
            remainder.terms[t.first] -= t.second;
        }

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
    {
        if (t.second != 0)
        {
            out.push_back(t);
        }
    }

    if (out.empty())
    {
        return {{0, 0}};
    }

    return out;
}
