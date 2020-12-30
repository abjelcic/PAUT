#ifndef UTILITY_H
#define UTILITY_H

#include <vector>
#include <initializer_list>
#include <fstream>
#include <cmath>
#include <exception>



template <typename T>
std::vector<T> uniformGrid1D( T start , T increment , T end , T tol )
{
    std::vector<T> ans;

    for( auto x = start ; x <= end+tol ; x += increment )
        ans.push_back(x);

    return ans;
}






template <typename F, int N>
class Vec
{
private:
    F m_coords[N];

public:

    static F dot( Vec<F,N> const & v1 , Vec<F,N> const & v2 )
    {
        F ans {0};
        for(unsigned i=0; i<N; ++i)
            ans += v1.m_coords[i] * v2.m_coords[i];
        return ans;
    }

    Vec()
    {
        for(unsigned i=0; i<N; ++i)
            m_coords[i] = 0.0;
        return;
    }

    Vec( std::initializer_list<F> const & v ) 
    {
        if( v.size() != N )
            throw std::logic_error("Vector and initializer list dimensions not compatible.");

        unsigned int i = 0;
        for( auto const & x : v )
            m_coords[i++] = x;

        return;
    }

    Vec( Vec<F,N> const & v )
    {
        for( unsigned i=0 ; i<N ; ++i )
            m_coords[i] = v.m_coords[i];

        return;
    }

    Vec<F,N> & operator += ( Vec<F,N> const & v )
    {
        for( unsigned i=0 ; i<N ; ++i )
            m_coords[i] += v.m_coords[i];

        return *this;
    }

    Vec<F,N> & operator -= ( Vec<F,N> const & v )
    {
        for( unsigned i=0 ; i<N ; ++i )
            m_coords[i] -= v.m_coords[i];

        return *this;
    }

    F norm() const
    {
        F ans {0};
        for(unsigned i=0; i<N; ++i)
            ans += m_coords[i]*m_coords[i];
        return std::sqrt(ans);    
    } 

    F operator () (unsigned i) const
    {
        #ifdef DEBUG
        if( i<=0 || i>N )
            throw std::out_of_range("Vector coordinate out of bounds.");
        #endif

        return m_coords[i-1];
    }

    F & operator () (unsigned i)
    {
        #ifdef DEBUG
        if( i<=0 || i>N )
            throw std::out_of_range("Vector coordinate out of bounds.");
        #endif

        return m_coords[i-1];
    }

    template <typename FF>
    Vec<F,N> & operator *= ( FF a )
    {
        for( unsigned i=0 ; i<N ; ++i )
            m_coords[i] *= static_cast<F>(a);
        return *this;
    }

    template <typename FF, int NN>
    friend std::ostream & operator << ( std::ostream & out , Vec<FF,NN> const & v );
};

template <typename F1, typename F2, int N>
Vec<F2,N> operator * ( F1 const & a , Vec<F2,N> const & v )
{
    auto v1 = v;
    v1 *= static_cast<F2>(a);
    return v1;
}

template <typename F1, typename F2, int N>
Vec<F2,N> operator * ( Vec<F1,N> const & v , F2 const & a )
{
    auto v1 = v;
    v1 *= static_cast<F1>(a);
    return v1;
}

template <typename F, int N>
Vec<F,N> operator + ( Vec<F,N> const & v1 , Vec<F,N> const & v2 )
{
    auto v = v1;
    v += v2;
    return v;
}

template <typename F, int N>
Vec<F,N> operator - ( Vec<F,N> const & v1 , Vec<F,N> const & v2 )
{
    auto v = v1;
    v -= v2;
    return v;
}

template <typename F, int N>
std::ostream & operator << ( std::ostream & out , Vec<F,N> const & v )
{
    out << "[";
    for(unsigned i=0; i<N; ++i)
        out << v.m_coords[i] << ( i<N-1 ? "," : "" );
    out << "]";
    return out;
}







#endif // UTILITY_H
