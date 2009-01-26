/**
 *  \file stuff.hh
 *  \brief  contains some stuff
 **/
#ifndef STUFF_HH_INCLUDED
#define STUFF_HH_INCLUDED

#define SEGFAULT int*i=0;*i=9;

#include <iomanip>


namespace Stuff
{

/**
 *  \todo doc me
 **/
template < class ReturnType >
ReturnType fromString(const std::string& s)
{
    std::stringstream ss;
    ss << s;
    ReturnType r;
    ss >> r;
    return r;
}

/**
 *  \todo doc
 **/
template < class ReturnType >
std::string toString(const ReturnType& s)
{
    std::stringstream ss;
    ss << s;
    std::string r;
    ss >> r;
    return r;
}

struct TexOutput
{
    void printTexInfo(std::ostream&) const
    {

    }
};

/**
 *  \brief Only free mem pointed to by valid pointer, log warning otherwise
 *
 **/

template < class T >
void safe_delete ( T t )
{
    if (t){
        delete t;
        t=0;
    }
    //else log warning
}

/**
 *  \brief prints a Dune::FieldVector
 *
 *  or anything compatible in terms of Iterators
 *  \tparam T
 *          should be Dune::FieldVector or compatible
 *  \tparam out
 *          std::ostream or compatible
 *  \param  arg
 *          Vector to be printed
 *  \param  name
 *          name to be printed along
 **/
template < class T, class stream >
void printFieldVector( T& arg, std::string name, stream& out, std::string prefix = "" )
{
    out << "\n" << prefix << "printing " << name << " (Dune::FieldVector)" << std::endl;
    typedef typename T::ConstIterator
        IteratorType;
    IteratorType itEnd = arg.end();
    for ( IteratorType it = arg.begin(); it != itEnd; ++it ) {
            out << prefix << std::setw( 7 ) << std::setprecision( 3 ) << *it;
    }
}

/**
 *  \brief prints a Dune::FieldMatrix
 *
 *  or anything compatible in terms of Iterators
 *  \tparam T
 *          should be Dune::FieldMatrix or compatible
 *  \tparam out
 *          std::ostream or compatible
 *  \param  arg
 *          Matrix to be printed
 *  \param  name
 *          name to be printed along
 **/
template < class T, class stream >
void printFieldMatrix( T& arg, std::string name, stream& out, std::string prefix = "" )
{
    out << "\n" << prefix << "printing " << name << " (Dune::FieldMatrix)";
    typedef typename T::ConstRowIterator
        RowIteratorType;
    typedef typename T::row_type::ConstIterator
        VectorInRowIteratorType;
    unsigned int row = 1;
    RowIteratorType rItEnd = arg.end();
    for ( RowIteratorType rIt = arg.begin(); rIt != rItEnd; ++rIt ) {
        out << "\n" << prefix << "row " << row << ":";
        VectorInRowIteratorType vItEnd = rIt->end();
        for (   VectorInRowIteratorType vIt = rIt->begin(); vIt != vItEnd; ++vIt ) {
            out << prefix << std::setw( 7 ) << std::setprecision( 3 ) << *vIt;
        }
        row += 1;
    }
}

template < class T, class stream >
void printSparseRowMatrixMatlabStyle( const T& arg, const std::string name, stream& out )
{
    out << "\n" << name << " = [ ";
    for ( int row = 0; row < arg.rows(); row++ ) {
        for ( int col = 0; col < arg.cols(); col++ ) {
            out << std::setw( 8 ) << std::setprecision( 2 ) << arg(row,col);
        }
        out << ";" << std::endl;
    }
    out << "];" << std::endl;
}

template < class T, class stream >
void printDiscreteFunctionMatlabStyle( const T& arg, const std::string name, stream& out )
{
    out << "\n" << name << " = [ ";
    typedef typename T::ConstDofIteratorType
        ConstDofIteratorType;
    ConstDofIteratorType itEnd = arg.dend();
    for ( ConstDofIteratorType it = arg.dbegin(); it != itEnd; ++it ) {
        out << std::setprecision( 2 ) << *it;
        out << ";" << std::endl;
    }
    out << "];" << std::endl;
}

template < class Matrix, class Function >
void DiagonalMult( const Matrix& matrix, Function& f )
{
    Function diag( "temp", f.space() );
    matrix.getDiag( diag );

    typedef typename Function::DofIteratorType DofIteratorType;
    DofIteratorType diag_it = diag.dbegin();
    DofIteratorType f_it = f.dbegin();

    for(int row=0; row< matrix.size(0); row++)
    {
        (*f_it) *= (*diag_it);
        ++f_it;
        ++diag_it;
    }
    return;
}

template < class Stream, class Type >
void printDoubleVec( Stream& stream, const Type * vec, const unsigned int N )
{
    stream << "\n [ " << std::setw(5);
    for ( unsigned int i = 0; i < N ; ++i )
        stream << vec[i] << " ";

    stream << " ] " << std::endl;
}

template < class Stream, class DiscFunc >
void oneLinePrint( Stream& stream, const DiscFunc& func )
{
    typedef typename DiscFunc::ConstDofIteratorType DofIteratorType;
    DofIteratorType it = func.dbegin();
    stream << "\n" << func.name() << ": [ ";
    for ( ; it != func.dend(); ++it )
        stream << std::setw(5) << *it << "  ";

    stream << " ] " << std::endl;
}

template < class Function >
void addScalarToFunc( Function& f, double sc )
{
    typedef typename Function::DofIteratorType DofIteratorType;
    DofIteratorType it = f.dbegin();
    for ( ; it != f.dend(); ++it )
        *it += sc;
    return;
}

} // end namepspace stuff


#endif // end of stuff.hh
