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

/**
 *  \brief Only free mem asc. to valid pointer, log warning otherwise
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

} // end namepspace stuff


#endif // end of stuff.hh
