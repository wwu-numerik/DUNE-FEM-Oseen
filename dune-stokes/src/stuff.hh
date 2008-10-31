/**
 *  \file stuff.hh
 *  \brief  contains some stuff
 **/
#ifndef STUFF_HH_INCLUDED
#define STUFF_HH_INCLUDED

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
 *  \todo doc me
 **/

template < class T >
void safe_delete ( T t )
{
    if (t){
        delete t;
        t=0;
    }
}

} //end namepsace stuff


#endif
