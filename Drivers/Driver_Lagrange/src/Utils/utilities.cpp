#include <iostream>
#include <mpi.h>

#include "utilities.hpp"

#define PRINTVALS true

namespace
{

bool compareChar( char const & c1, char const & c2 )
{
  if( c1 == c2 )
  {
    return true;
  }
  else if( std::toupper( c1 ) == std::toupper( c2 ) )
  {
    return true;
  }
  return false;
}

void mpi_print( char const * const str )
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if( rank == 0 )
  {
    std::cout << std::string( str ) << std::endl;
  }
}

};

bool stringCompare( std::string const & str1, std::string const & str2)
{
  return ( ( str1.size() == str2.size() ) && std::equal( str1.begin(), str1.end(), str2.begin(), &compareChar ) );
}

int get_int_value( pugi::xml_node const & node, std::string const & name, int & val )
{
  pugi::xml_attribute const attr = node.attribute( name.c_str() );
  if( attr != NULL )
  {
    val = attr.as_int();
#if PRINTVALS
    char buf[200];
    std::snprintf( buf, 200, "%20s %15i", name.c_str(), val );
    mpi_print( buf );
#endif
    return 0;
  }
  else
  {
    return -1;
  }
};

int get_dbl_value( pugi::xml_node const & node, std::string const & name, double & val )
{
  pugi::xml_attribute const attr = node.attribute( name.c_str() );
  if( attr != NULL )
  {
    val = attr.as_double();
#if PRINTVALS
    char buf[200];
    std::snprintf( buf, 200, "%20s %15.6e", name.c_str(), val );
    mpi_print( buf );
#endif
    return 0;
  }
  else
  {
    return -1;
  }
};

int get_lgl_value( pugi::xml_node const & node, std::string const & name, bool & val )
{
  pugi::xml_attribute const attr = node.attribute( name.c_str() );
  if( attr != NULL )
  {
    val = attr.as_bool();
#if PRINTVALS
    char buf[200];
    std::snprintf( buf, 200, "%20s %15i", name.c_str(), val );
    mpi_print( buf );
#endif
    return 0;
  }
  else
  {
    return -1;
  }
};

int get_string( pugi::xml_node const & node, std::string const & name, std::string & val )
{
  pugi::xml_attribute const attr = node.attribute( name.c_str() );
  if( attr != NULL )
  {
    val = attr.as_string();
#if PRINTVALS
    char buf[200];
    std::snprintf( buf, 200, "%20s %50s", name.c_str(), val.c_str() );
    mpi_print( buf );
#endif
    return 0;
  }
  else
  {
    return -1;
  }
};

void warning::push( std::string const & level )
{
  m_levels.push_back( level );
}

void warning::print()
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if( rank == 0 )
  {
    std::cout << "Warning: ";
    for( size_t i = 0; i < m_levels.size()-1; ++i )
    {
      std::cout << m_levels[i] << ":";
    }
    std::cout << m_levels[m_levels.size()-1];
    std::cout << " not provided. Default used.\n";
  }
}

void warning::print( std::string const & level )
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if( rank == 0 )
  {
    std::cout << "Warning: ";
    for( size_t i = 0; i < m_levels.size(); ++i )
    {
      std::cout << m_levels[i] << ":";
    }
    std::cout << level;
    std::cout << " not provided. Default used.\n";
  }
}

void errorMsg( std::string const & message )
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if( rank == 0 )
  {
    std::cout << message << std::endl;
  }
}
