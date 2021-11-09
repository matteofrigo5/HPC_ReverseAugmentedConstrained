#include <vector>
#include <string>
#include "pugixml.hpp"

bool stringCompare( std::string const & str1, std::string const & str2);

int get_int_value( pugi::xml_node const & node, std::string const & name, int & val );

int get_dbl_value( pugi::xml_node const & node, std::string const & name, double & val );

int get_lgl_value( pugi::xml_node const & node, std::string const & name, bool & val );

int get_string( pugi::xml_node const & node, std::string const & name, std::string & val );

void errorMsg( std::string const & message );

class warning
{
  public:

  warning()
  {};

  warning( std::string const & level )
  {
    push( level );
  };

  warning( warning const & warn )
  {
    m_levels = warn.m_levels;
  }

  ~warning() = default;

  void push( std::string const & level );

  void print();

  void print( std::string const & level );

  private:

  std::vector< std::string > m_levels;
};
