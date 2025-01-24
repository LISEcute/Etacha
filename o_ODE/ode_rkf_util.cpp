//# include <stdlib.h>
//# include <iostream.h>
//# include <fstream.h>
//# include <iomanip.h>
#include <QDebug>
# include <math.h>
# include <time.h>


#include "o_ODE/ode.hpp"

//****************************************************************************80

double r8_epsilon ( )
{
  const double value = 2.220446049250313E-016;
  return value;
}
//****************************************************************************80
double r8_max ( double x, double y )
{
  double value;

 if ( y < x )  {    value = x;  }
 else  {    value = y;  }
 return value;
}
//****************************************************************************80

double r8_min ( double x, double y )
{  double value;
  if ( y < x )  {    value = y;  }
  else          {    value = x;  }
  return value;
}
//****************************************************************************80
double r8_sign ( double x )
{
 double value;
  if ( x < 0.0 )  {    value = -1.0;  }
  else  {    value = 1.0;  }
  return value;
}
//****************************************************************************80

double r8_abs ( double x )
{
  if ( 0.0 <= x )  {    return x;  }
  else             {    return ( -x );
  }
}
//****************************************************************************80
float r4_sign ( float x )
{
  if ( x < 0.0 )  {    return ( -1.0 );  }
  else  {    return ( +1.0 );  }
}
//****************************************************************************80
//****************************************************************************80
void timestamp ( )
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm_ptr;
//  size_t len;
  time_t now;

  now = time ( NULL );
  tm_ptr = localtime ( &now );

//  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );
  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

//  qDebug() << time_buffer << "\n";
  return;
# undef TIME_SIZE
}
