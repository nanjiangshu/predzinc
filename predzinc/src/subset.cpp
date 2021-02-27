# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "subset.h"

void ksub_next ( int n, int k, int a[], bool *more )/*{{{*/

//******************************************************************************
//
//  Purpose:
//
//    KSUB_NEXT generates the subsets of size K from a set of size N.
//
//  Parameters:
//
//    Input, int N, the size of the set from which subsets are drawn.
//
//    Input, int K, the desired size of the subsets.  K must
//    be between 0 and N.
//
//    Output, int A[K].  A[I] is the I-th element of the
//    subset.  Thus A[I] will be an integer between 1 and N.
//    Note that the routine will return the values in A
//    in sorted order: 1 <= A[0] < A[1] < ... < A[K-1] <= N
//
//    Input/output, bool *MORE.  Set MORE = FALSE before first call
//    for a new sequence of subsets.  It then is set and remains
//    TRUE as long as the subset computed on this call is not the
//    final one.  When the final subset is computed, MORE is set to
//    FALSE as a signal that the computation is done.
//
{
  int j;
  static int m = 0;
  static int m2 = 0;
//
  if ( k < 0 || n < k )
  {
    cout << "\n";
    cout << "KSUB_NEXT - Fatal error!\n";
    cout << "N = " << n << "\n";
    cout << "K = " << k << "\n";
    cout << "but 0 <= K <= N is required!\n";
    exit ( 1 );
  }

  if ( !(*more) )
  {
    m2 = 0;
    m = k;
  }
  else
  {
    if ( m2 < n-m )
    {
      m = 0;
    }
    m = m + 1;
    m2 = a[k-m];
  }

  for ( j = 1; j <= m; j++ )
  {
    a[k+j-m-1] = m2 + j;
  }

  *more = ( a[0] != (n-k+1) );

  return;
}
/*}}}*/
void ksub_next2 ( int n, int k, int a[], int *in, int *iout )/*{{{*/

//*******************************************************************************
//
//  Purpose:
//
//    KSUB_NEXT2 generates the subsets of size K from a set of size N.
//
//  Discussion:
//
//    This routine uses the revolving door method.  It has no "memory".
//    It simply calculates the successor of the input set,
//    and will start from the beginning after the last set.
//
//  Modified:
//
//    29 May 2003
//
//  Reference:
//
//    Albert Nijenhuis and Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the size of the set from which subsets are drawn.
//    N must be positive.
//
//    Input, int K, the size of the desired subset.  K must be
//    between 0 and N.
//
//    Input/output, int A[K].  On input, the user must
//    supply a subset of size K in A.  That is, A must
//    contain K unique numbers, in order, between 1 and N.  On
//    output, A(I) is the I-th element of the output subset.
//    The output array is also in sorted order.
//
//    Output, int *IN, the element of the output subset which
//    was not in the input set.  Each new subset differs from the
//    last one by adding one element and deleting another.
//
//    Output, int *IOUT, the element of the input subset which
//    is not in the output subset.
//
{
  int j;
  int m;

  if ( n <= 0 )
  {
    cout << "\n";
    cout << "KSUB_NEXT2 - Fatal error!\n";
    cout << "  N = " << n << "\n";
    cout << "  but 0 < N is required!\n";
    exit ( 1 );
  }

  if ( k < 0 || n < k )
  {
    cout << "\n";
    cout << "KSUB_NEXT2 - Fatal error!\n";
    cout << "  N = " << n << "\n";
    cout << "  K = " << k << "\n";
    cout << "  but 0 <= K <= N is required!\n";
    exit ( 1 );
  }

  j = 0;

  for ( ; ; )
  {
    if ( 0 < j || ( k % 2 ) == 0 )
    {
      j = j + 1;

      if ( k < j )
      {
        a[k-1] = k;
        *in = k;
        *iout = n;
        return;
      }

      if ( a[j-1] != j )
      {
        *iout = a[j-1];
        *in = *iout - 1;
        a[j-1] = *in;

        if ( j != 1 )
        {
          *in = j - 1;
          a[j-2] = *in;
        }

        return;

      }

    }

    j = j + 1;
    m = n;

    if ( j < k )
    {
      m = a[j] - 1;
    }

    if ( m != a[j-1] )
    {
      break;
    }

  }

  *in = a[j-1] + 1;
  a[j-1] = *in;
  *iout = *in - 1;

  if ( j != 1 )
  {
    a[j-2] = *iout;
    *iout = j - 1;
  }

  return;
}
/*}}}*/
void ksub_next3 ( int n, int k, int a[], bool *more, int *in, int *iout )/*{{{*/

//*******************************************************************************
//
//  Purpose:
//
//    KSUB_NEXT3 generates the subsets of size K from a set of size N.
//
//  Discussion:
//
//    The routine uses the revolving door method.
//
//  Modified:
//
//    29 May 2003
//
//  Reference:
//
//    Albert Nijenhuis and Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the size of the set from which subsets are drawn.
//    N must be positive.
//
//    Input, int K, the size of the desired subsets.  K must be
//    between 0 and N.
//
//    Output, int A[K].  A(I) is the I-th element of the
//    output subset.  The elements of A are sorted.
//
//    Input/output, bool *MORE.  On first call, set MORE = FALSE
//    to signal the beginning.  MORE will be set to TRUE, and on
//    each call, the routine will return another K-subset.
//    Finally, when the last subset has been returned,
//    MORE will be set FALSE and you may stop calling.
//
//    Output, int *IN, the element of the output subset which
//    was not in the input set.  Each new subset differs from the
//    last one by adding one element and deleting another.  IN is not
//    defined the first time that the routine returns, and is
//    set to zero.
//
//    Output, int *IOUT, the element of the input subset which is
//    not in the output subset.  IOUT is not defined the first time
//    the routine returns, and is set to zero.
//
{
  int j;
  int m;

  if ( n <= 0 )
  {
    cout << "\n";
    cout << "KSUB_NEXT3 - Fatal error!\n";
    cout << "  N = " << n << "\n";
    cout << "  but 0 < N is required!\n";
    exit ( 1 );
  }

  if ( k < 0 || n < k )
  {
    cout << "\n";
    cout << "KSUB_NEXT3 - Fatal error!\n";
    cout << "  N = " << n << "\n";
    cout << "  K = " << k << "\n";
    cout << "  but 0 <= K <= N is required!\n";
    exit ( 1 );
  }

  if ( !(*more) )
  {
    *in = 0;
    *iout = 0;
    ivec_indicator ( k, a );
    *more = ( k != n );
    return;
  }

  j = 0;

  for ( ; ; )
  {
    if ( 0 < j || ( k % 2 ) == 0 )
    {
      j = j + 1;

      if ( a[j-1] != j )
      {
        *iout = a[j-1];
        *in = *iout - 1;
        a[j-1] = *in;

        if ( j != 1 )
        {
          *in = j - 1;
          a[j-2] = *in;
        }

        if ( k != 1 )
        {
          *more = ( a[k-2] == k-1 );
        }

        *more = ( !(*more) ) || ( a[k-1] != n );

        return;

      }

    }

    j = j + 1;
    m = n;

    if ( j < k )
    {
      m = a[j] - 1;
    }

    if ( m != a[j-1] )
    {
      break;
    }

  }

  *in = a[j-1] + 1;
  a[j-1] = *in;
  *iout = *in - 1;

  if ( j != 1 )
  {
    a[j-2] = *iout;
    *iout = j - 1;
  }

  if ( k != 1 )
  {
    *more = ( a[k-2] == k-1 );
  }

  *more = ( !(*more) ) || ( a[k-1] != n );

  return;
}
/*}}}*/
void ivec_indicator ( int n, int a[] )/*{{{*/

//******************************************************************************
//
//  Purpose:
//
//    IVEC_INDICATOR sets an integer vector to the indicator vector.
//
//  Modified:
//
//    25 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, int A[N], the initialized array.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = i + 1;
  }

}
/*}}}*/
void ksub_next4 ( int n, int k, int a[], bool *done )/*{{{*/

//*******************************************************************************
//
//  Purpose:
//
//    KSUB_NEXT4 generates the subsets of size K from a set of size N.
//
//  Discussion:
//
//    The subsets are generated one at a time.
//
//    The routine should be used by setting DONE to TRUE, and then calling
//    repeatedly.  Each call returns with DONE equal to FALSE, the array
//    A contains information defining a new subset.  When DONE returns
//    equal to TRUE, there are no more subsets.
//
//    There are ( N*(N-1)*...*(N+K-1)) / ( K*(K-1)*...*2*1) such subsets.
//
//  Modified:
//
//    03 June 2003
//
//  Parameters:
//
//    Input, int N, the size of the entire set.
//
//    Input, int K, the size of the desired subset.  K must be
//    between 0 and N.
//
//    Input/output, int A[K], contains information about
//    the subsets.  On the first call with DONE = TRUE, the input contents
//    of A don't matter.  Thereafter, the input value of A
//    should be the same as the output value of the previous call.
//    In other words, leave the array alone!
//    On output, as long as DONE is returned FALSE, A contains
//    information defining a subset of K elements of a set of N elements.
//    In other words, A will contain K distinct numbers (in order)
//    between 1 and N.
//
//    Input/output, bool *DONE.
//    On the first call, DONE is an input quantity with a value
//    of TRUE which tells the program to initialize data and
//    return the first subset.
//    On return, DONE is an output quantity that is TRUE as long as
//    the routine is returning another subset, and FALSE when
//    there are no more.
//
{
  int j;
  int jsave;

  if ( k < 0 || n < k )
  {
    cout << "\n";
    cout << "KSUB_NEXT4 - Fatal error!\n";
    cout << "  N = " << n << "\n";
    cout << "  K = " << k << "\n";
    cout << "  but 0 <= K <= N is required!\n";
    exit ( 1 );
  }
//
//  First call:
//
  if ( *done )
  {
    ivec_indicator ( k, a );

    if ( 0 < n )
    {
      *done = false;
    }
    else
    {
      *done = true;
    }
  }
//
//  Next call.
//
  else
  {
    if ( a[0] < n-k+1 )
    {
      *done = false;

      jsave = k-1;

      for ( j = 0; j < k-1; j++ )
      {
        if ( a[j] + 1 < a[j+1] )
        {
          jsave = j;
          break;
        }

      }

      ivec_indicator ( jsave, a );
      a[jsave] = a[jsave] + 1;
    }
    else
    {
      *done = true;
    }

  }

  return;
}
/*}}}*/
int i_min ( int i1, int i2 )/*{{{*/

//******************************************************************************
//
//  Purpose:
//
//    I_MIN returns the smaller of two integers.
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I_MIN, the smaller of i1 and i2.
//
{
  if ( i1 < i2 ) 
  {
    return i1;
  }
  else 
  {
    return i2;
  }

}
/*}}}*/
int i_max ( int i1, int i2 )/*{{{*/

//******************************************************************************
//
//  Purpose:
//
//    I_MAX returns the maximum of two integers.
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I_MAX, the larger of i1 and i2.
//
{
  if ( i2 < i1 ) 
  {
    return i1;
  }
  else 
  {
    return i2;
  }

}
/*}}}*/
int combin2 ( int n, int k )/*{{{*/

//******************************************************************************
//
//  Purpose:
//
//    COMBIN2 computes the binomial coefficient C(N,K).
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in integer arithmetic.
//
//  Formula:
//
//    ICNK = C(N,K) = N! / ( K! * (N-K)! )
//
//  Modified:
//
//    07 May 2003
//
//  Reference:
//
//    M L Wolfson and H V Wright,
//    Combinatorial of M Things Taken N at a Time,
//    ACM algorithm 160,
//    Communications of the ACM,
//    April, 1963.
//
//  Parameters:
//
//    Input, int N, K, are the values of N and K.
//
//    Output, int COMBIN2, the number of combinations of N
//    things taken K at a time.
//
{
  int cnk;
  int i;
  int mn;
  int mx;
//
  mn = i_min ( k, n-k );

  if ( mn < 0 )
  {
    cnk = 0;
  }
  else if ( mn == 0 )
  {
    cnk = 1;
  }
  else
  {
    mx = i_max ( k, n-k );
    cnk = mx + 1;

    for ( i = 2; i <= mn; i++ )
    {
      cnk = ( cnk * ( mx + i ) ) / i;
    }

  }

  return cnk;
}
/*}}}*/
