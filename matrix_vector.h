/*

           C++ Library for MATRIX and VECTOR operations
 
   developed by : G R Krishna Chand Avatar, MTech (Aerospace Engineering)

*/

#include<iostream>
#include<cmath>
#include<vector>
#include<cassert>
#include<algorithm>

using namespace std;


// INDICES START FROM 1
// Class definition for a 2D matrix using vector template class
class vec{
  private:
    vector<double> v;
    int size;
  public:
    // Constructor
    vec(int n){
        size = n;
        for(int i=1;i<=size;i++)
            v.push_back(0.0);
    } 
    
    // Operator overloading

    double & operator()(size_t j){ // "&" indicates reference to the current functor
     // so that this could be used as l-value
        return v[j-1];
    } 

    vec operator+(vec v1){
        vec temp(size);
        // Check whether same dimensions
        for(int j=1; j<=size; j++)
           temp(j) = v[j-1] + v1(j); 
        return temp;
    }

    vec operator-(vec v1){
        vec temp(size);
        // Check whether same dimensions
        for(int j=1; j<=size; j++)
           temp(j) = v[j-1] - v1(j); 
        return temp;
    }

    // Vector double operation
    vec operator*(double d){
        vec temp(size);
        for(int j=1; j<=size; j++)
            temp(j) = v[j-1]*d;
        return temp;
    }

    vec operator/(double d){
        vec temp(size);
        for(int j=1; j<=size; j++)
            temp(j) = v[j-1]/d;
        return temp;
    }
    
    // VECTOR FUNCTIONS 

    // Return size
    int length(){
        return size;
    }

    // Norm
    double norm(){
      double sum = 0.0;
      for(int j=0;j<size;j++)
          sum += v[j]*v[j];
      return pow(sum,0.5);
    }
    
    // Dot product
    double dot(vec v1){
      double sum = 0.0;
      for(int j=1;j<=size;j++)
          sum += v[j-1]*v1(j);
      return sum;
    }

    friend double dot(vec v1, vec v2);

    // Normalizing
    vec normalize(){
        vec temp(length());
        double v_norm = norm();
        for(int j=1;j<=size;j++)
           temp(j) = v[j-1]/v_norm;
        return temp;
    }
 
    void display(){

        cout<<"Current vector: "<<endl;
        for(vector<double>::iterator j = v.begin(); j < v.end();j++)
          cout<<*j<<"\t"; // using iterator to scan through the elements
          cout<<endl;
    }
    
    void roundoffzeros(){
        for(int i=0;i<length();i++)
          if(fabs(v[i])<1.0e-15)
            v[i] = 0.0;
    }

    // maximum or minimum

    double maximum(){

      double max = *max_element(v.begin(), v.end());
      return max;
    }

    // Destructor
    ~vec(){}
};

// Class definition for a 2D matrix using vector template class
class matrix2D{
  private: 
    vector< vector<double> >  matrix;
  public:
    // Defining constructor using initialisation list
    matrix2D(size_t rows, size_t cols):matrix(rows, vector<double>(cols)){}

    // Using operator overloading to directly access (i,j)-th element of the matrix
    double & operator()(size_t i, size_t j){ // "&" indicates reference to the current functor
     // so that this could be used as l-value
        return matrix[i-1][j-1];
    }

    // Determining number of rows and cols

    size_t nrows() { return matrix.size();}
    size_t ncols() { return matrix[0].size();}
    
    void display(){

      cout<<"Current matrix: "<<endl;
      for(int i=1;i<=nrows();i++){
        for(int j=1;j<=ncols();j++)
          cout<< matrix[i-1][j-1] <<"\t";
          cout<< endl << endl;
      }       
    }

    // MATRIX FUNCTIONS

    // maximum or minimum

    double maximum(){

      double max;

      max = *max_element(matrix[0].begin(), matrix[0].end());

      for(int i=1; i< nrows(); i++ )
      {
         if (max < *max_element(matrix[i].begin(), matrix[i].end()))
            max = *max_element(matrix[i].begin(), matrix[i].end());   
      }

      return max;
    }

    double minimum(){

      double min;

      min = *min_element(matrix[0].begin(), matrix[0].end());

      for(int i=1; i< nrows(); i++ )
      {
         if (min > *min_element(matrix[i].begin(), matrix[i].end()))
            min = *min_element(matrix[i].begin(), matrix[i].end());   
      }

      return min;
    }

    // transpose of the matrix

    matrix2D transpose(){
       
      matrix2D temp(ncols(), nrows());

      for(int i=1; i<= nrows(); i++)
        for(int j=1; j<= ncols(); j++)
          temp(j,i) = matrix[i-1][j-1];

      return temp;

    }

    // OPERATOR OVERLOADING 

    matrix2D operator+(matrix2D B){
        matrix2D temp(nrows(), ncols());
        int check = (nrows()==B.nrows() && ncols()==B.ncols())? 1 : 0;
        if(check){
            for(int i=1; i<=nrows(); i++)
                for(int j=1; j<=ncols(); j++)
                    temp(i,j) = matrix[i-1][j-1] + B(i,j);
        }
        else
            cout<<"Matrix dimensions are inconsistent"<< endl;
        return temp;
    }

    matrix2D operator-(matrix2D B){
        matrix2D temp(nrows(), ncols());
        int check = (nrows()==B.nrows() && ncols()==B.ncols())? 1 : 0;
        if(check){
            for(int i=1; i<=nrows(); i++)
                for(int j=1; j<=ncols(); j++)
                    temp(i,j) = matrix[i-1][j-1] - B(i,j);
        }

        else
            cout<<"Matrix dimensions are inconsistent"<< endl;
        return temp;
    }

    matrix2D operator*(matrix2D B){
        matrix2D temp(nrows(), B.ncols());
        int check = (ncols()==B.nrows())? 1 : 0;
        if(check){
            for(int i=1; i<=nrows(); i++)
                for(int j=1; j<=B.ncols(); j++)
                    for(int k=1; k<=ncols(); k++)
                       temp(i,j) += matrix[i-1][k-1]*B(k,j);
        }

        else
            cout<<"Matrix dimensions are inconsistent"<< endl;
        return temp;
    }
   
    // Matrix vector multiplication

    vec operator*(vec y){
        vec x(y.length());
        if(nrows()==ncols() && y.length() == nrows())
        {
            for(int i = 1; i<=nrows(); i++)
                for(int j = 1; j<=nrows(); j++)
                    x(i) += matrix[i-1][j-1]*y(j);
        }
        else
            cout<< " Matrix and vector dimensions inconsistent"<<endl;
        
        x.roundoffzeros();  
        return x;
    }

    void roundoffzeros(){
        for(int i=0;i<nrows();i++)
          for(int j=0; j<ncols();j++)
            if(fabs(matrix[i][j])<1.0e-15)
               matrix[i][j] = 0.0;
    }

    // Destructor
    ~matrix2D(){}
};

double dot(vec v1, vec v2){

    double product = 0.0;
    int size = v1.length();

    if(v1.length() == v2.length())  
        for(int i=1; i<=size; i++)
           product += v1(i)*v2(i);
    else 
    {
        cout<<"Dimensions not consistent"<<endl;
        return -1;
    }

    return product;
}

//********************* Identity matrix ******************

matrix2D eye(int n){

    matrix2D identity(n,n);
    for(int i=1; i<=n; i++)
        identity(i,i) = 1.0;
    return identity;
}

matrix2D resize(matrix2D A, int m, int n){
    
    matrix2D resized(m,n);
    int p,q; // matrix size indicators

    if(m <= A.nrows())
        p = m;
    else
        p = A.nrows();

    if(n <= A.ncols())
        q = n;
    else
        q = A.ncols();

    for(int i=1; i<=p; i++)
        for(int j=1; j<=q; j++ )
            resized(i,j) = A(i,j);
   
    return resized; 
}

// *********  Solving Ax = b using Gauss elimination procedure *******

void swap_rows(matrix2D& Ab, int n, int pivot_row, int j){
    double temp;
    for(int i=1;i<=(n+1);i++){
        temp = Ab(pivot_row,i);
        Ab(pivot_row,i) = Ab(j,i);
        Ab(j,i) = temp;
    }
}

vec operator/(vec b, matrix2D A){

    int n = A.nrows();

    matrix2D Ab = resize(A, n, n+1); // Concatenated matrix A|b
    for(int i=1; i<=n; i++)
        Ab(i,n+1) = b(i);
    vec x(n);
    double temp;

    //Ab.display();

    for(int j=1; j<n; j++){

        // Find pivot
       int pivot_row = j;
       double pivot_value = fabs(Ab(j,j));
       for(int i=j;i<=n;i++)
          if(fabs(Ab(i,j))>pivot_value){
            pivot_row = i;
            pivot_value = fabs(Ab(i,j));
          }
       if(pivot_row != j)
          swap_rows(Ab,n,pivot_row,j);
       //Ab.display();
       assert(fabs(Ab(j,j))>1.0e-015);

       for(int i=j+1;i<=n;i++){
         temp = Ab(i,j)/Ab(j,j);
         for(int k=j;k<=(n+1);k++)
            Ab(i,k) = Ab(i,k) - temp*Ab(j,k);
       }
           
    }

    x(n) = Ab(n,n+1)/Ab(n,n);
    for(int i=n-1;i>=1;i--){
        temp = Ab(i,n+1);
        for(int j=i+1;j<=n;j++)
            temp -= x(j)*Ab(i,j);
        x(i) = temp/Ab(i,i);
    }
    
    //Ab.display();
    x.roundoffzeros();
    return x;
}


// *********************** Thomas Algorithm ************************

void  tdma( matrix2D & A, vec b, vec & x) // using referenced value of x
{

  int n = b.length();

  /*
      The tridiagonal form of the matrix equation is:

      alpha(i)*x(i-1) + beta(i)*x(i) + gamma(i)*x(i+1) = b(i)

                             for i = 2 to n-1

      BCs : x(i=1) = x1, x(i=n) = xn
 
  */

  vec alpha(n), beta(n), gamma(n);

  // Filling alpha, beta, gamma vectors
  
  for(int i=2;i<n; i++)
  {
    alpha(i) = A(i,i-1);
    beta(i)  = A(i,i);
    gamma(i) = A(i,i+1);
  }

  alpha(1) = 0.0; beta(1) = A(1,1); gamma(1) = A(1,2);
  alpha(n) = A(n,n-1); beta(n) = A(n,n);


   // Forward elimination

  double temp;

  for(int i=3; i<n; i++)
  {
    temp = alpha(i)/beta(i-1);
    beta(i) = beta(i) - temp*gamma(i-1);
    b(i) = b(i) - temp*b(i-1);
  }

  // Backward substitution
   
  x(n-1) = b(n-1)/beta(n-1);

  for(int i=n-2; i>1; i--)
  {
    x(i) = (b(i) - gamma(i)*x(i+1))/beta(i);
  }
}


// ********************* End Thomas Algorithm ***********************


// **********************  GMRES algorithm ************************

void GMRES(matrix2D A, vec b, vec & x0, double tol) // using reference
{
  // USING GIVENS ROTATION to reduce Heissenberg matrix

  int n = b.length();
  //vec x(n);  // Solution vector
  double residual = 1.0;

  // Initial residual
  vec r0(n);
  r0 = b - A*x0;

  // Normalize residual vector
  vec v = r0.normalize();

  // Initialize variables (with least storage space; to be resized)
  // H = Heissenberg matrix

  int k=1;
  matrix2D J = eye(1), Jtotal=eye(2), H(1,1), Htemp(1,1), HH(1,1); 
  matrix2D bb(1,1), c(1,1), cc(1,1), tempMat(1,1), V(n,1), Vold(1,1), hNewCol(1,1);
  vec w(n), vj(n);

  bb(1,1) = r0.norm();

  // Initialise matrix V (matrix of orthogonal basis vectors)
  for(int i=1; i<=n; i++)
     V(i,1) = v(i);  // Filling first column with first basis vector

  // Computation

  while(residual>tol){

     H = resize(H,k+1,k);

     // Arnoldi iteration using Gram-Schmidt process
     w = A*v;

     for(int j=1; j<=k; j++){
        for(int i=1; i<= V.nrows(); i++){
            // setting vector vj to be jth column of V
            vj(i) = V(i,j);
        }

        // calculate inner product
        H(j,k) = dot(vj,w);
        w = w - vj*H(j,k);
     }

     // Gram-Schmidt orthogonalization step

     H(k+1,k) = w.norm();
     v = w.normalize();

     // Appending an additional column to matrix V
     V = resize(V,V.nrows(),k+1);

     for(int i=1;i<=V.nrows();i++){
        // copying centries of v to new column of V
        V(i,k+1) = v(i);
     }

     //////// Least squares step ////////////

     if(k==1){
        // First pass through, Htemp = H
        Htemp = H;
     }
     else{
        // For subsequent passes: Htemp = Jtotal*H
        Jtotal = resize(Jtotal,k+1,k+1);
        Jtotal(k+1,k+1) = 1.0;
        Htemp = Jtotal*H;
        // form next Givens rotation matrix

        J = eye(k-1);
     }

     
     J = resize(J,k+1,k+1);

     // Reducing H-matrix from (k+1,k) order to (k,k)
     J(k,k) = Htemp(k,k)/pow(pow(Htemp(k,k),2)+pow(Htemp(k+1,k),2),0.5);
     J(k,k+1) = Htemp(k+1,k)/pow(pow(Htemp(k,k),2)+pow(Htemp(k+1,k),2),0.5);
     J(k+1,k) = -Htemp(k+1,k)/pow(pow(Htemp(k,k),2)+pow(Htemp(k+1,k),2),0.5);
     J(k+1,k+1) = Htemp(k,k)/pow(pow(Htemp(k,k),2)+pow(Htemp(k+1,k),2),0.5); 

     // Combine together with previous Givens rotations

     Jtotal = J*Jtotal;
     HH = Jtotal*H;

     bb = resize(bb,k+1,1);
     c = Jtotal*bb;

     residual = fabs(c(k+1,1));
     k++; // Increment in k
  }

  cout<<"GMRES procedure converged in "<<k-1<<" iterations"<<endl;
  
  // Extracting upper triangular square matrix
  HH = resize(HH, HH.nrows()-1, HH.ncols());
  cc = resize(c, HH.nrows(),1);

  vec temp(HH.nrows());
  for(int i=1;i<=HH.nrows();i++)
    temp(i) = cc(i,1);

  // Solve linear system
  vec y = temp/HH;

  // Remove the newest column of matrix V
  V = resize(V,V.nrows(),V.ncols()-1);
  //y.display();
  //V.display();
  //x0.display();
  x0 = x0 + V*y;
}

// ********************** END GMRES algorithm ************************

// ********************* Start BiCG algorithm ************************

void BICG(matrix2D A, vec b, vec & x0, double tol) // using reference
{

  


}

// *********************** End BiCG algorithm ************************



// *******************  SPATIAL DERIVATIVE OPERATOR FUNCTIONs  *******

// CD2 scheme for first order spatial derivative
vec first_CD2(vec x, double step_size) {
   int n = x.length();
   vec CD2_x(n);
   for(int i=2; i<=n-1; i++)
     CD2_x(i) = (x(i+1) - x(i-1))/(2.0*step_size);

   return CD2_x;
}

// CD2 scheme for second order spatial derivative
vec second_CD2(vec x, double step_size){
   int n = x.length();
   vec CD2_x(n);
   for(int i=2; i<=n-1; i++)
     CD2_x(i) = (x(i+1) - 2*x(i) + x(i-1))/(step_size*step_size);

   return CD2_x;
}
