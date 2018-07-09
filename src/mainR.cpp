#include <Rcpp.h>

int main (const int argc, const char* argv[]);

void exit(int status) {throw(Rcpp::exception("Exiting"));}

using namespace Rcpp;

// [[Rcpp::export]]
void mainR(const StringVector & arguments)
{
  const int argc = arguments.size();
  const char *argv[argc];
  for (int i=0;i<argc;i++) argv[i]=arguments(i).begin();
  try {
  main(argc,argv);
  } catch(Rcpp::exception &ex) {}
}
